"""Engine checks beyond the sample-exact comparison: aliasing, switching without clicks, wow & flutter.

  engine_checks.py <tt_render> <model.ttm> [<second_model.ttm>]
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile

import numpy as np
import scipy.io.wavfile as wf
import scipy.signal as ss

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'Tools', 'pipeline'))
import tt_common as tc  # noqa: E402
import ttm  # noqa: E402

FS = 48000


def render(binary, models, x, *opts):
    with tempfile.TemporaryDirectory() as d:
        a, b = os.path.join(d, 'in.wav'), os.path.join(d, 'out.wav')
        wf.write(a, FS, np.asarray(x, np.float32))
        out = subprocess.run([binary, ','.join(models), a, b, *map(str, opts)], capture_output=True, text=True,
                             check=True)
        lat = int([l for l in out.stdout.splitlines() if l.startswith('latency')][0].split()[1])
        _, y = wf.read(b)
        return y.astype(float), lat


def spur_db(y, f0, skip=0.5):
    y = y[int(skip * FS):]
    n = 1 << int(np.floor(np.log2(len(y))))
    y = y[:n]
    w = ss.windows.blackmanharris(n)
    Y = np.abs(np.fft.rfft(y * w)) / np.sum(w) * 2
    f = np.fft.rfftfreq(n, 1 / FS)
    fund = Y[np.abs(f - f0) < 20].max()
    mask = (np.abs(f - f0) > 60) & (f > 30)
    k = np.argmax(Y * mask)
    return 20 * np.log10(Y[k] / fund + 1e-30), f[k]


def aliasing(binary, model):
    print('Aliasing: 15 kHz sine at +6 VU, INPUT +12 dB (an in-band spur can only be an alias)')
    t = np.arange(int(2.5 * FS)) / FS
    x = tc.vu_amp(6) * np.sin(2 * np.pi * 15000 * t)
    for os_ in (1, 2, 4, 8, 16):
        y, _ = render(binary, [model], x, '--os', os_, '--drive', 12, '--no-flutter', '--no-autogain')
        db, f = spur_db(y, 15000)
        print(f'  OS {os_:2d}x  worst spur {db:7.1f} dB re fundamental at {f:7.0f} Hz')


def clicks(binary, models):
    print('Switching: 1 kHz sine at 0 VU; energy above 14 kHz (10 ms windows), dB re the sine')
    t = np.arange(int(3.0 * FS)) / FS
    x = tc.vu_amp(0) * np.sin(2 * np.pi * 1000 * t)
    sos = ss.butter(10, 14000, 'highpass', fs=FS, output='sos')
    cases = [('OS 8x -> 2x', [models[0]], ('--switch-os', 2, '--at', 1.0)),
             ('OS 8x -> 1x', [models[0]], ('--switch-os', 1, '--at', 1.0)),
             ('bypass on', [models[0]], ('--bypass-at', 1.0))]
    if len(models) > 1:
        cases.append(('model A -> B', models[:2], ('--switch-model', 1, '--at', 1.0)))
    for name, ms, opts in cases:
        y, lat = render(binary, ms, x, '--no-flutter', '--no-autogain', *opts)
        h = ss.sosfilt(sos, y)
        win = int(0.01 * FS)
        r = np.sqrt(np.convolve(h ** 2, np.ones(win) / win, 'valid'))
        ref = tc.vu_amp(0) / np.sqrt(2)
        steady = 20 * np.log10(np.median(r[int(0.4 * FS):int(0.9 * FS)]) / ref + 1e-30)
        around = 20 * np.log10(r[int(0.9 * FS):int(1.6 * FS)].max() / ref + 1e-30)
        print(f'  {name:14s} steady {steady:7.1f} dB, worst around the switch {around:7.1f} dB')


def flutter(binary, model):
    m = ttm.read(model)
    print('Wow & flutter: 3150 Hz tone, W&F on; speed deviation from the instantaneous frequency')
    t = np.arange(int(20 * FS)) / FS
    x = tc.vu_amp(-10) * np.sin(2 * np.pi * 3150 * t)
    for amt in (1.0, 0.0):
        y, _ = render(binary, [model], x, '--wf', amt, '--no-autogain')
        a = ss.hilbert(y[FS:])
        ph = np.unwrap(np.angle(a))
        fi = np.diff(ph) * FS / (2 * np.pi)
        dec = 48
        b, a2 = ss.butter(4, 400, fs=FS)
        fi = ss.filtfilt(b, a2, fi)[::dec]
        fe = FS / dec
        dev = fi / 3150.0 - 1.0
        bu, au = ss.butter(1, [0.5, 30], 'bandpass', fs=fe)
        un = 200 * np.std(ss.filtfilt(bu, au, dev)[int(fe):-int(fe)])
        bw, aw = ss.butter(1, [2.0, 8.0], 'bandpass', fs=fe)
        we = 100 * np.std(ss.filtfilt(bw, aw, dev)[int(fe):-int(fe)])
        print(f'  amount {amt:.0f}: unweighted 2-sigma {un:.3f} %, weighted rms {we:.3f} %')
    meta = m.meta if hasattr(m, 'meta') else {}
    print('  model (measured on the capture):', meta.get('wow_flutter', meta) if isinstance(meta, dict) else meta)


def main():
    binary, models = sys.argv[1], sys.argv[2:]
    aliasing(binary, models[0])
    clicks(binary, models)
    flutter(binary, models[0])


if __name__ == '__main__':
    main()
