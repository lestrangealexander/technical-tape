#!/usr/bin/env python3
"""How closely does a model reproduce the machine it was captured from? Real captures: the reference is
the capture itself.

  validate_capture.py <model.ttm> <normal.wav> [<hot.wav>] [--png out.png] [--renderer path/to/tt_render]

The capture's steady tones are played through the model at the level they went onto the tape, and the
two come out side by side:
  - the 1 kHz tone (0 VU in the normal pass; in the hot pass, 0 VU plus the measured input offset)
  - the 3.15 kHz tone (-10 VU, and -10 VU plus the offset)
For each: the output level, the harmonics H2..H9 (dB re the fundamental), and how well one period of the
waveform matches (the residual of the model's period against the capture's, the fundamental's phase
lined up, every harmonic below 20 kHz included). Wow and flutter are taken out by analysing short
stretches, each lined up on its own fundamental, and averaging. The model renders without flutter.
"""
import json
import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'Tools', 'pipeline'))
import analyse  # noqa: E402
import refengine  # noqa: E402
import tt_common as tc  # noqa: E402
import ttm  # noqa: E402

FS = tc.FS
WIN = 0.05            # seconds per analysis stretch
F_MAX = 20000.0       # harmonics compared up to here


def locate(path):
    """The capture at 48 kHz, on the test file's clock, and the offset of the test file in it."""
    fs, y = tc.read_wav(path, None)
    if fs != FS:
        y = tc.to_48k(y, fs)
    test = tc.test_signal()
    warn = []
    bulk, ratio, _ = analyse.align(y, test, warn)
    if abs(ratio - 1.0) > 20e-6:
        n = int((len(y) - 1) / ratio)
        y = tc.frac_read(y, np.arange(n) * ratio)
        bulk, _, _ = analyse.align(y, test, warn)
    return y, float(bulk)


def peak_freq(x, f_guess):
    win = np.blackman(len(x))
    spec = np.abs(np.fft.rfft(x * win, 8 * len(x)))
    f = np.fft.rfftfreq(8 * len(x), 1 / FS)
    k = int(np.argmax(spec * (np.abs(f - f_guess) < 0.02 * f_guess)))
    a, b, c = np.log(spec[k - 1:k + 2] + 1e-30)
    return float(f[k] + 0.5 * (a - c) / (a - 2 * b + c) * (f[1] - f[0]))


def harmonic_shape(x, f0):
    """Complex amplitudes of harmonics 1..K (K f0 < F_MAX), averaged over short stretches, each stretch's
    harmonics rotated so its fundamental has zero phase (a delay, like flutter, then cancels)."""
    K = int(F_MAX // f0)
    n = int(WIN * FS)
    t = np.arange(n) / FS
    basis = np.concatenate([np.cos(2 * np.pi * f0 * np.outer(t, np.arange(1, K + 1))),
                            np.sin(2 * np.pi * f0 * np.outer(t, np.arange(1, K + 1))), np.ones((n, 1))], axis=1)
    pinv = np.linalg.pinv(basis)
    acc = np.zeros(K, complex)
    count = 0
    for s in range(0, len(x) - n, n):
        co = pinv @ x[s:s + n]
        A = co[:K] - 1j * co[K:2 * K]
        A *= np.exp(-1j * np.arange(1, K + 1) * np.angle(A[0]))
        acc += A
        count += 1
    return acc / max(count, 1)


def period(A, npts=400):
    th = np.linspace(0, 2 * np.pi, npts, endpoint=False)
    return np.real(np.exp(1j * np.outer(th, np.arange(1, len(A) + 1))) @ A)


def render(model, path, x, drive, renderer):
    if renderer is None:
        return refengine.process(model, x, FS, drive_db=drive, factor=8), refengine.latency(model)
    import scipy.io.wavfile as wf
    with tempfile.TemporaryDirectory() as d:
        a, b = os.path.join(d, 'in.wav'), os.path.join(d, 'out.wav')
        wf.write(a, FS, x.astype(np.float32))
        out = subprocess.run([renderer, path, a, b, '--os', '8', '--drive', str(drive), '--no-flutter', '--no-autogain'],
                             capture_output=True, text=True, check=True)
        lat = int([ln for ln in out.stdout.splitlines() if ln.startswith('latency')][0].split()[1])
        _, y = tc.read_wav(b)
        return y, lat


def main():
    args = [a for a in sys.argv[1:]]
    png = renderer = None
    if '--png' in args:
        i = args.index('--png'); png = args[i + 1]; del args[i:i + 2]
    if '--renderer' in args:
        i = args.index('--renderer'); renderer = args[i + 1]; del args[i:i + 2]
    mpath, normal = args[0], args[1]
    hot = args[2] if len(args) > 2 else None
    model = ttm.read(mpath)
    meta = json.loads(model.meta_json) if hasattr(model, 'meta_json') else model.meta
    passes = meta.get('passes', {})
    g0 = 10 ** (passes['cal_gain_db'] / 20.0)
    test = tc.test_signal()
    cal, flt = test.of('cal')[0], test.of('flutter')[0]

    cases = []
    for path, which in ((normal, 'normal'), (hot, 'hot')):
        if path is None:
            continue
        if which == 'hot' and 'hot_in_db' not in passes:
            continue
        y, bulk = locate(path)
        drive = passes['hot_in_db'] if which == 'hot' else 0.0
        gain = g0 * (10 ** (-passes['hot_out_db'] / 20.0) if which == 'hot' else 1.0)
        for sec, fnom, start_s, dur_s in ((cal, 1000.0, 6.0, 4.0), (flt, 3150.0, 6.0, 4.0)):
            a = int(sec.start + start_s * FS)
            cap = tc.frac_read(y, a + bulk + np.arange(int(dur_s * FS))) * gain
            f0 = peak_freq(cap, fnom)
            vu = sec.level + drive
            amp = tc.vu_amp(sec.level)
            n = int((dur_s + 1.0) * FS)
            x = amp * np.sin(2 * np.pi * f0 * np.arange(n) / FS)
            ym, lat = render(model, mpath, x, drive, renderer)
            mod = ym[lat + int(0.5 * FS): lat + int(0.5 * FS) + int(dur_s * FS)]
            Ac, Am = harmonic_shape(cap, f0), harmonic_shape(mod, f0)
            K = min(len(Ac), len(Am))
            Ac, Am = Ac[:K], Am[:K]
            lvl_c, lvl_m = 20 * np.log10(abs(Ac[0])), 20 * np.log10(abs(Am[0]))
            hc = 20 * np.log10(np.abs(Ac[1:9]) / abs(Ac[0]) + 1e-12)
            hm = 20 * np.log10(np.abs(Am[1:9]) / abs(Am[0]) + 1e-12)
            # waveform match: the model's period scaled to the capture's fundamental
            Am_n = Am * abs(Ac[0]) / abs(Am[0])
            resid = 10 * np.log10(np.sum(np.abs(Am_n - Ac) ** 2) / np.sum(np.abs(Ac) ** 2))
            thd_c = 10 * np.log10(np.sum(np.abs(Ac[1:]) ** 2) / abs(Ac[0]) ** 2)
            thd_m = 10 * np.log10(np.sum(np.abs(Am[1:]) ** 2) / abs(Am[0]) ** 2)
            cases.append(dict(label=f'{fnom / 1000:g} kHz at {vu:+.1f} VU ({which} pass)', vu=vu, f0=f0,
                              level_db=lvl_m - lvl_c, thd_capture=thd_c, thd_model=thd_m, h_capture=hc, h_model=hm,
                              waveform_db=resid, Pc=period(Ac), Pm=period(Am_n)))
    print(f'{os.path.basename(mpath)} against {os.path.basename(normal)}' + (f' + {os.path.basename(hot)}' if hot else ''))
    for c in cases:
        print(f'  {c["label"]:30s} level {c["level_db"]:+.2f} dB   THD capture {c["thd_capture"]:6.1f} dB, model '
              f'{c["thd_model"]:6.1f} dB   waveform match {c["waveform_db"]:6.1f} dB')
        print('      H2..H9 capture ' + ' '.join(f'{v:6.1f}' for v in c['h_capture']))
        print('             model   ' + ' '.join(f'{v:6.1f}' for v in c['h_model']))
    if png:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, len(cases), figsize=(3.4 * len(cases), 3.0), squeeze=False)
        for ax, c in zip(axs[0], cases):
            th = np.linspace(0, 2, len(c['Pc']) * 2, endpoint=False)
            ax.plot(th, np.tile(c['Pc'], 2), color='#3A302C', lw=2.2, label='capture')
            ax.plot(th, np.tile(c['Pm'], 2), color='#EE7B5F', lw=1.2, label='model')
            ax.set_title(f'{c["label"]}\nmatch {c["waveform_db"]:.1f} dB', fontsize=8)
            ax.set_xticks([]); ax.set_yticks([])
        axs[0][0].legend(fontsize=7, loc='lower left')
        fig.tight_layout()
        fig.savefig(png, dpi=110)
    return cases


if __name__ == '__main__':
    main()
