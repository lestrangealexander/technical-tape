"""The C++ engine (tt_render) against the Python reference engine (refengine.process), sample by sample.

  compare_engine.py <model.ttm> <tt_render> [seconds]

Both run the same model on the same test signal (sines at several levels and frequencies, a pink-noise
ramp from -20 to +18 VU, drum-like transients). For every oversampling factor and two drive settings
the C++ output, aligned by the reported latencies, is subtracted from the reference output. The
residual should sit near the float32 noise floor (below -90 dB re the output).
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'Tools', 'pipeline'))
import tt_common as tc  # noqa: E402
import refengine  # noqa: E402
import ttm  # noqa: E402

FS = 48000


def test_signal(seconds=6.0, seed=1):
    rng = np.random.default_rng(seed)
    n = int(seconds * FS)
    t = np.arange(n) / FS
    ref = tc.vu_amp(0.0)                          # 0 VU sine peak, full-scale units
    x = np.zeros(n)
    seg = n // 4
    # 1: stepped sines (-20, -6, 0, +6, +12 VU) at 60 Hz, 1 kHz, 7 kHz
    k = seg // 5
    for i, vu in enumerate((-20, -6, 0, 6, 12)):
        sl = slice(i * k, (i + 1) * k)
        f = (60.0, 1000.0, 7000.0)[i % 3]
        x[sl] = tc.vu_amp(vu) * np.sin(2 * np.pi * f * t[sl])
    # 2: pink noise ramp
    w = rng.standard_normal(seg)
    W = np.fft.rfft(w)
    fr = np.fft.rfftfreq(seg, 1 / FS)
    W[1:] /= np.sqrt(fr[1:])
    W[0] = 0
    p = np.fft.irfft(W, seg)
    p /= np.sqrt(np.mean(p ** 2))
    ramp = 10 ** (np.linspace(-20, 18, seg) / 20)
    x[seg:2 * seg] = p * ramp * ref / np.sqrt(2)
    # 3: transients (decaying noise bursts and clicks)
    for j in range(12):
        s0 = 2 * seg + j * (seg // 12)
        L = seg // 12
        env = np.exp(-np.arange(L) / (0.02 * FS))
        x[s0:s0 + L] = rng.standard_normal(L) * env * ref * 10 ** ((j % 4 * 6 - 6) / 20)
    # 4: two-tone at +6 VU
    sl = slice(3 * seg, n)
    x[sl] = tc.vu_amp(6) * (0.8 * np.sin(2 * np.pi * 60 * t[sl]) + 0.2 * np.sin(2 * np.pi * 7000 * t[sl]))
    return x


def write_f32(path, x):
    import scipy.io.wavfile as wf
    wf.write(path, FS, np.asarray(x, np.float32))


def run_cpp(binary, model_path, x, factor, drive):
    with tempfile.TemporaryDirectory() as d:
        a, b = os.path.join(d, 'in.wav'), os.path.join(d, 'out.wav')
        write_f32(a, x)
        out = subprocess.run([binary, model_path, a, b, '--os', str(factor), '--drive', str(drive), '--no-flutter',
                              '--no-autogain', '--random-blocks', '--block', '700'],
                             capture_output=True, text=True, check=True)
        lat = int([l for l in out.stdout.splitlines() if l.startswith('latency')][0].split()[1])
        rt = [l for l in out.stdout.splitlines() if l.startswith('realtime')]
        import scipy.io.wavfile as wf
        _, y = wf.read(b)
        return y.astype(float), lat, (rt[0] if rt else '')


def main():
    model_path, binary = sys.argv[1], sys.argv[2]
    seconds = float(sys.argv[3]) if len(sys.argv) > 3 else 6.0
    model = ttm.read(model_path)
    x = test_signal(seconds)
    worst = -999.0
    for drive in (0.0, 12.0):
        for factor in (1, 2, 4, 8, 16):
            yr = refengine.process(model, x, FS, drive_db=drive, factor=factor)
            yc, lat_c, rt = run_cpp(binary, model_path, x, factor, drive)
            shift = lat_c - refengine.latency(model)
            yc = yc[shift:]
            n = min(len(yc), len(yr))
            e = yc[:n] - yr[:n]
            sig = np.sqrt(np.mean(yr[:n] ** 2))
            res = 20 * np.log10(np.sqrt(np.mean(e ** 2)) / sig + 1e-30)
            pk = 20 * np.log10(np.max(np.abs(e)) / np.max(np.abs(yr[:n])) + 1e-30)
            worst = max(worst, res)
            print(f'drive {drive:+5.1f} dB  OS {factor:2d}x  residual rms {res:7.1f} dB  peak {pk:7.1f} dB   {rt}')
    print(f'worst residual {worst:.1f} dB re output')
    return 0 if worst < -80 else 1


if __name__ == '__main__':
    sys.exit(main())
