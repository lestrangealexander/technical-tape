#!/usr/bin/env python3
"""How close is a model to the machine it was captured from? (Virtual machine only: needs its clean path.)

  validate_model.py <model.ttm> <tape> <speed> [renderer]
  renderer: 'ref' (default, the Python reference engine) or a path to the C++ tt_render binary

Compares steady sines (gain, H2, H3), pink noise spectra, a two-tone IMD test and a transient null.
"""
import os
import subprocess
import sys
import tempfile

import numpy as np
import scipy.signal as ss

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'Tools', 'pipeline'))
sys.path.insert(0, os.path.join(HERE, '..', 'Tools', 'sim'))
import tt_common as tc  # noqa: E402  (puts Tools/capture on the path)
import make_capture as mc  # noqa: E402
import refengine  # noqa: E402
import tt_common as tc  # noqa: E402
import ttm  # noqa: E402
import virtual_tape as vt  # noqa: E402

FS = tc.FS


def render(model, path, x, renderer, factor=8, drive=0.0):
    if renderer == 'ref':
        return refengine.process(model, x, FS, drive_db=drive, factor=factor), refengine.latency(model)
    with tempfile.TemporaryDirectory() as d:
        a, b = os.path.join(d, 'in.wav'), os.path.join(d, 'out.wav')
        tc.write_wav24(a, x * 0.5) if False else ss_write_f32(a, x)
        out = subprocess.run([renderer, path, a, b, '--os', str(factor), '--drive', str(drive), '--no-flutter',
                              '--no-autogain'], capture_output=True, text=True, check=True)
        lat = int([l for l in out.stdout.splitlines() if l.startswith('latency')][0].split()[1])
        fs, y = tc.read_wav(b)
        return y, lat


def ss_write_f32(path, x):
    import scipy.io.wavfile as wf
    wf.write(path, FS, np.asarray(x, np.float32))


def align(y, lat, n):
    return np.concatenate([y[lat:], np.zeros(lat)])[:n]


def tone_levels(y, f, n0):
    seg = y[n0:]
    w = np.blackman(len(seg))
    Y = np.abs(np.fft.rfft(seg * w))
    fr = np.fft.rfftfreq(len(seg), 1 / FS)
    def a(fq):
        if fq >= FS / 2 - 500:
            return np.nan
        k = int(round(fq * len(seg) / FS))
        return np.sqrt(np.sum(Y[k - 3:k + 4] ** 2))
    return a(f), a(2 * f), a(3 * f)


def main():
    path, tape, speed = sys.argv[1:4]
    renderer = sys.argv[4] if len(sys.argv) > 4 else 'ref'
    model = ttm.read(path)
    m = vt.Machine(tape, speed)
    print(f'model {os.path.basename(path)} vs virtual {tape} {speed} ips, renderer {renderer}')

    # ---- delay between model and machine, from broadband noise at low level
    x = mc.pink(2.0, -20.0, 11)
    ym, lat = render(model, path, x, renderer)
    yt = m.record_play(x)
    lag, _ = tc.frac_delay_xcorr(yt, ym, 2000)
    print(f'latency reported {lat}, measured {lag:.2f} samples')
    shift = lag

    def model_out(x, **kw):
        y, _ = render(model, path, x, renderer, **kw)
        return tc.frac_read(np.concatenate([y, np.zeros(4096)]), np.arange(len(x)) + shift)

    # ---- steady sines
    print('\nSteady sines: model minus machine (fundamental dB | H2 dB re fund | H3 dB re fund)')
    freqs = (50, 100, 300, 1000, 3000, 6000, 10000, 15000)
    hdr = ''.join(f'{f:>22}' for f in freqs)
    print(f'{"VU":>5}{hdr}')
    worst = 0.0
    for vu in (-20, -10, 0, 6, 12, 18):
        row = []
        for f in freqs:
            n = int(0.6 * FS)
            t = np.arange(n) / FS
            x = tc.vu_amp(vu) * np.sin(2 * np.pi * f * t)
            a = tone_levels(model_out(x), f, n // 3)
            b = tone_levels(m.record_play(x), f, n // 3)
            g = 20 * np.log10(a[0] / b[0])
            h2 = 20 * np.log10(a[1] / a[0]) - 20 * np.log10(b[1] / b[0])
            h3 = 20 * np.log10(a[2] / a[0]) - 20 * np.log10(b[2] / b[0])
            hb2, hb3 = 20 * np.log10(b[1] / b[0]), 20 * np.log10(b[2] / b[0])
            worst = max(worst, abs(g))
            def fmt(d, ref):
                return '   .  ' if np.isnan(d) or ref < -90 else f'{d:+6.1f}'
            row.append(f'{g:+6.2f}|{fmt(h2, hb2)}|{fmt(h3, hb3)}')
        print(f'{vu:>+5}' + ''.join(f'{r:>22}' for r in row))
    print(f'worst fundamental error {worst:.2f} dB')

    # ---- pink noise spectra
    print('\nPink noise: model minus machine, 1/3-octave bands (dB)')
    cen = 1000 * 2 ** (np.arange(-14, 13) / 3)
    print(f'{"VU":>5}' + ''.join(f'{c:>7.0f}' for c in cen[::2]))
    for vu in (-10, 0, 6, 12):
        x = mc.pink(4.0, vu, 21)
        a, b = model_out(x)[FS // 2:], m.record_play(x)[FS // 2:]
        f, pa = ss.welch(a, FS, nperseg=8192)
        _, pb = ss.welch(b, FS, nperseg=8192)
        d = [10 * np.log10(np.sum(pa[(f >= c * 2 ** -(1 / 6)) & (f < c * 2 ** (1 / 6))])
                           / np.sum(pb[(f >= c * 2 ** -(1 / 6)) & (f < c * 2 ** (1 / 6))])) for c in cen]
        print(f'{vu:>+5}' + ''.join(f'{v:+7.2f}' for v in d[::2]) + f'   rms {np.sqrt(np.mean(np.square(d))):.2f}')

    # ---- two-tone IMD (SMPTE-like 60 Hz + 7 kHz, 4:1)
    print('\nTwo-tone 60 Hz + 7 kHz (4:1): sidebands at 7k +/- 60, 120 Hz, dB re 7k')
    for vu in (0, 6, 12):
        n = FS
        t = np.arange(n) / FS
        amp = tc.vu_amp(vu)
        x = amp * (0.8 * np.sin(2 * np.pi * 60 * t) + 0.2 * np.sin(2 * np.pi * 7000 * t))
        res = []
        for y in (model_out(x), m.record_play(x)):
            seg = y[n // 4:]
            Y = np.abs(np.fft.rfft(seg * np.blackman(len(seg))))
            def at(fq):
                k = int(round(fq * len(seg) / FS))
                return np.sqrt(np.sum(Y[k - 2:k + 3] ** 2))
            c = at(7000)
            res.append([20 * np.log10(at(7000 + d) / c) for d in (-120, -60, 60, 120)])
        print(f'{vu:+d} VU  model ' + ' '.join(f'{v:6.1f}' for v in res[0]) + '   machine ' + ' '.join(f'{v:6.1f}' for v in res[1]))

    # ---- transient null (drum-like hits)
    print('\nTransient null: residual after subtracting the machine output (dB re machine output)')
    rng = np.random.default_rng(5)
    n = 2 * FS
    x = np.zeros(n)
    for k in range(8):
        p = int((0.05 + 0.24 * k) * FS)
        L = int(0.2 * FS)
        tt = np.arange(L) / FS
        f0 = rng.uniform(60, 300)
        hit = np.sin(2 * np.pi * f0 * tt * (1 + 0.5 * np.exp(-tt * 30))) * np.exp(-tt * 18)
        hit += 0.4 * rng.standard_normal(L) * np.exp(-tt * 60)
        x[p:p + L] += hit
    for vu in (0, 6, 12):
        xs = x / np.sqrt(np.mean(x ** 2)) * tc.vu_amp(vu) / np.sqrt(2)
        a, b = model_out(xs), m.record_play(xs)
        r = 10 * np.log10(np.sum((a - b) ** 2) / np.sum(b ** 2))
        print(f'{vu:+d} VU RMS: {r:.1f} dB')


if __name__ == '__main__':
    main()
