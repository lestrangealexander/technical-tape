#!/usr/bin/env python3
"""A virtual 3-head tape machine, for testing the capture pipeline and the plugin before the real
captures exist.

Record EQ (pre-emphasis + HF record boost) -> 8x-oversampled memory nonlinearity (tanh with low-passed
feedback, even-order asymmetry, dynamic self-erasure) -> playback (de-emphasis, gap loss, head bump,
LF roll-off) -> transport (record-to-repro transit with wow, flutter and drift) -> hiss.

  virtual_tape.py <outdir>        writes the 12 captures the real machine will give:
                                  {ATR,Maxell}_{3.75,7.5,15}ips_{0dB,+6dB}.wav
"""
from __future__ import annotations

import ctypes
import math
import os
import subprocess
import sys

import numpy as np
import scipy.signal as ss

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'pipeline'))
import tt_common as tc  # noqa: E402

FS = tc.FS
OS = 8

SPEEDS = {   # transit s, pre-emphasis us, record boost Hz, gap loss Hz, head bump Hz, hiss dB re 0 VU, flutter scale
    '15':   dict(T0=0.090, tau=35e-6, boost=16000, gap=20500, bump=95.0, hiss=-64.0, flut=1.0),
    '7.5':  dict(T0=0.180, tau=50e-6, boost=12000, gap=15500, bump=55.0, hiss=-61.0, flut=1.6),
    '3.75': dict(T0=0.360, tau=90e-6, boost=7000, gap=9500, bump=32.0, hiss=-57.0, flut=2.6),
}
TAPES = {    # saturation reference (VU), feedback, asymmetry, self-erasure onset (VU of HF), brightness
    'ATR':    dict(sat=12.0, kappa=0.25, a2=0.020, erase=18.0, bright=1.0),
    'Maxell': dict(sat=8.0, kappa=0.35, a2=0.040, erase=12.0, bright=0.85),
    # test-only: a memory-free machine (record EQ -> static saturation -> playback), for checking the method
    'Static': dict(sat=10.0, kappa=0.0, a2=0.030, erase=200.0, bright=1.0),
}
IFACE_LATENCY = 0.006


def _lib():
    so = os.path.join(HERE, 'vt_core.so')
    src = os.path.join(HERE, 'vt_core.c')
    if not os.path.exists(so) or os.path.getmtime(so) < os.path.getmtime(src):
        subprocess.check_call(['cc', '-O2', '-shared', '-fPIC', '-o', so, src, '-lm'])
    lib = ctypes.CDLL(so)
    d = ctypes.c_double
    p = ctypes.POINTER(ctypes.c_double)
    lib.vt_process.argtypes = [p, p, p, ctypes.c_long, d, d, d, d, d, d, d, d]
    return lib


def _ptr(a):
    return a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


def shelf(tz, tp, fs):
    return ss.bilinear([tz, 1.0], [tp, 1.0], fs)


def peaking(f0, gain_db, q, fs):
    A = 10 ** (gain_db / 40)
    w = 2 * np.pi * f0 / fs
    al = np.sin(w) / (2 * q)
    b = np.array([1 + al * A, -2 * np.cos(w), 1 - al * A])
    a = np.array([1 + al / A, -2 * np.cos(w), 1 - al / A])
    return b / a[0], a / a[0]


class Machine:
    def __init__(self, tape='ATR', speed='15'):
        self.tape, self.speed = tape, speed
        self.sp, self.tp = SPEEDS[speed], TAPES[tape]
        self.lib = _lib()
        sp = self.sp
        self.pre = shelf(sp['tau'], sp['tau'] / 5.0, FS)
        self.deemph = shelf(sp['tau'] / 5.0, sp['tau'], FS)
        self.rboost = peaking(sp['boost'], 4.0, 1.1, FS)
        self.bump = peaking(sp['bump'], 2.5, 1.2, FS)
        self.dip = peaking(2.2 * sp['bump'], -1.0, 1.5, FS)
        self.lf = ss.butter(2, max(18.0, sp['bump'] / 3.2), 'highpass', fs=FS)
        g = sp['gap'] * self.tp['bright']
        self.gapfir = ss.firwin2(63, [0, 0.55 * g, 0.85 * g, g, 24000],
                                 [1, 0.96, 0.75, 0.5, 0.0], fs=FS)
        self.gain = 1.0
        self.gain = self._unity_gain()

    # ---------------------------------------------------------------- the deterministic path
    def record_play(self, x):
        """Machine output with no transport effects and no hiss (sample-aligned with x)."""
        sp, tp = self.sp, self.tp
        r = ss.lfilter(*self.pre, x)
        r = ss.lfilter(*self.rboost, r)
        u = ss.resample_poly(r, OS, 1) / tc.vu_amp(tp['sat'])
        uhf = ss.sosfilt(ss.butter(2, 6000, 'highpass', fs=FS * OS, output='sos'), u)
        m = np.empty_like(u)
        fs_os = FS * OS
        self.lib.vt_process(_ptr(u), _ptr(uhf), _ptr(m), len(u),
                            1.0, tp['kappa'], 1.0 - math.exp(-2 * math.pi * 2500.0 / fs_os), tp['a2'],
                            60000.0, 10 ** ((tp['erase'] - tp['sat']) / 20.0), math.exp(-1.0 / (0.004 * fs_os)), fs_os)
        y = ss.resample_poly(m, 1, OS) * tc.vu_amp(tp['sat'])
        y = ss.lfilter(*self.deemph, y)
        y = np.convolve(y, self.gapfir)[31:31 + len(y)]
        for b, a in (self.bump, self.dip, self.lf):
            y = ss.lfilter(b, a, y)
        return y * self.gain

    def _unity_gain(self):
        t = np.arange(int(1.0 * FS)) / FS
        x = tc.vu_amp(0.0) * np.sin(2 * np.pi * 1000.0 * t)
        y = self.record_play(x)[FS // 2:]
        return (tc.vu_amp(0.0) / math.sqrt(2)) / math.sqrt(np.mean(y ** 2))

    # ---------------------------------------------------------------- transport and hiss
    def speed_deviation(self, dur, rng):
        fe = 1000
        n = int(dur * fe) + 1
        t = np.arange(n) / fe
        k = self.sp['flut']
        b, a = ss.butter(1, 0.08, fs=fe)
        d = ss.filtfilt(b, a, rng.standard_normal(n))
        d *= 0.6e-4 * k / d.std()
        b, a = ss.butter(2, [3, 80], 'bandpass', fs=fe)
        fl = ss.filtfilt(b, a, rng.standard_normal(n))
        fl *= 0.8e-4 * k / fl.std()
        e = d + fl
        cap = {'15': 9.6, '7.5': 4.8, '3.75': 2.4}[self.speed]          # capstan rotation
        for f, amp in ((0.62, 2.5e-4), (1.35, 1.2e-4), (cap, 1.6e-4), (19.0, 0.8e-4)):
            e += k * amp * np.sin(2 * np.pi * f * t + rng.uniform(0, 2 * np.pi))
        E = np.concatenate([[0.0], np.cumsum(e[:-1]) / fe])
        return t, e, E

    def capture(self, x, seed, pre_roll=0.5, tail=1.0):
        """What the DAW records: repro output with transit delay, wow/flutter and hiss."""
        rng = np.random.default_rng(seed)
        y = self.record_play(x)
        T0 = self.sp['T0']
        n_out = len(x) + int((T0 + IFACE_LATENCY + tail) * FS)
        dur = n_out / FS + 1.0
        te, e, E = self.speed_deviation(dur, rng)
        t_out = np.arange(n_out) / FS
        Et = np.interp(t_out, te, E)
        tau = t_out - T0 - IFACE_LATENCY
        for _ in range(4):
            tau = t_out + Et - T0 - IFACE_LATENCY - np.interp(tau, te, E)
        out = tc.frac_read(np.concatenate([y, np.zeros(FS)]), tau * FS)
        out[tau < 0] = 0.0
        hiss = rng.standard_normal(n_out) * tc.vu_amp(0.0) / math.sqrt(2) * 10 ** (self.sp['hiss'] / 20)
        hiss = ss.lfilter(*self.deemph, hiss)
        return out + hiss, dict(e=e, te=te, tau=tau)


# The +6 dB pass can't be made with digital gain (the +18 VU sweep would peak at +3 dBFS), so it is
# made in the analog domain: about +6 dB more into the machine and about 6 dB less back into the
# interface. The pipeline measures both offsets; these deliberately aren't exactly 6.
HOT_IN_DB, HOT_OUT_DB = 5.8, -6.3


def main(outdir, only=None):
    os.makedirs(outdir, exist_ok=True)
    x = tc.test_signal().x
    seed = 1
    for tape in TAPES:
        if tape == 'Static' and only is None:
            continue
        for speed in SPEEDS:
            m = Machine(tape, speed)
            for gin, gout, tag in ((0.0, 0.0, '0dB'), (HOT_IN_DB, HOT_OUT_DB, '+6dB')):
                seed += 1
                path = os.path.join(outdir, f'{tape}_{speed}ips_{tag}.wav')
                if only and only not in os.path.basename(path):
                    continue
                y, _ = m.capture(x * 10 ** (gin / 20), seed)
                y *= 10 ** (gout / 20)
                peak = np.max(np.abs(y))
                tc.write_wav24(path, y)
                print(f'{os.path.basename(path)}: {len(y)/FS:.1f} s, peak {20*np.log10(peak):.1f} dBFS', flush=True)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
