"""Shared pieces of the Technical Tape capture pipeline: the test signal's exact layout, WAV I/O,
the synchronized sweep, band-limited fractional reads and a few numeric helpers.

The capture file is regenerated sample-for-sample from Tools/capture/make_capture.py, so every
section (calibration tone, sync markers, the 15 sweeps, pink noise, flutter tone) is known exactly.
"""
from __future__ import annotations

import math
import os
import sys
from dataclasses import dataclass, field

import numpy as np
import scipy.io.wavfile as wavfile
import scipy.signal as ss

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'capture'))
import make_capture as mc  # noqa: E402

FS = mc.FS                                   # 48 kHz
REF_PEAK = mc.REF_PEAK                       # 0 VU sine peak (-15 dBFS)
F1, F2 = 10.0, 22000.0
SWEEP_L, SWEEP_SEC = mc.sync_sweep_params(F1, F2, 12.0)
N_SWEEP = int(round(SWEEP_SEC * FS))
SWEEP_LEVELS = list(range(-30, 13, 3))       # VU, as recorded in the normal pass
NOISE_LEVELS = (-20, -10, 0, 3)              # VU RMS
NOISE_SEEDS = (100, 101, 102, 103)


def vu_amp(vu):
    """Sine peak amplitude (full scale = 1) for a level in VU."""
    return REF_PEAK * 10.0 ** (np.asarray(vu, float) / 20.0)


@dataclass
class Section:
    kind: str           # silence, cal, sync, sweep, noise, flutter
    start: int          # sample index in the test file
    length: int
    level: float | None = None
    index: int = 0      # sweep / noise number


@dataclass
class TestSignal:
    x: np.ndarray
    sections: list = field(default_factory=list)

    def of(self, kind):
        return [s for s in self.sections if s.kind == kind]


_TEST = None


def test_signal() -> TestSignal:
    """The capture test file, rebuilt exactly as make_capture.main() lays it out."""
    global _TEST
    if _TEST is not None:
        return _TEST
    parts, secs = [], []
    pos = 0

    def add(x, kind=None, level=None, index=0):
        nonlocal pos
        if kind:
            secs.append(Section(kind, pos, len(x), level, index))
        parts.append(x)
        pos += len(x)

    add(mc.silence(2.0), 'silence')
    add(mc.tone(1000, 30.0, 0.0), 'cal', 0.0)
    add(mc.silence(2.0))
    add(mc.marker_burst(), 'sync', -6.0, 0)
    add(mc.silence(2.0))
    for i, vu in enumerate(SWEEP_LEVELS):
        add(mc.sweep(F1, F2, 12.0, vu), 'sweep', float(vu), i)
        add(mc.silence(3.0))
    for i, vu in enumerate(NOISE_LEVELS):
        add(mc.pink(8.0, vu, NOISE_SEEDS[i]), 'noise', float(vu), i)
        add(mc.silence(2.0))
    add(mc.tone(3150, 20.0, -10.0), 'flutter', -10.0)
    add(mc.silence(2.0))
    add(mc.marker_burst(), 'sync', -6.0, 1)
    add(mc.silence(3.0), 'silence')
    _TEST = TestSignal(np.concatenate(parts), secs)
    return _TEST


def sweep_shape():
    """Unit-amplitude synchronized sweep with make_capture's fades (sample-exact)."""
    t = np.arange(N_SWEEP) / FS
    x = np.sin(2 * np.pi * F1 * SWEEP_L * np.exp(t / SWEEP_L))
    return mc.fade(x, 0.10, 0.005)


def sweep_phase():
    t = np.arange(N_SWEEP) / FS
    return 2 * np.pi * F1 * SWEEP_L * np.exp(t / SWEEP_L), F1 * np.exp(t / SWEEP_L)


# ------------------------------------------------------------------ WAV I/O
def read_wav(path, channel=None):
    """Any PCM / float WAV -> float64 mono in [-1, 1). A stereo file uses the louder channel
    unless `channel` is given."""
    fs, d = wavfile.read(path)
    if d.dtype == np.int16:
        d = d / 32768.0
    elif d.dtype == np.int32:
        d = d / 2147483648.0
    elif d.dtype == np.uint8:
        d = (d.astype(np.float64) - 128.0) / 128.0
    else:
        d = d.astype(np.float64)
    if d.ndim == 2:
        if channel is None:
            channel = int(np.argmax(np.sqrt(np.mean(d ** 2, axis=0))))
        d = d[:, channel]
    return fs, d


def to_48k(y, fs_in):
    """A capture at another rate onto the test file's 48 kHz clock: linear-phase Kaiser FIR, flat to
    22 kHz, -130 dB wherever something would fold into 0-22 kHz (so the tape's ultrasonic harmonics,
    hiss and any bias leak captured at 96/192 kHz can't alias into the measurement)."""
    fs_in = int(round(fs_in))
    if fs_in == FS:
        return np.asarray(y, float)
    g = math.gcd(fs_in, FS)
    up, down = FS // g, fs_in // g
    rate = fs_in * up                                  # the rate the filter runs at
    if up == 1:
        # 96 / 192 kHz: pass 22 kHz, stop from 26 kHz (what folds lands above 22 kHz)
        n = int(np.ceil(122.0 * rate / (14.36 * 3500.0))) | 1
        h = ss.firwin(2 * n + 1, 24000.0, window=('kaiser', 13.4), fs=rate)
    else:
        # 44.1 / 88.2 / 176.4 kHz: as steep as the input allows, cutoff between 22 kHz and Nyquist
        top = min(fs_in / 2.0, FS / 2.0)
        cut = 0.5 * (22000.0 + top) if top > 22100.0 else 0.995 * top
        width = max(2.0 * (top - cut), 100.0)
        n = int(np.ceil(122.0 * rate / (14.36 * width))) | 1
        h = ss.firwin(min(n, 64 * max(up, down) * 8 + 1), cut, window=('kaiser', 13.4), fs=rate)
    return ss.resample_poly(np.asarray(y, float), up, down, window=h)        # (it applies the x up gain)


def write_wav24(path, x, fs=FS):
    pcm = np.clip(np.round(np.asarray(x) * 8388607.0), -8388608, 8388607).astype('<i4')
    raw = pcm.view(np.uint8).reshape(-1, 4)[:, :3].tobytes()
    import struct
    with open(path, 'wb') as fh:
        fh.write(b'RIFF' + struct.pack('<I', 36 + len(raw)) + b'WAVE')
        fh.write(b'fmt ' + struct.pack('<IHHIIHH', 16, 1, 1, fs, fs * 3, 3, 24))
        fh.write(b'data' + struct.pack('<I', len(raw)) + raw)


# ------------------------------------------------------------------ band-limited fractional reads
_K = np.arange(-31, 33)
_Q = 2048
_mu = np.arange(_Q + 1)[:, None] / _Q
_xx = _K[None, :] - _mu
_TAB = (np.sinc(_xx) * np.i0(9.0 * np.sqrt(np.clip(1 - (_xx / 32.5) ** 2, 0, None))) / np.i0(9.0))


def frac_read(r, pos, chunk=40000):
    """r evaluated at fractional sample positions `pos` (64-tap Kaiser sinc, 2048 phases)."""
    out = np.empty(len(pos))
    for s in range(0, len(pos), chunk):
        p = pos[s:s + chunk]
        n = np.floor(p).astype(np.int64)
        q = np.round((p - n) * _Q).astype(np.int64)
        idx = np.clip(n[:, None] + _K[None, :], 0, len(r) - 1)
        out[s:s + chunk] = np.sum(r[idx] * _TAB[q], axis=1)
    return out


# ------------------------------------------------------------------ helpers
def db(x, floor=1e-30):
    return 20.0 * np.log10(np.maximum(np.abs(x), floor))


def next_pow2(n):
    return 1 << int(math.ceil(math.log2(max(1, n))))


def frac_delay_xcorr(a, b, max_lag):
    """Lag (samples, sub-sample by parabolic peak) that best aligns b to a (b[n+lag] ~ a[n])."""
    n = next_pow2(len(a) + len(b))
    c = np.fft.irfft(np.fft.rfft(b, n) * np.conj(np.fft.rfft(a, n)), n)
    c = np.concatenate([c[-max_lag:], c[:max_lag + 1]])
    i = int(np.argmax(np.abs(c)))
    if 0 < i < len(c) - 1:
        y0, y1, y2 = np.abs(c[i - 1]), np.abs(c[i]), np.abs(c[i + 1])
        off = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2) if (y0 - 2 * y1 + y2) != 0 else 0.0
    else:
        off = 0.0
    return i - max_lag + off, c[i]


def cheb_coef(n, m):
    """Coefficient of cos(m t) in cos^n(t) (m >= 1, same parity as n)."""
    if m > n or (n - m) % 2:
        return 0.0
    return 2.0 ** (1 - n) * math.comb(n, (n - m) // 2)


def kaiser_lowpass(numtaps, cutoff, beta):
    """Linear-phase low-pass, cutoff as a fraction of the sample rate (0..0.5)."""
    return ss.firwin(numtaps, cutoff * 2.0, window=('kaiser', beta))
