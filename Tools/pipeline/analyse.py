"""Per-capture analysis: find the test signal in a recording of the tape machine and measure everything
the model needs from it.

  alignment     sync markers at both ends -> bulk delay and any constant speed ratio
  calibration   1 kHz 0 VU tone -> output level, speed check, THD
  sweeps        de-warped (transport jitter removed), deconvolved, harmonic responses H_m(F) for
                m = 1..M_MEAS per unit input amplitude, with the noise floor of every window
  flutter       jitter estimates from all sweeps and the 3150 Hz tone
  pink noise    the four sections, input and output, for checks and the detector fit
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import scipy.signal as ss

import tt_common as tc

FS = tc.FS
L = tc.SWEEP_L
NK = 32768                      # kernel spectra: 1.46 Hz bins at 48 kHz
NF = 1 << 21                    # deconvolution FFT
PRE = int(4.3 * FS)             # output segment starts this long before each sweep: the noise window
                                # (orders 18..26, around -5 s) must see real output down to 20 Hz
POST = int(1.25 * FS)           # ... and runs this long after it
M_MEAS = 16                     # harmonic orders measured
WIN_PRE = int(0.005 * FS)
NOISE_ORDERS = (18, 26)         # the deconvolution between these orders' slots holds only noise
F_LO_OK, F_HI_OK = 11.0, 21800.0
FK = np.fft.rfftfreq(NK, 1.0 / FS)


def harmonic_post(m):
    if m == 1:
        return int(0.35 * FS)
    gap = L * math.log(m / (m - 1.0))
    return int(min(0.15, 0.8 * gap) * FS)


@dataclass
class SweepResult:
    level: float                 # nominal VU in this file
    H: np.ndarray                # (M_MEAS, NK/2+1) complex, per unit nominal input amplitude
    noise: np.ndarray            # (M_MEAS, NK/2+1) expected noise power in each H_m bin
    jitter: np.ndarray           # de-warp correction actually applied over the sweep, s (48 kHz grid)
    valid: np.ndarray            # bool mask of where `jitter` was measured (40 Hz..21 kHz)


@dataclass
class CaptureResult:
    path: str
    n_samples: int
    channel: int | None
    bulk_delay: float            # samples (fractional)
    speed_ratio: float
    marker_quality: tuple
    cal_dbfs: float              # 1 kHz tone RMS in the capture, dBFS
    cal_freq: float
    cal_thd_db: float
    noise_floor_dbfs: float
    sweeps: list = field(default_factory=list)
    flutter_sweeps: list = field(default_factory=list)   # (t0 s, 1 kHz-rate delay s) per sweep, hi band
    flutter_tone: np.ndarray | None = None               # 1 kHz-rate delay s from the 3150 Hz tone
    pink: list = field(default_factory=list)              # (level VU, input, output)
    warnings: list = field(default_factory=list)


def _find_marker(y, template, lo, hi):
    """Best match of `template` in y[lo:hi] by normalised cross-correlation (level-independent, so the
    loud sweeps and noise nearby can't win), refined to a fraction of a sample."""
    lo, hi = max(0, int(lo)), min(len(y), int(hi))
    seg = y[lo:hi]
    m = len(template)
    n = tc.next_pow2(len(seg) + m)
    c = np.fft.irfft(np.fft.rfft(seg, n) * np.conj(np.fft.rfft(template, n)), n)[:len(seg) - m]
    e = np.concatenate([[0.0], np.cumsum(seg ** 2)])
    local = np.sqrt(np.maximum(e[m:len(seg)] - e[:len(seg) - m], 1e-30))
    ncc = c / (local * np.linalg.norm(template))
    i = int(np.argmax(ncc))
    off = 0.0
    if 0 < i < len(c) - 1:
        y0, y1, y2 = c[i - 1], c[i], c[i + 1]
        den = y0 - 2 * y1 + y2
        if den != 0:
            off = 0.5 * (y0 - y2) / den
    return lo + i + off, float(ncc[i])


MARKER_OFFSETS = (0.0, 0.5, 1.25)          # the three chirps of a sync burst (make_capture.marker_burst)


def _ncc(y, template, lo, hi, sos=None):
    """Normalised cross-correlation of `template` against y[lo:hi] (optionally both band-limited):
    returns (ncc array, start offset)."""
    lo, hi = max(0, int(lo)), min(len(y), int(hi))
    seg = y[lo:hi]
    if sos is not None:
        seg = ss.sosfiltfilt(sos, seg)
        template = ss.sosfiltfilt(sos, template)
    m = len(template)
    if len(seg) <= m:
        return np.zeros(1), lo
    n = tc.next_pow2(len(seg) + m)
    c = np.fft.irfft(np.fft.rfft(seg, n) * np.conj(np.fft.rfft(template, n)), n)[:len(seg) - m]
    e = np.concatenate([[0.0], np.cumsum(seg ** 2)])
    local = np.sqrt(np.maximum(e[m:len(seg)] - e[:len(seg) - m], 1e-30))
    return c / (local * np.linalg.norm(template)), lo


def _find_burst(y, test, burst, lo, hi):
    """Start of a sync burst in y[lo:hi], robust to wow and to a dull top end: each 20 ms chirp is
    matched on its own (band-limited to 150 Hz - 5 kHz), and the candidate whose three chirps sit at
    0, 0.5 and 1.25 s (+/- 15 ms, so wow can't break the pattern) wins. Refined to a fraction of a sample
    on the full-band chirp. Returns (position, match quality 0..1)."""
    n = int(0.020 * FS)
    chirp = test.x[burst.start:burst.start + n]
    sos = ss.butter(4, [150.0, 5000.0], 'bandpass', fs=FS, output='sos')
    ncc, off = _ncc(y, chirp, lo, hi + int(1.3 * FS), sos)
    if len(ncc) < 2:
        return float(lo), 0.0
    tol = int(0.015 * FS)
    # local maxima of the first chirp's match, best first
    cand = np.flatnonzero((ncc[1:-1] > ncc[:-2]) & (ncc[1:-1] >= ncc[2:]) & (ncc[1:-1] > 0.3)) + 1
    cand = cand[np.argsort(ncc[cand])[::-1]][:40]
    best, best_score = None, -1.0
    for c in cand:
        score = ncc[c]
        for o in MARKER_OFFSETS[1:]:
            k = c + int(round(o * FS))
            a, b = max(0, k - tol), min(len(ncc), k + tol + 1)
            score += ncc[a:b].max() if b > a else 0.0
        if score > best_score:
            best, best_score = c, score
    if best is None:
        return float(lo), 0.0
    # refine on the full band, around the band-limited estimate
    c2, off2 = _ncc(y, chirp, off + best - 64, off + best + 64 + n)
    i = int(np.argmax(c2))
    frac = 0.0
    if 0 < i < len(c2) - 1:
        y0, y1, y2 = c2[i - 1], c2[i], c2[i + 1]
        den = y0 - 2 * y1 + y2
        if den != 0:
            frac = 0.5 * (y0 - y2) / den
    return off2 + i + frac, float(best_score / 3.0)


def align(y, test, warnings):
    """Bulk delay (samples) and speed ratio from the two sync bursts."""
    s0, s1 = test.of('sync')
    p0, q0 = _find_burst(y, test, s0, 0, min(len(y), s0.start + 60 * FS))
    exp1 = p0 + (s1.start - s0.start)             # where the end burst would be at the right speed
    p1, q1 = _find_burst(y, test, s1, exp1 - 12 * FS, exp1 + 12 * FS)
    if min(q0, q1) < 0.5:
        warnings.append(f'Sync markers are weak (match {q0:.2f} / {q1:.2f}); check that this is the right file '
                        'and that it was recorded off the repro head.')
    ratio = (p1 - p0) / float(s1.start - s0.start)
    return p0 - s0.start, ratio, (q0, q1)


def _demod_phase(z, sos):
    zl = ss.sosfiltfilt(sos, z.real) + 1j * ss.sosfiltfilt(sos, z.imag)
    return np.unwrap(np.angle(zl)), np.abs(zl)


def _smooth(v, width):
    if width < 3:
        return v
    w = np.hanning(width)
    w /= w.sum()
    pad = width // 2
    vp = np.concatenate([np.full(pad, v[0]), v, np.full(pad, v[-1])])
    return ss.fftconvolve(vp, w, mode='same')[pad:pad + len(v)]


def _fill_invalid(d, ok):
    """Hold the nearest valid value outside the valid stretch of each row."""
    n = d.shape[1]
    idx = np.arange(n)
    for k in range(d.shape[0]):
        v = np.flatnonzero(ok[k])
        if len(v) == 0:
            d[k] = 0.0
            continue
        fwd = np.maximum.accumulate(np.where(ok[k], idx, -1))
        fwd[:v[0]] = v[0]
        d[k] = d[k][fwd]
    return d


def dewarp_corrections(segs, levels):
    """Transport jitter of every sweep of one pass, from the phase of its fundamental.

    The machine's own phase response (common to all sweeps, or varying smoothly with level: a
    weighted quadratic in level, smoothed along the sweep) is estimated across the sweeps and kept;
    what is left is jitter, converted to time. Stretches where a sweep's fundamental is more than
    40 dB below its own maximum (deep in the noise) are held at the nearest good value."""
    theta, f_t = tc.sweep_phase()
    bands = {'hi': (ss.butter(4, 60, fs=FS, output='sos'), (f_t >= 250) & (f_t <= 21000)),
             'lo': (ss.butter(4, 12, fs=FS, output='sos'), (f_t >= 40) & (f_t < 400))}
    lev = np.asarray(levels, float)
    lc = (lev - lev.mean()) / 10.0
    w_mean = int(L * math.log(2) / 6 * FS) | 1                           # 1/6 octave along the sweep
    w_lev = int(L * math.log(2) / 3 * FS) | 1                            # 1/3 octave
    rot = np.exp(-1j * theta)
    delta, oks = {}, {}
    for name, (sos, band) in bands.items():
        pa = [_demod_phase(s[PRE:PRE + tc.N_SWEEP] * rot, sos) for s in segs]
        psi = np.stack([p for p, _ in pa])
        amp = np.stack([a for _, a in pa])
        ok = band[None, :] & (amp > amp.max(axis=1, keepdims=True) * 0.01)
        ia = int(np.argmax(band))
        psi -= 2 * np.pi * np.round((psi[:, ia] - psi[0, ia]) / (2 * np.pi))[:, None]
        W = ok.astype(float)
        S = [np.sum(W * lc[:, None] ** p, axis=0) for p in range(5)]
        T = [np.sum(W * lc[:, None] ** p * psi, axis=0) for p in range(3)]
        A = np.empty((psi.shape[1], 3, 3))
        for i in range(3):
            for j in range(3):
                A[:, i, j] = S[i + j]
        A[:, 1, 1] += 1e-3
        A[:, 2, 2] += 1e-3
        A[:, 0, 0] += 1e-9
        c = np.linalg.solve(A, np.stack(T, axis=1)[:, :, None])[:, :, 0]
        # hold the fit outside the measured stretch, then smooth along the sweep with a local cubic
        # (Savitzky-Golay), which follows the curved phase trend of the machine without bias
        c = _fill_invalid(c.T.copy(), np.repeat((S[0] > 0.5)[None, :], 3, axis=0)).T
        c[:, 0] = ss.savgol_filter(c[:, 0], w_mean, 3, mode='interp')
        c[:, 1] = ss.savgol_filter(c[:, 1], w_lev, 3, mode='interp')
        c[:, 2] = ss.savgol_filter(c[:, 2], w_lev, 3, mode='interp')
        ref = (np.stack([np.ones_like(lc), lc, lc ** 2], axis=1) @ c.T)
        d = (ref - psi) / (2 * np.pi * f_t[None, :])
        delta[name] = _fill_invalid(d, ok)
        oks[name] = ok
    wgt = np.clip(np.log(f_t / 250.0) / np.log(400.0 / 250.0), 0.0, 1.0)
    d = (1 - wgt)[None, :] * delta['lo'] + wgt[None, :] * delta['hi']
    valid = (f_t >= 40) & (f_t <= 21000)
    return d, delta['hi'], oks['hi'], valid


def _harmonic_spectra(spec, amp):
    """H_m(F) for m = 1..M_MEAS from one deconvolution spectrum, and the noise power each window collects.

    The sweep is a sine, so the m-th harmonic of the system's response comes back rotated by (-j)^(m-1)
    relative to the cosine (Chebyshev) convention the model uses: sin^n has alternating signs and swaps
    sin/cos for even orders. Rotating the deconvolution by j^(m-1) *before* cutting out order m's slot
    puts each harmonic back in the cosine convention with a compact impulse response (rotating after
    windowing would truncate the long Hilbert tails of the even orders)."""
    hq = [np.fft.irfft(spec * (1j ** q), NF) for q in range(4)]
    H = np.zeros((M_MEAS, NK // 2 + 1), complex)
    tapers = []
    for m in range(1, M_MEAS + 1):
        h = hq[(m - 1) % 4]
        c = PRE - L * math.log(m) * FS
        start = int(round(c)) - WIN_PRE
        n = WIN_PRE + harmonic_post(m)
        w = np.take(h, np.arange(start, start + n), mode='wrap').copy()
        t = np.ones(n)
        a = 96
        t[:a] = 0.5 - 0.5 * np.cos(np.linspace(0, np.pi, a))
        b = max(96, n // 5)
        t[-b:] = 0.5 + 0.5 * np.cos(np.linspace(0, np.pi, b))
        tapers.append(t)
        H[m - 1] = np.fft.rfft(w * t, NK) * np.exp(-2j * np.pi * FK * (start - c) / FS) / amp
    # noise: two stretches of the deconvolution where nothing but noise lives - between the slots of
    # orders NOISE_ORDERS (sees real output down to ~20 Hz) and 1.5..2.1 s after the linear response
    # (sees it up to ~10 kHz). The larger of the two is used, so noise is never mistaken for signal.
    ests = []
    for a, b in ((int(round(PRE - L * math.log(NOISE_ORDERS[1]) * FS)), int(round(PRE - L * math.log(NOISE_ORDERS[0]) * FS))),
                 (PRE + int(1.5 * FS), PRE + int(2.1 * FS))):
        seg = np.take(hq[0], np.arange(a, b), mode='wrap')[:NK]
        tn = np.hanning(len(seg))
        ests.append(smooth_octave(np.abs(np.fft.rfft(seg * tn, NK)) ** 2 / np.sum(tn ** 2), 1 / 3))
    p = np.maximum(ests[0], ests[1])
    noise = np.stack([p * np.sum(t ** 2) for t in tapers]) / amp ** 2
    return H, noise


def smooth_octave(p, frac):
    """Smooth a spectrum on the FK grid over +/- frac/2 octave (power average)."""
    out = np.empty_like(p)
    cs = np.concatenate([[0.0], np.cumsum(p)])
    lo = np.searchsorted(FK, FK * 2 ** (-frac / 2))
    hi = np.searchsorted(FK, FK * 2 ** (frac / 2), side='right')
    hi = np.maximum(hi, lo + 1)
    out[:] = (cs[hi] - cs[lo]) / (hi - lo)
    return out


def analyse(path, channel=None, log=print):
    test = tc.test_signal()
    fs, y = tc.read_wav(path, channel)
    warnings = []
    if fs != FS:
        y = tc.to_48k(y, fs)
        log(f'  recorded at {fs / 1000:g} kHz: converted to 48 kHz (flat to 22 kHz, nothing folds back)')
    bulk, ratio, q = align(y, test, warnings)
    if abs(ratio - 1.0) > 20e-6:
        # constant speed / clock difference: resample the capture onto the test file's clock
        log(f'  speed ratio {ratio:.6f}: resampling')
        n = int((len(y) - 1) / ratio)
        y = tc.frac_read(y, np.arange(n) * ratio)
        bulk, ratio2, q = align(y, test, warnings)
        warnings.append(f'Playback ran {100 * (ratio - 1):+.3f} % off the recording clock (two-pass capture or '
                        'clock mismatch); corrected.')
    pos = lambda start, n: tc.frac_read(y, start + bulk + np.arange(n))   # noqa: E731

    res = CaptureResult(path, len(y), channel, float(bulk), float(ratio), q, 0.0, 0.0, 0.0, 0.0, warnings=warnings)

    # ---- calibration tone (centre 26 s)
    cal = test.of('cal')[0]
    c = pos(cal.start + 2 * FS, 26 * FS)
    res.cal_dbfs = float(10 * np.log10(np.mean(c ** 2)))
    win = np.blackman(len(c))
    spec = np.abs(np.fft.rfft(c * win, 4 * len(c)))
    fgrid = np.fft.rfftfreq(4 * len(c), 1 / FS)
    k = int(np.argmax(spec * (np.abs(fgrid - 1000) < 50)))
    a, b, cc = np.log(spec[k - 1:k + 2] + 1e-30)
    res.cal_freq = float(fgrid[k] + 0.5 * (a - cc) / (a - 2 * b + cc) * (fgrid[1] - fgrid[0]))
    fund = np.sum(spec[np.abs(fgrid - res.cal_freq) < 5] ** 2)
    harm = sum(np.sum(spec[np.abs(fgrid - m * res.cal_freq) < 5] ** 2) for m in range(2, 10))
    res.cal_thd_db = float(10 * np.log10(harm / fund + 1e-30))
    if abs(res.cal_dbfs + 18.0) > 3.0:
        warnings.append(f'The 1 kHz 0 VU tone came back at {res.cal_dbfs:.1f} dBFS RMS (expected about -18). '
                        'The model is normalised to it, but check the record/playback calibration.')
    if abs(res.cal_freq / 1000.0 - 1) > 0.002:
        warnings.append(f'The 1 kHz tone plays back at {res.cal_freq:.2f} Hz: a speed error of '
                        f'{100 * (res.cal_freq / 1000 - 1):+.2f} %.')

    # ---- noise floor (start silence)
    sil = test.of('silence')[0]
    nz = pos(sil.start + int(0.3 * FS), int(1.4 * FS))
    res.noise_floor_dbfs = float(10 * np.log10(np.mean(nz ** 2) + 1e-30))

    # ---- sweeps: read, de-warp, deconvolve
    sweeps = test.of('sweep')
    seg_len = PRE + tc.N_SWEEP + POST
    raw = []
    for s in sweeps:
        raw.append(pos(s.start - PRE, seg_len))
    clip = max(np.max(np.abs(r)) for r in raw)
    if clip > 0.98:
        warnings.append(f'The capture peaks at {20*np.log10(clip):.1f} dBFS during the sweeps: reduce the playback '
                        'level into the interface, the loudest sweeps may be clipped by the converter.')
    levels = [s.level for s in sweeps]
    d, d_hi, valid_hi, valid = dewarp_corrections(raw, levels)
    log(f'  de-warp: largest correction {1e6 * np.max(np.abs(d)):.0f} us')
    X = np.fft.rfft(tc.sweep_shape(), NF)
    fx = np.fft.rfftfreq(NF, 1 / FS)
    mask = np.clip((fx - 9.0) / 2.0, 0, 1) * np.clip((22200.0 - fx) / 400.0, 0, 1)
    inv = np.conj(X) / (np.abs(X) ** 2 + 1e-9 * np.max(np.abs(X)) ** 2) * mask
    for i, s in enumerate(sweeps):
        full = np.concatenate([np.full(PRE, d[i, 0]), d[i], np.full(POST, d[i, -1])])
        yd = tc.frac_read(y, s.start - PRE + bulk + np.arange(seg_len) + full * FS)
        H, noise = _harmonic_spectra(np.fft.rfft(yd, NF) * inv, tc.vu_amp(s.level))
        res.sweeps.append(SweepResult(s.level, H.astype(np.complex64), noise.astype(np.float32),
                                      d[i].astype(np.float32), valid))
        # flutter samples from the hi band (60 Hz demodulator), 1 kHz rate, over its valid stretch
        v = np.flatnonzero(valid_hi[i])
        if len(v) > 2 * FS:
            res.flutter_sweeps.append(ss.resample_poly(d_hi[i][v[0]:v[-1] + 1], 1, 48))

    # ---- 3150 Hz flutter tone
    ft = test.of('flutter')[0]
    seg = pos(ft.start + FS, 18 * FS)
    t = np.arange(len(seg)) / FS
    z = seg * np.exp(-2j * np.pi * 3150.0 * t)
    psi, amp = _demod_phase(z, ss.butter(4, 200, fs=FS, output='sos'))
    psi -= np.polyval(np.polyfit(t, psi, 1), t)
    res.flutter_tone = ss.resample_poly(-psi / (2 * np.pi * 3150.0), 1, 48)

    # ---- pink noise sections
    for s in test.of('noise'):
        n = s.length
        res.pink.append((s.level, test.x[s.start:s.start + n].astype(np.float32), pos(s.start, n).astype(np.float32)))
    return res
