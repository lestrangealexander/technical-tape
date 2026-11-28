"""Music captures (Tools/capture/make_music_capture.py): a mix through the machine at three levels.

For each recording: line it up with the bounce that was played (sync markers at 0:27 and 1:32), take the
machine's wow and flutter out, measure how far above or below calibration it went onto the tape and how
far the playback was turned, check the session against the main capture (the same 0 VU pink noise), and
keep the music for the fit and the checks: how the machine treats a real mix as the level changes, which
pink noise and sweeps can't pin down (the highs of a mix, and its loudness).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import scipy.signal as ss

import analyse as an
import make_capture as mc
import tt_common as tc

FS = tc.FS
T_CAL, T_LEVEL, T_PINK = (2.0, 8.0), (11.0, 6.0), (18.0, 8.0)     # start, length (s) in the template
T_SYNC = (27.0, 92.0)
T_MUSIC = (30.0, 60.0)
LEVEL_TONE_VU = -20.0


@dataclass
class MusicPass:
    path: str
    level_in_db: float           # dB above calibration onto the tape (music and tones alike)
    level_out_db: float          # dB the playback was turned up (+) or down (-) against the main capture
    x: np.ndarray                # the music as played, at calibration level (0 VU = the template's 0 VU)
    y: np.ndarray                # the machine's output, lined up, flutter out, at calibration level
    active: tuple                # (start, end) samples of the music that has signal
    pink_check_db: np.ndarray | None = None   # 0 VU pink, this session minus the main capture (0 dB pass only)
    info: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)


def _find_burst(y, lo, hi):
    """Where the sync burst (three chirps in a fixed pattern) starts in y[lo:hi], by normalised
    cross-correlation of the whole burst, refined to a fraction of a sample. A floor under the local energy
    keeps digital silence (a bounce has some) from looking like a match."""
    tpl = mc.marker_burst()
    lo, hi = max(0, int(lo)), min(len(y), int(hi))
    seg = y[lo:hi]
    m = len(tpl)
    n = tc.next_pow2(len(seg) + m)
    c = np.fft.irfft(np.fft.rfft(seg, n) * np.conj(np.fft.rfft(tpl, n)), n)[:len(seg) - m]
    e = np.concatenate([[0.0], np.cumsum(seg ** 2)])
    local = e[m:len(seg)] - e[:len(seg) - m]
    local = np.sqrt(np.maximum(local, 1e-6 * np.max(local)))
    ncc = c / (local * np.linalg.norm(tpl))
    i = int(np.argmax(ncc))
    off = 0.0
    if 0 < i < len(c) - 1:
        y0, y1, y2 = c[i - 1], c[i], c[i + 1]
        den = y0 - 2 * y1 + y2
        if den != 0:
            off = 0.5 * (y0 - y2) / den
    return lo + i + off, float(ncc[i])


def _sync(y):
    """Sample positions of the two sync bursts: the first anywhere in the first 45 s (a recording starts
    early), the second where the template puts it after the first (65 s later, +-2 s)."""
    p1, q1 = _find_burst(y, 0, 45 * FS)
    gap = (T_SYNC[1] - T_SYNC[0]) * FS
    p2, q2 = _find_burst(y, p1 + gap - 2 * FS, p1 + gap + 2 * FS + len(mc.marker_burst()))
    return [(p1, q1), (p2, q2)]


def dewarp_music(x, y):
    """The machine's output with its wow and flutter taken out, against the music as played. As
    build.dewarp_pink (the delay tracked in the 1-5 kHz band), but music has quiet stretches: the phase is
    taken only where that band carries signal (20 ms averages), unwrapped across the gaps and interpolated
    through them."""
    sos = ss.butter(4, [1000.0, 5000.0], 'bandpass', fs=FS, output='sos')
    xa, ya = ss.hilbert(ss.sosfiltfilt(sos, x)), ss.hilbert(ss.sosfiltfilt(sos, y))
    k = int(0.020 * FS)
    box = np.ones(k) / k
    step = FS // 1000
    c = np.convolve(ya * np.conj(xa), box, 'same')[::step]
    w = np.convolve(np.abs(xa) ** 2, box, 'same')[::step]
    ok = w > 1e-3 * np.percentile(w, 99)
    if ok.sum() < 100:
        return y
    td = np.flatnonzero(ok)
    ph = np.unwrap(np.angle(c[ok]))
    tau = -(ph - np.median(ph)) / (2 * np.pi * 2485.0)
    tau = np.clip(tau, -0.002, 0.002)
    tau_ms = np.interp(np.arange(len(c)), td, tau)
    tau_ms = ss.sosfiltfilt(ss.butter(2, 60.0, fs=1000.0, output='sos'), tau_ms)
    t = np.interp(np.arange(len(y)) / step, np.arange(len(c)), tau_ms)
    return tc.frac_read(y, np.arange(len(y)) + t * FS)


def _read48(path, channel):
    fs, v = tc.read_wav(path, channel)
    return v if fs == FS else tc.to_48k(v, fs)


def _tone_level(v, f0=1000.0):
    """Level (dB re a full-scale sine's rms... any fixed reference) of the f0 component and of its 3rd
    harmonic, Blackman-windowed FFT."""
    w = np.blackman(len(v))
    spec = np.abs(np.fft.rfft(v * w, 4 * len(v)))
    fr = np.fft.rfftfreq(4 * len(v), 1 / FS)
    p = lambda f: np.sum(spec[np.abs(fr - f) < 6.0] ** 2)        # noqa: E731
    return 10 * np.log10(p(f0) + 1e-30), 10 * np.log10(p(3 * f0) + 1e-30) - 10 * np.log10(p(f0) + 1e-30)


def _h3_curve(r0, r6, info):
    """The main capture's 1 kHz 3rd harmonic re the fundamental against level (merged passes)."""
    import build as B
    merged, _ = B.merged_measurements(r0, r6, lambda s='': None)
    lv = np.array([m[0] for m in merged])
    H = np.stack([m[1] for m in merged])
    k = int(np.argmin(np.abs(B.FK - 1000.0)))
    return lv, 20 * np.log10(np.abs(H[:, 2, k]) / np.abs(H[:, 0, k]) + 1e-12)


def analyse_pass(source_path, rec_path, r0, r6, channel=None, log=print, level_out_db=None):
    """One music recording against the bounce that was played, on the main capture's calibration.
    level_out_db: how far the return was turned down for this pass, if known (otherwise measured)."""
    import build as B
    xs = _read48(source_path, None)
    y = _read48(rec_path, channel)
    warns = []
    sx = _sync(xs)
    sy = _sync(y)
    for (p, q), which in zip(sx + sy, ('source', 'source', 'recording', 'recording')):
        if q < 0.6:
            warns.append(f'sync markers weak in the {which} (match {q:.2f})')
    ratio = (sy[1][0] - sy[0][0]) / (sx[1][0] - sx[0][0])
    if abs(ratio - 1.0) > 20e-6:
        y = tc.frac_read(y, np.arange(int((len(y) - 1) / ratio)) * ratio)
        sy = _sync(y)
        warns.append(f'playback ran {100 * (ratio - 1):+.3f} % off the recording clock; corrected')
    shift = sy[0][0] - sx[0][0]                                        # y sample = x sample + shift
    yy = tc.frac_read(y, np.arange(len(xs)) + shift)
    g0 = 10 ** ((-18.0 - r0.cal_dbfs) / 20.0)                          # the main capture's playback calibration
    yy *= g0
    seg = lambda v, t: v[int(t[0] * FS) + FS // 2:int((t[0] + t[1]) * FS) - FS // 2]   # noqa: E731
    lx_lvl, _ = _tone_level(seg(xs, T_LEVEL))
    ly_lvl, _ = _tone_level(seg(yy, T_LEVEL))
    gain = ly_lvl - lx_lvl                                             # level in + level out (+ drift)
    lx_cal, _ = _tone_level(seg(xs, T_CAL))
    _, h3 = _tone_level(seg(yy, T_CAL))
    played_db = lx_cal - _tone_level(mc.tone(1000, T_CAL[1], 0.0)[FS // 2:-FS // 2])[0]   # the bounce's own level
    lv, h3c = _h3_curve(r0, r6, None)
    # The level onto the tape. A pass turned down digitally (the low one) has its return untouched: the
    # level tone gives it. A hot pass may have had its return turned down: there the calibration tone's
    # 3rd harmonic, which rises steeply with level, places it on the main capture's own curve. Where the
    # two agree to 1.5 dB the return was untouched and the tone (more precise) is used.
    vu_from_h3 = None
    if gain > -5.0:
        sel = (lv >= -6.0) & np.isfinite(h3c)
        order = np.argsort(h3c[sel])
        vu_from_h3 = float(np.interp(h3, h3c[sel][order], lv[sel][order])) - played_db
    if level_out_db is not None:
        level_in = gain - level_out_db
    elif vu_from_h3 is None or abs(vu_from_h3 - gain) <= 3.0:
        level_in = gain
    else:
        level_in = vu_from_h3
        warns.append(f'the return seems to have been turned by {gain - level_in:+.1f} dB (from the tone\'s 3rd '
                     'harmonic, +-1.5 dB): say how much, if it was')
    level_out = gain - level_in
    xm = xs[int(T_MUSIC[0] * FS):int((T_MUSIC[0] + T_MUSIC[1]) * FS)] * 1.0
    ym = yy[int(T_MUSIC[0] * FS):int((T_MUSIC[0] + T_MUSIC[1]) * FS)] * 10 ** (-level_out / 20.0)
    e = np.convolve(xm ** 2, np.ones(FS // 10) / (FS // 10), 'same')
    act = np.flatnonzero(e > 1e-3 * np.max(e))
    active = (int(act[0]), int(act[-1])) if len(act) else (0, len(xm))
    ym = dewarp_music(xm, ym)
    pink_check = None
    if abs(level_in) < 2.0:
        a, b = int((T_PINK[0] + 1.0) * FS), int((T_PINK[0] + 7.0) * FS)
        main = [p_ for p_ in r0.pink if abs(p_[0]) < 0.1]
        if main:
            pm = B._band_levels(np.asarray(main[0][2], float)[FS:7 * FS] * g0)
            pink_check = B._band_levels(yy[a:b] * 10 ** (-level_out / 20.0)) - pm
    log(f'  {rec_path.rsplit("/", 1)[-1]}: onto the tape {level_in:+.1f} dB, playback {level_out:+.1f} dB '
        f'(tone {gain:+.1f} dB, its 3rd harmonic {h3:.1f} dB), music {(active[1] - active[0]) / FS:.0f} s'
        + (f'; session vs main capture (0 VU pink), 1/3 octaves: mean {np.mean(pink_check):+.2f}, '
           f'worst {pink_check[np.argmax(np.abs(pink_check))]:+.1f} dB' if pink_check is not None else ''))
    return MusicPass(rec_path, float(level_in), float(level_out), xm, ym, active, pink_check,
                     {'tone_gain_db': gain, 'cal_h3_db': h3, 'speed_ratio': ratio, 'played_db': played_db},
                     warns)


def bands(y):
    import build as B
    return B._band_levels(y)


def compare(model, passes, factor=2, seconds=None, release_ms=None):
    """A model against the machine on the music of each pass: per 1/3 octave (63 Hz-16 kHz) model minus
    machine, total power and the part that follows the input; the top end (6-13 kHz) and the presence
    (2.5-5 kHz) re the mids (250 Hz-2 kHz), output re input, machine and model; and the loudness
    (K-weighted) of both, dB re the music as played."""
    import build as B
    import dsp_shared as ds
    import refengine as re_
    cent = np.array(B.BB_BANDS)
    mid = (cent >= 250) & (cent <= 2000)
    top = (cent >= 6000) & (cent <= 13000)
    pres = (cent >= 2500) & (cent <= 5000)
    lat = re_.latency(model)
    out = []
    for p in passes:
        a, b = p.active
        if seconds:
            b = min(b, a + int(seconds * FS))
        x = p.x[a:b]
        ys = re_.process(model, np.concatenate([x, np.zeros(lat + 64)]), FS, drive_db=p.level_in_db, factor=factor,
                         auto_gain=False, release_ms=release_ms)[lat:lat + len(x)]
        ym = p.y[a:b]
        xin = x * 10 ** (p.level_in_db / 20.0)
        bx, bm, bs = B._band_levels(xin), B._band_levels(ym), B._band_levels(ys)
        cm, cs = B.coherent_levels(xin, ym), B.coherent_levels(xin, ys)
        bal = lambda v, sel: float(np.mean(v[sel] - bx[sel]) - np.mean(v[mid] - bx[mid]))   # noqa: E731
        out.append(dict(level=p.level_in_db, total=bs - bm, coherent=cs - cm,
                        top=(bal(bm, top), bal(bs, top)), presence=(bal(bm, pres), bal(bs, pres)),
                        top_coherent=(bal(cm, top), bal(cs, top)),
                        loudness=(ds.k_loudness_db(ym) - ds.k_loudness_db(x), ds.k_loudness_db(ys) - ds.k_loudness_db(x))))
    return out
