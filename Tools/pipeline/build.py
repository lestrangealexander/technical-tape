"""Model building: two analysed passes of one tape at one speed -> a .ttm model.

  1. place the +6 dB pass (input and output offsets measured, not assumed) and merge it with the normal pass
  2. gate each harmonic against its own noise floor; extrapolate below with H_m ~ A^(m-1)
  3. fit the record pre-emphasis P(f) from how compression onset shifts with frequency
  4. re-index to the pre-emphasised level, solve the power-series system per level
  5. kernels, then a low-rank basis per branch across levels (SVD)
  6. wow & flutter spectrum, auto-gain table
"""
from __future__ import annotations

import datetime
import math

import numpy as np
import scipy.signal as ss

import analyse as an
import dsp_shared as ds
import make_capture as mc
import refengine as re_
import tt_common as tc
import ttm

FS = tc.FS
FK = an.FK
NK = an.NK
PRE_DELAY = 64                       # kernel pre-delay (1.3 ms): room for converter pre-ringing
LEN_LINEAR = 8192                    # 171 ms
LEN_HARM = 4096                      # 85 ms
LEVEL_STEP = 1.5
K_MAX = 15
SNR_MIN_DB = 12.0


def band_power(H, centres, frac=1 / 6):
    """Mean |H|^2 in fractional-octave bands, H (..., bins)."""
    out = []
    for fc in centres:
        sel = (FK >= fc * 2 ** (-frac / 2)) & (FK < fc * 2 ** (frac / 2))
        out.append(np.mean(np.abs(H[..., sel]) ** 2, axis=-1))
    return np.stack(out, axis=-1)


def db10(p):
    return 10 * np.log10(np.maximum(p, 1e-30))


# ------------------------------------------------------------------ 1. passes
HOT_SEARCH_DB = (-3.0, 36.0)   # the hot pass's input offset is searched over this range: it need not be +6 dB
HOT_MARGIN_DB = 1.0            # the hot pass is used only this far above the normal pass's top level


def _ratio_curves(res, orders, fin):
    """Harmonic-to-fundamental ratio (dB) per level and band for each order; NaN where the harmonic is
    within 15 dB of its noise."""
    out = []
    for m in orders:
        arr = []
        for s in res.sweeps:
            h1 = band_power(s.H[0], fin)
            hm = band_power(s.H[m - 1], m * fin)
            nm = band_power(np.sqrt(s.noise[m - 1]), m * fin)
            arr.append(np.where(hm > 30 * nm, db10(hm) - db10(h1), np.nan))
        out.append(np.array(arr))
    return out


def _compression(res, fin):
    """Gain per level and band relative to the low-level gain (dB)."""
    g = np.stack([db10(band_power(s.H[0], fin)) for s in res.sweeps])
    return g - np.median(g[:5], axis=0)


def _curve_mismatch(lev0, lev6, A0, A6, d):
    """rms (dB) of the hot pass's curves against the normal pass's with the hot levels moved up by d, and
    how many points overlapped."""
    err, cnt = 0.0, 0
    for R0, R6 in zip(A0, A6):
        for b in range(R0.shape[1]):
            ok0 = ~np.isnan(R0[:, b])
            if ok0.sum() < 2:
                continue
            x = lev6 + d
            inside = (x >= lev0[ok0].min()) & (x <= lev0[ok0].max()) & ~np.isnan(R6[:, b])
            if not inside.any():
                continue
            ri = np.interp(x[inside], lev0[ok0], R0[ok0, b])
            err += float(np.sum((R6[inside, b] - ri) ** 2))
            cnt += int(inside.sum())
    return (math.sqrt(err / cnt) if cnt else math.inf), cnt


def pass_offsets(r0, r6, g0):
    """Input offset, total gain offset and time offset of the hot pass, all measured.

    The input offset is the level shift that lines up how the odd harmonics (3rd and 5th, the
    saturation's own signature) and the compression grow with level in the hot pass against the normal
    pass. Both are independent of the output gain. The second harmonic is left out: part of it comes from
    after the tape and follows the output level, which the hot pass may change as well, so it need not
    line up (its mismatch is reported). The search is wide: a hot pass is not necessarily +6 dB."""
    lev0 = np.array([s.level for s in r0.sweeps])
    lev6 = np.array([s.level for s in r6.sweeps])
    fin = np.array([250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000])
    odd0, odd6 = _ratio_curves(r0, (3, 5), fin), _ratio_curves(r6, (3, 5), fin)
    c0, c6 = _compression(r0, fin), _compression(r6, fin)
    best = None
    for d in np.arange(HOT_SEARCH_DB[0], HOT_SEARCH_DB[1] + 1e-9, 0.05):
        eo, no = _curve_mismatch(lev0, lev6, odd0, odd6, d)
        ec, _ = _curve_mismatch(lev0, lev6, [c0], [c6], d)
        if no < 20 or not math.isfinite(eo):
            continue
        e = math.sqrt(eo ** 2 + (8.0 * ec) ** 2)
        if best is None or e < best[0]:
            best = (e, float(d), eo, ec, no)
    if best is None:
        raise RuntimeError('the hot pass does not overlap the normal pass anywhere: its level offset cannot be measured')
    _, d_in, eo, ec, no = best
    e2, n2 = _curve_mismatch(lev0, lev6, _ratio_curves(r0, (2,), fin), _ratio_curves(r6, (2,), fin), d_in)
    fit = dict(odd_rms_db=round(eo, 2), compression_rms_db=round(ec, 3), points=no,
               h2_rms_db=round(e2, 1) if n2 and math.isfinite(e2) else None)
    # total gain offset: the two passes compared at the same level on the tape (so the tape's own
    # compression cancels), 200 Hz - 5 kHz
    fb = np.geomspace(200, 5000, 14)
    m0 = db10(np.stack([band_power(s.H[0] * g0, fb) for s in r0.sweeps]))
    m6 = db10(np.stack([band_power(s.H[0] * g0, fb) for s in r6.sweeps]))
    diffs = []
    for k, l6 in enumerate(lev6 + d_in):
        if lev0.min() <= l6 <= lev0.max():
            diffs.append(np.mean(m6[k] - np.array([np.interp(l6, lev0, m0[:, b]) for b in range(len(fb))])))
    d_sum = float(np.median(diffs))
    # time offset from the closest pairs at low level
    taus = []
    sel = (FK > 300) & (FK < 8000)
    for k, l6 in enumerate(lev6 + d_in):
        j = int(np.argmin(np.abs(lev0 - l6)))
        if abs(lev0[j] - l6) < 1.6 and l6 < 6:
            ph = np.unwrap(np.angle(r6.sweeps[k].H[0][sel] * np.conj(r0.sweeps[j].H[0][sel])))
            w = np.abs(r0.sweeps[j].H[0][sel])
            ph -= ph[np.argmin(np.abs(FK[sel] - 1000))]
            slope = np.sum(w * FK[sel] * ph) / np.sum(w * FK[sel] ** 2)
            taus.append(-slope / (2 * np.pi))
    tau = float(np.median(taus)) if taus else 0.0
    return d_in, d_sum, tau, fit


def merged_measurements(r0, r6, log):
    """List of (level VU, H (M, bins), noise (M, bins)) on one calibrated scale, normal pass + hot pass."""
    g0 = 10 ** ((-18.0 - r0.cal_dbfs) / 20.0)
    entries = [(s.level, s.H.astype(complex) * g0, s.noise.astype(float) * g0 ** 2) for s in r0.sweeps]
    info = {'cal_gain_db': 20 * math.log10(g0)}
    if r6 is not None:
        d_in, d_sum, tau, fit = pass_offsets(r0, r6, g0)
        d_out = d_sum - d_in
        top0 = max(s.level for s in r0.sweeps)
        log(f'  hot pass: input {d_in:+.2f} dB, output {d_out:+.2f} dB, time {tau * 1e6:+.1f} us (odd harmonics '
            f'line up to {fit["odd_rms_db"]:.2f} dB rms, compression to {fit["compression_rms_db"]:.2f} dB, '
            f'{fit["points"]} points; 2nd harmonic {fit["h2_rms_db"]} dB); it carries the model from '
            f'{top0:+.0f} to {max(s.level for s in r6.sweeps) + d_in:+.1f} VU')
        info.update(hot_in_db=d_in, hot_out_db=d_out, hot_tau_us=tau * 1e6, hot_fit=fit)
        if abs(d_in - 6.0) > 3.0:
            info.setdefault('warnings', []).append(
                f'The hot pass went onto the tape {d_in:+.1f} dB above the normal pass (+6 was planned). That is '
                f'fine: it was measured and used as such, and the model is captured up to '
                f'{max(s.level for s in r6.sweeps) + d_in:+.1f} VU.')
        if fit['odd_rms_db'] > 3.0:
            info.setdefault('warnings', []).append(
                f'The two passes agree only to {fit["odd_rms_db"]:.1f} dB rms in their odd harmonics where they '
                'overlap: check that nothing else changed between them (tape, bias, calibration).')
        g6 = g0 * 10 ** (-d_sum / 20.0)
        rot = np.exp(2j * np.pi * FK * tau)
        hot = [(s.level + d_in, s.H.astype(complex) * g6 * rot[None, :], s.noise.astype(float) * g6 ** 2) for s in r6.sweeps]
        # Each pass was de-warped against its own reference, which keeps a little slow wow (a smooth phase
        # error, common to that pass). Align the hot pass to the normal one with that smooth phase
        # difference before averaging, or the two would partly cancel at high frequencies. A time warp at
        # input frequency f moves order m's response (at F = m f) by m times the phase at f.
        lev0 = np.array([s.level for s in r0.sweeps])
        cross = np.zeros(len(FK), complex)
        for lv, H6, _ in hot:
            j = int(np.argmin(np.abs(lev0 - lv)))
            if abs(lev0[j] - lv) < 1.6 and lv <= 6.0:
                cross += H6[0] * np.conj(r0.sweeps[j].H[0].astype(complex) * g0)
        if np.any(cross):
            cs = an.smooth_octave(cross.real, 1 / 6) + 1j * an.smooth_octave(cross.imag, 1 / 6)
            dphi = np.unwrap(np.angle(cs))
            info['hot_phase_align_deg_10k'] = float(np.degrees(np.interp(10000, FK, dphi)))
            for k, (lv, H6, nz) in enumerate(hot):
                for m in range(1, H6.shape[0] + 1):
                    H6[m - 1] *= np.exp(-1j * m * np.interp(FK / m, FK, dphi))
        # The normal pass stands up to its own top level; the hot pass carries on from there. Where both
        # cover a level they are not averaged: the hot pass was played back at a different output level,
        # so what comes from after the tape (the 2nd harmonic, mostly) differs, and the normal pass's
        # setting is the machine as calibrated.
        entries.extend(e for e in hot if e[0] > top0 + HOT_MARGIN_DB)
    entries.sort(key=lambda e: e[0])
    merged, consistency = [], []
    i = 0
    while i < len(entries):
        grp = [entries[i]]
        while i + len(grp) < len(entries) and entries[i + len(grp)][0] - entries[i][0] < 1.0:
            grp.append(entries[i + len(grp)])
        i += len(grp)
        if len(grp) == 1:
            merged.append(grp[0])
            continue
        w = [1.0 / np.maximum(g[2], 1e-30) for g in grp]
        # magnitude and phase averaged separately: a small phase disagreement between passes must not
        # turn into a level loss
        mag = sum(np.abs(g[1]) * wi for g, wi in zip(grp, w)) / sum(w)
        H = mag * np.exp(1j * np.angle(sum(g[1] * wi for g, wi in zip(grp, w))))
        nz = 1.0 / sum(w)
        lv = float(np.mean([g[0] for g in grp]))
        fb = np.array([100, 1000, 5000, 10000])
        a, b = db10(band_power(grp[0][1][0], fb)), db10(band_power(grp[1][1][0], fb))
        consistency.append((lv, (b - a).tolist()))
        merged.append((lv, H, nz))
    info['pass_consistency'] = consistency
    return merged, info


# ------------------------------------------------------------------ 2. gating
def gate(levels, H, noise):
    """H, noise: (J, M, bins). Each harmonic is kept where it clears its noise floor by SNR_MIN_DB; below
    its lowest reliable level it follows H_m ~ A^(m-1); where it never clears, it fades out."""
    J, M, nb = H.shape
    Hg = H.copy()
    presence = np.ones((M, nb))
    for m in range(2, M + 1):
        p = np.stack([an.smooth_octave(np.abs(H[j, m - 1]) ** 2, 1 / 6) for j in range(J)])
        V = p > noise[:, m - 1] * 10 ** (SNR_MIN_DB / 10)
        V &= ((FK >= 1.15 * m * tc.F1) & (FK <= 21800.0))[None, :]   # where the m-th harmonic sweep exists
        anyv = V.any(axis=0)
        jhi = np.where(anyv, J - 1 - np.argmax(V[::-1], axis=0), -1)
        jlo = np.where(anyv, np.argmax(V, axis=0), J)        # the lowest level where it clears
        cols = np.arange(nb)
        # A harmonic can sink into the noise between two levels where it clears it (two passes whose
        # playback levels differ, or a null where two mechanisms cancel). Such a gap is filled between
        # the valid levels on either side (dB and phase), never above what would have cleared the noise:
        # it did not, so it was smaller than that. The measurements on both sides are kept.
        idx = np.arange(J)[:, None]
        prev_v = np.maximum.accumulate(np.where(V, idx, -1), axis=0)
        next_v = np.minimum.accumulate(np.where(V, idx, J)[::-1], axis=0)[::-1]
        gap = ~V & (prev_v >= 0) & (next_v < J)
        if gap.any():
            jj, bb = np.nonzero(gap)
            jp, jn = prev_v[jj, bb], next_v[jj, bb]
            t = (levels[jj] - levels[jp]) / (levels[jn] - levels[jp])
            hp, hn = H[jp, m - 1, bb], H[jn, m - 1, bb]
            lmag = (1 - t) * np.log(np.abs(hp) + 1e-30) + t * np.log(np.abs(hn) + 1e-30)
            cap = 0.5 * np.log(noise[jj, m - 1, bb] * 10 ** (SNR_MIN_DB / 10) + 1e-60)
            dph = np.angle(hn * np.conj(hp))
            Hg[jj, m - 1, bb] = np.exp(np.minimum(lmag, cap)) * np.exp(1j * (np.angle(hp) + t * dph))
        for j in range(J):
            below = anyv & (j < jlo)
            above = anyv & (j > jhi)
            if below.any():
                ref = H[jlo[below], m - 1, cols[below]]
                Hg[j, m - 1, below] = ref * 10 ** ((m - 1) * (levels[j] - levels[jlo[below]]) / 20.0)
            if above.any():
                Hg[j, m - 1, above] = H[jhi[above], m - 1, cols[above]]
        Hg[:, m - 1, ~anyv] = 0.0                              # never above the noise: nothing there
        pr = an.smooth_octave(anyv.astype(float), 1 / 3)
        pr = np.clip((pr - 0.25) / 0.5, 0.0, 1.0)
        presence[m - 1] = pr
        Hg[:, m - 1, :] *= pr[None, :]
    # below where each order was excited (the sweep starts at 10 Hz) there is no data: continue every
    # response from just above that edge with a minimum-phase roll-off (2nd-order Butterworth an octave
    # down) instead of a cut-off, whose ringing would be longer than the kernels
    for m in range(1, M + 1):
        fe = (1.4 if m == 1 else 1.4 * m) * tc.F1
        ie = int(np.searchsorted(FK, fe))
        lo = FK < FK[ie]
        if not lo.any():
            continue
        fc = 0.5 * FK[ie]
        s_lo = 1j * FK[lo] / fc
        s_e = 1j * FK[ie] / fc
        bw = lambda s: s ** 2 / (s ** 2 + math.sqrt(2.0) * s + 1.0)   # noqa: E731
        Hg[:, m - 1, lo] = Hg[:, m - 1, ie][:, None] * (bw(s_lo) / bw(s_e))[None, :]
    # the sweep stops at 22 kHz: nothing above 21 kHz is trusted, for any order
    edge = np.clip((22000.0 - FK) / 1000.0, 0.0, 1.0)
    edge = 0.5 - 0.5 * np.cos(np.pi * edge)
    Hg *= edge[None, None, :]
    return Hg, presence


def choose_order(levels, Hg, presence):
    """Highest harmonic order that matters: above -80 dB relative to the fundamental somewhere."""
    J, M, _ = Hg.shape
    K = 3
    fin = np.geomspace(60, 6000, 20)
    for m in range(2, M + 1):
        if presence[m - 1].max() < 0.5:
            continue
        hm = band_power(Hg[-3:, m - 1], np.minimum(m * fin, 20000))
        h1 = band_power(Hg[-3:, 0], fin)
        rel = db10(hm) - db10(h1)
        if np.nanmax(rel) > -80:
            K = m
    return int(min(K, K_MAX))


# ------------------------------------------------------------------ 3. pre-emphasis fit
def _pava_decreasing(y):
    """Closest non-increasing sequence (least squares, pool-adjacent-violators)."""
    vals, cnt = [], []
    for v in y:
        vals.append(float(v)); cnt.append(1)
        while len(vals) > 1 and vals[-2] < vals[-1]:
            v2, c2 = vals.pop(), cnt.pop()
            v1, c1 = vals.pop(), cnt.pop()
            vals.append((v1 * c1 + v2 * c2) / (c1 + c2)); cnt.append(c1 + c2)
    return np.repeat(vals, cnt)


def fit_preemphasis(levels, Hg):
    """Record pre-emphasis estimate |P(f)|: how much earlier (in input level) each frequency band reaches
    the same compression as 1 kHz. Under record EQ -> nonlinearity -> playback, compression curves at
    different frequencies are shifted copies of each other, and the shift is P. Measured at several
    iso-compression thresholds (-0.25 .. -4 dB) and averaged; no global curve fit, so no smoothing bias."""
    centres = 1000.0 * 2 ** (np.arange(-30, 26) / 6.0)          # 31 Hz .. 18 kHz, 1/6 octave
    centres = centres[(centres > 30) & (centres < 18500)]
    Mdb = db10(band_power(Hg[:, 0], centres))                   # (J, bands)
    c = Mdb - Mdb[:2].mean(axis=0, keepdims=True)
    lv = np.asarray(levels, float)
    thresholds = (-0.25, -0.5, -1.0, -2.0, -4.0)

    def crossing(curve, th):
        """Input level where a compression curve first falls to th (linear interpolation), or None."""
        below = np.flatnonzero(curve <= th)
        if len(below) == 0 or below[0] == 0:
            return None
        i = below[0]
        a0, a1, y0, y1 = lv[i - 1], lv[i], curve[i - 1], curve[i]
        return float(a0 + (th - y0) * (a1 - a0) / (y1 - y0)) if y1 != y0 else float(a1)

    # monotone (non-increasing in level) version of every band's curve, so a stray point can't fake an
    # early crossing
    for b in range(c.shape[1]):
        c[:, b] = _pava_decreasing(c[:, b])
    i1k = int(np.argmin(np.abs(centres - 1000)))
    ref = {th: crossing(c[:, i1k], th) for th in thresholds}
    P = np.full(len(centres), np.nan)
    for b in range(len(centres)):
        sh = []
        for th in thresholds:
            cb = crossing(c[:, b], th)
            if ref[th] is not None and cb is not None:
                sh.append(ref[th] - cb)
        if sh:
            P[b] = float(np.mean(sh))
    informative = ~np.isnan(P)
    info = dict(bands_used=int(informative.sum()))
    if informative.sum() < 3 or ref[-0.5] is None:
        info.update(used=False, residual_db=float('nan'), residual_flat_db=float('nan'))
        return centres, np.zeros(len(centres)), info
    idx = np.flatnonzero(informative)
    P = np.interp(np.log(centres), np.log(centres[idx]), P[idx])
    # light smoothing that keeps peaks (a record EQ's HF boost) and edges: local quadratic fits, over
    # 2/3 octave above 1 kHz and 1.5 octaves below (where onsets are noisier and record EQ is flat)
    Pn = ss.savgol_filter(P, 5, 2, mode='interp')
    Pw = ss.savgol_filter(P, 9, 2, mode='interp') if len(P) > 9 else Pn
    wlo = np.clip(np.log2(1000.0 / centres) / 1.0, 0.0, 1.0)      # 1 below 500 Hz, 0 above 1 kHz
    P = wlo * Pw + (1 - wlo) * Pn
    P = np.clip(P - P[i1k], -12.0, 24.0)

    def resid(Pv):                     # shifted 1 kHz curve vs each band
        e = []
        for b in idx:
            pred = np.interp(lv + Pv[b], lv, c[:, i1k], left=0.0, right=np.nan)
            ok = ~np.isnan(pred)
            e.append(pred[ok] - c[ok, b])
        e = np.concatenate(e)
        return float(np.sqrt(np.mean(e ** 2)))
    res_p, res_0 = resid(P), resid(np.zeros_like(P))
    used = res_p < 0.8 * res_0
    info.update(residual_db=res_p, residual_flat_db=res_0, used=bool(used))
    return centres, (P if used else np.zeros_like(P)), info


# ------------------------------------------------------------------ 3b. above the top capture
EXTRAP_TO_VU = 45.0                 # synthetic levels above the hottest capture, up to at least this level...
EXTRAP_MIN_DB = 9.0                 # ...and at least this far above it...
EXTRAP_STEP = 3.0                   # ...this far apart
FIT_LEVELS = 5                      # top captured levels that fit the continuation (12 dB)
_TH = (np.arange(256) + 0.5) * np.pi / 256


def _static(s0, k, x, deriv=False):
    """Saturating static curve x / (1 + |x/s0|^k)^(1/k) (k: knee, from soft to hard), or its slope."""
    r = np.abs(x / s0) ** k
    return (1.0 + r) ** (-1.0 / k - 1.0) if deriv else x / (1.0 + r) ** (1.0 / k)


def _gtab(s0, k, M, sdb):
    """g_m(s) = c_m(s) / s, m = 1..M, for cosine amplitudes s = 10^(sdb/20): the cosine coefficient m of the
    static curve driven by s cos(theta) for odd m, of its slope for even m (the harmonics of a small
    asymmetry, d * f'(x))."""
    s = 10 ** (np.asarray(sdb, float) / 20.0)
    xs = s[:, None] * np.cos(_TH)[None, :]
    cosm = np.cos(np.arange(1, M + 1)[:, None] * _TH[None, :])
    co = (2.0 / len(_TH)) * _static(s0, k, xs) @ cosm.T
    ce = (2.0 / len(_TH)) * _static(s0, k, xs, deriv=True) @ cosm.T
    odd = (np.arange(1, M + 1) % 2 == 1)[None, :]
    return np.where(odd, co, ce) / s[:, None]


def extrapolate_levels(levels, Hg, K, pf, pdb):
    """Synthetic levels above the hottest capture, so the model keeps saturating when it is driven harder
    than the machine was (holding the top level instead would freeze the curve: no further compression).

    Under record EQ -> static saturation -> playback, harmonic m at output frequency F for an input of
    level a is S_m(F) g_m(a |P(F/m)|), with g_m the describing functions of the saturation curve. A
    two-parameter saturator (onset s0, knee k) is fitted to how the fundamental and the 3rd harmonic change
    over the top captured levels, every frequency at once (S_m(F) solved per bin). Each response is then
    continued from its top captured value by the fitted curve's increments: continuous at the top, and
    the measured frequency shapes carry on."""
    lv = np.asarray(levels, float)
    J, M, nb = Hg.shape
    top = lv[-1]
    jf = np.arange(max(0, J - FIT_LEVELS), J)
    sgrid = np.arange(-40.0, 100.0, 0.25)

    def pdb_at(f):
        return np.interp(np.log(np.maximum(f, 1.0)), np.log(pf), pdb)

    fb = 100.0 * 2 ** (np.arange(0, int(12 * np.log2(12000 / 100.0)) + 1) / 12.0)
    bi = np.unique(np.clip(np.searchsorted(FK, fb), 0, nb - 1))
    orders = [m for m in (1, 3) if m <= max(K, 1) and np.any(Hg[jf, m - 1][:, bi])]

    def fit_err(s0db, k):
        g = _gtab(10 ** (s0db / 20.0), k, 3, sgrid)
        err = 0.0
        for m in orders:
            Hm = Hg[jf, m - 1][:, bi]
            G = np.interp(lv[jf][:, None] + pdb_at(FK[bi] / m)[None, :], sgrid, g[:, m - 1])
            S = np.sum(Hm * G, axis=0) / np.maximum(np.sum(G ** 2, axis=0), 1e-30)
            e = np.abs(Hm - S[None, :] * G) ** 2
            p = np.sum(np.abs(Hm) ** 2, axis=0)
            ok = p > 0
            err += (1.0 if m == 1 else 0.5) * float(np.mean(np.sum(e, axis=0)[ok] / p[ok]))
        return err

    best = None
    for k in (1.0, 1.5, 2.0, 3.0, 4.0, 6.0):
        for s0db in np.arange(top - 12.0, top + 18.01, 0.5):
            e = fit_err(s0db, k)
            if best is None or e < best[0]:
                best = (e, float(s0db), k)
    err, s0db, k = best
    g = _gtab(10 ** (s0db / 20.0), k, M, sgrid)
    steps = max(math.ceil(EXTRAP_MIN_DB / EXTRAP_STEP - 1e-9), math.ceil((EXTRAP_TO_VU - top) / EXTRAP_STEP - 1e-9))
    new = top + EXTRAP_STEP * np.arange(1, steps + 1)
    span_db = float(new[-1] - top)
    Hx = np.zeros((len(new), M, nb), complex)
    for m in range(1, M + 1):
        Ht = Hg[-1, m - 1]
        if not np.any(Ht):
            continue
        pd = pdb_at(FK / m)
        Gt = np.interp(top + pd, sgrid, g[:, m - 1])
        # the fitted curve's own growth from the top capture, at each bin's pre-emphasised level: the
        # fundamental keeps compressing (to a ceiling where it is already saturated), harmonics follow
        # their describing functions; where a harmonic's function is near a zero the ratio is
        # meaningless, so that bin holds
        Gall = np.stack([np.interp(a + pd, sgrid, g[:, m - 1]) for a in new])
        scale = np.max(np.abs(Gall), axis=0)
        ok = np.abs(Gt) > 0.05 * np.maximum(scale, 1e-30)
        R = np.where(ok[None, :], Gall / np.where(ok, Gt, 1.0)[None, :], 1.0)
        R = np.clip(R, -10 ** (span_db / 20.0), 10 ** (span_db / 20.0))
        Hx[:, m - 1] = Ht[None, :] * R
        Hx[:, m - 1, ~np.any(Hg[:, m - 1] != 0, axis=0)] = 0.0
    # what the fitted curve says the fundamental does at 1 kHz from the top capture upwards (for the log)
    g1 = np.interp(new, sgrid, g[:, 0]) / np.interp(top, sgrid, g[:, 0])
    info = dict(onset_vu=s0db, knee=k, fit_residual_db=float(db10(err)), levels_added=len(new),
                top_vu=float(new[-1]), gain_1k_db=[round(float(v), 2) for v in 20 * np.log10(np.maximum(g1, 1e-9))])
    return np.concatenate([lv, new]), np.concatenate([Hg, Hx]), info


# ------------------------------------------------------------------ 4. re-index and solve
def _interp_levels(levels, Hm, target, m):
    """Hm (J, bins) complex at measured levels -> values at target level (bins,) per bin.
    Cubic Hermite in dB and phase inside; H_m ~ A^(m-1) below; held above."""
    lv = np.asarray(levels, float)
    J = len(lv)
    mag = db10(np.abs(Hm) ** 2)
    ph = np.unwrap(np.angle(Hm), axis=0)
    t = np.clip(target, lv[0], lv[-1])
    i = np.clip(np.searchsorted(lv, t) - 1, 0, J - 2)
    cols = np.arange(Hm.shape[1])

    def herm(Y):
        y0, y1 = Y[i, cols], Y[i + 1, cols]
        ym, yp = Y[np.maximum(i - 1, 0), cols], Y[np.minimum(i + 2, J - 1), cols]
        x0, x1 = lv[i], lv[i + 1]
        xm, xp = lv[np.maximum(i - 1, 0)], lv[np.minimum(i + 2, J - 1)]
        h = x1 - x0
        d0 = np.where(i > 0, (y1 - ym) / np.maximum(x1 - xm, 1e-9), (y1 - y0) / h)
        d1 = np.where(i + 2 < J, (yp - y0) / np.maximum(xp - x0, 1e-9), (y1 - y0) / h)
        s = (t - x0) / h
        h00 = 2 * s ** 3 - 3 * s ** 2 + 1
        h10 = s ** 3 - 2 * s ** 2 + s
        h01 = -2 * s ** 3 + 3 * s ** 2
        h11 = s ** 3 - s ** 2
        return h00 * y0 + h10 * h * d0 + h01 * y1 + h11 * h * d1

    mdb, p = herm(mag), herm(ph)
    below = target < lv[0]
    mdb = np.where(below, mdb + (m - 1) * (target - lv[0]), mdb)
    out = 10 ** (mdb / 20.0) * np.exp(1j * p)
    out[np.all(Hm == 0, axis=0)] = 0.0                  # nothing measured there: exactly nothing
    return out


MARGIN = 1.02                        # Chebyshev interval = margin x detected peak (headroom for the envelope)


def cheb_rescale(K, mu):
    """Matrix R (K x K, 1-based orders) with T_n(x/mu) = sum_k R[n,k] T_k(x)."""
    from numpy.polynomial import chebyshev as C
    R = np.zeros((K + 1, K + 1))
    for n in range(1, K + 1):
        c = np.zeros(n + 1); c[n] = 1.0
        p = C.cheb2poly(c)                                  # power coefficients of T_n
        p = p / mu ** np.arange(n + 1)                      # T_n(x / mu)
        R[n, :n + 1] = C.poly2cheb(p)
    return R


def branch_spectra(levels, Hg, K, pemp_freqs, pemp_db, Lgrid):
    """Kernel spectra of the Chebyshev branches, (len(Lgrid), K, bins).

    Branch n's input is mu*L*T_n(vp / (mu*L)) with L the look-ahead peak envelope of the pre-emphasised
    signal vp. For a sine of pre-emphasised amplitude L this is a mix of pure harmonics, so the kernels
    are the measured harmonic responses taken at the input level that puts the pre-emphasised amplitude
    at L, divided by the pre-emphasis, then corrected for the margin mu (well conditioned, mu ~ 1)."""
    hp = ds.design_minphase(pemp_freqs, pemp_db, FS, 256)
    Pm = [ss.freqz(hp, worN=np.maximum(FK / m, 0.01), fs=FS)[1] for m in range(1, K + 1)]
    R = cheb_rescale(K, MARGIN)
    T = np.zeros((K, K))                                    # Hhat_k = sum_n mu R[n,k] K_n
    for n in range(1, K + 1):
        for k in range(1, n + 1):
            T[k - 1, n - 1] = MARGIN * R[n, k]
    Tinv = np.linalg.inv(T)
    out = np.zeros((len(Lgrid), K, len(FK)), complex)
    for li, L in enumerate(Lgrid):
        Hh = np.zeros((K, len(FK)), complex)
        for m in range(1, K + 1):
            P = Pm[m - 1]
            A = L - 20 * np.log10(np.maximum(np.abs(P), 1e-9))
            Hh[m - 1] = _interp_levels(levels, Hg[:, m - 1, :], A, m) / (np.abs(P) * np.exp(1j * m * np.angle(P)))
        out[li] = Tinv @ Hh
    return out, hp


# ------------------------------------------------------------------ 5. kernels and basis
def kernels(Kspec, K):
    """Windowed time-domain kernels per order: list of (len(Lgrid), length) arrays, and the worst
    fraction of energy left outside the window (dB, at the level where it is worst, re the fundamental)."""
    res, trunc = [], []
    e_fund = None
    for n in range(1, K + 1):
        ln = LEN_LINEAR if n == 1 else LEN_HARM
        kk = np.fft.irfft(Kspec[:, n - 1, :], NK, axis=1)
        kk = np.roll(kk, PRE_DELAY, axis=1)
        e_all = np.sum(kk ** 2, axis=1)
        w = np.ones(ln)
        a = 12
        w[:a] = 0.5 - 0.5 * np.cos(np.linspace(0, np.pi, a))
        b = ln // 8
        w[-b:] = 0.5 + 0.5 * np.cos(np.linspace(0, np.pi, b))
        kw = kk[:, :ln] * w[None, :]
        if e_fund is None:
            e_fund = e_all
        trunc.append(float(db10(np.max((e_all - np.sum(kw ** 2, axis=1)) / np.maximum(e_fund, 1e-30)))))
        res.append(kw)
    return res, trunc


def basis(kern, n, ref_norm, max_rank=4):
    """A low-rank basis across levels if one is accurate enough (errors judged against the fundamental
    kernel's size at the same level: -70 dB), otherwise None (the branch keeps one kernel per level)."""
    if n == 1:
        return None
    U, S, Vt = np.linalg.svd(kern, full_matrices=False)
    for R in range(1, min(max_rank, len(S)) + 1):
        approx = (U[:, :R] * S[:R]) @ Vt[:R]
        rel = np.linalg.norm(kern - approx, axis=1) / np.maximum(ref_norm, 1e-30)
        if db10(np.max(rel) ** 2) < -70.0:
            B = Vt[:R].copy()
            W = (U[:, :R] * S[:R]).copy()
            sc = np.max(np.abs(B), axis=1)
            B /= sc[:, None]
            W *= sc[None, :]
            return B, W, float(db10(np.max(rel) ** 2))
    return None


# ------------------------------------------------------------------ 6. flutter, auto gain
def flutter_model(results):
    fe = 1000.0
    trim = int(0.3 * fe)              # the demodulators' start-up and end transients are not transport
    segs = [s[trim:-trim] for r in results for s in r.flutter_sweeps if len(s) > 2 * trim + 1024]
    tones = [r.flutter_tone[trim:-trim] for r in results
             if r.flutter_tone is not None and len(r.flutter_tone) > 2 * trim + 1024]
    f = np.fft.rfftfreq(8192, 1.0 / fe)                       # common grid, 0.12 Hz

    def mean_psd(tracks, nper):
        acc, n = np.zeros(len(f)), 0
        for t in tracks:
            ff, p = ss.welch(t, fe, nperseg=min(nper, len(t)), detrend='linear')
            acc += np.interp(f, ff, p)
            n += 1
        return acc / n if n else None

    psd_s, psd_t = mean_psd(segs, 4096), mean_psd(tones, 8192)
    if psd_s is None and psd_t is None:
        return np.zeros((0, 2)), np.zeros((0, 2)), 0.0, {}
    if psd_s is None:
        psd_s = psd_t
    if psd_t is None:
        psd_t = psd_s
    psd = np.where(f <= 40, 0.5 * (psd_s + psd_t), psd_t)
    use = (f >= 0.2) & (f <= 200)
    f, psd = f[use], psd[use]
    # continuum: running median over 1/3 octave
    cont = np.array([np.median(psd[(f >= x * 2 ** -(1 / 6)) & (f <= x * 2 ** (1 / 6))]) for x in f])
    df = f[1] - f[0]
    peaks = []
    ex = psd - cont
    for i in range(2, len(f) - 2):
        if psd[i] > 4 * cont[i] and psd[i] == psd[i - 2:i + 3].max() and 0.3 <= f[i] <= 150:
            pw = np.sum(np.maximum(ex[i - 3:i + 4], 0)) * df
            peaks.append((float(f[i]), float(math.sqrt(pw))))
    peaks.sort(key=lambda p: -p[1])
    peaks = peaks[:10]
    grid = 0.2 * 2 ** (np.arange(0, 60) / 6.0)
    grid = grid[grid <= 200]
    cont_g = np.exp(np.interp(np.log(grid), np.log(f), np.log(np.maximum(cont, 1e-30))))
    allv = np.concatenate([np.abs(s) for s in segs + tones]) if segs or tones else np.zeros(1)
    max_delay = float(np.percentile(allv, 99.99))
    # monitored wow & flutter of the 3150 Hz tone: speed deviation seen at the repro head, %
    wf = {}
    if tones:
        t = tones[0]
        dev = -np.gradient(t) * fe
        e = int(fe)                                                # filter edges left out of the statistics
        b, a = ss.butter(1, [0.5, 30], 'bandpass', fs=fe)
        wf['unweighted_2sigma_pct'] = float(200 * np.std(ss.filtfilt(b, a, dev)[e:-e]))
        bw, aw = ss.butter(1, [2.0, 8.0], 'bandpass', fs=fe)     # rough stand-in for the IEC 4 Hz weighting
        wf['weighted_rms_pct'] = float(100 * np.std(ss.filtfilt(bw, aw, dev)[e:-e]))
        wf['delay_rms_us'] = float(1e6 * np.std(t))
    return np.array(peaks).reshape(-1, 2), np.stack([grid, cont_g], axis=1), max_delay, wf


def auto_gain_table(model, log):
    x = mc.pink(3.0, 0.0, 7)
    lin = ds.k_loudness_db(x[FS // 2:])
    drives = np.arange(-18.0, 18.01, 1.5)
    g = []
    for d in drives:
        y = re_.process(model, x, FS, drive_db=float(d), factor=4)
        g.append(lin - ds.k_loudness_db(y[FS // 2:]))
    return float(drives[0]), 1.5, np.array(g)


# ------------------------------------------------------------------ 7. broadband calibration
BB_RELEASES = (20.0, 40.0, 60.0, 100.0) # detector release times tried against the capture's pink noise, ms
BB_WINDOW = (1.5, 6.5)                  # seconds of each 8 s pink section that are compared
BB_ALPHAS = (1.0, 0.5, 0.0)             # shares of P(f) tried first as shared pre-emphasis, then +-0.25 around the best
BB_SEARCH_WINDOW = (1.5, 4.0)           # shorter stretches while searching
BB_SEARCH_RELEASE = 20.0
BB_BANDS = 1000.0 * 2 ** (np.arange(-13, 11) / 3.0)   # 1/3 octaves, 63 Hz .. 10 kHz


def _band_levels(y):
    f, p = ss.welch(y, FS, nperseg=8192)
    return np.array([10 * np.log10(np.sum(p[(f >= c * 2 ** (-1 / 6)) & (f < c * 2 ** (1 / 6))]) + 1e-30)
                     for c in BB_BANDS])


def broadband_sections(r0, r6, info, window=BB_WINDOW):
    """The captured pink-noise sections that load the tape: the normal pass's 0 and +3 VU, and the hot
    pass's sections from 3 dB under the normal pass's top level up (input, drive dB, calibrated output
    band levels)."""
    g0 = 10 ** (info['cal_gain_db'] / 20.0)
    a, b = int(window[0] * FS), int(window[1] * FS)
    out = []
    for lv, xin, yout in r0.pink:
        if lv >= 0:
            out.append((lv, np.asarray(xin, float), 0.0, _band_levels(np.asarray(yout, float)[a:b] * g0)))
    if r6 is not None and 'hot_in_db' in info:
        gh = g0 * 10 ** (-info['hot_out_db'] / 20.0)
        top0 = max(s.level for s in r0.sweeps)
        for lv, xin, yout in r6.pink:
            if lv + info['hot_in_db'] >= top0 - 3.0:
                out.append((lv + info['hot_in_db'], np.asarray(xin, float), info['hot_in_db'],
                            _band_levels(np.asarray(yout, float)[a:b] * gh)))
    return out


def pick_release(errors):
    """The shortest release whose pink-noise error is within 0.02 dB (or 5 %) of the best: a longer
    release holds the level up after transients, so it is only worth it when it clearly fits better."""
    best = min(errors.values())
    tol = max(0.02, 0.05 * best)
    return min(r for r, e in errors.items() if e <= best + tol)


def broadband_error(model, sections, release_ms, window=BB_WINDOW):
    """rms (dB) over bands and sections of model minus capture, and the mean error per section."""
    a, b = int(window[0] * FS), int(window[1] * FS)
    errs = []
    for lv, xin, drive, ref in sections:
        y = re_.process(model, xin[:b + 4096], FS, drive_db=drive, factor=2, auto_gain=False, release_ms=release_ms)
        errs.append(_band_levels(y[a:b]) - ref)
    e = np.array(errs)
    return float(np.sqrt(np.mean(e ** 2))), [float(np.mean(x)) for x in e]


# ------------------------------------------------------------------ all of it
def build(tape, speed, r0, r6, out_path, log=print, extra_meta=None, broadband=True):
    merged, info = merged_measurements(r0, r6, log)
    levels = np.array([m[0] for m in merged])
    H = np.stack([m[1] for m in merged])
    noise = np.stack([m[2] for m in merged])
    Hg, presence = gate(levels, H, noise)
    K = choose_order(levels, Hg, presence)
    log(f'  levels {levels[0]:+.1f}..{levels[-1]:+.1f} VU ({len(levels)}), harmonic orders 1..{K}')
    pf, pdb_iso, pinfo = fit_preemphasis(levels, Hg)
    log(f'  pre-emphasis fit: {"used" if pinfo["used"] else "not used (flat)"}, residual {pinfo["residual_db"]:.2f} dB '
        f'vs {pinfo["residual_flat_db"]:.2f} dB flat; P(10k) {np.interp(10000, pf, pdb_iso):+.1f} dB')

    def assemble(pdb, quiet=False):
        """Everything from the pre-emphasis on: continuation, re-indexing, kernels, branches."""
        levels_x, Hg_x, xinfo = extrapolate_levels(levels, Hg, K, pf, pdb)
        Lmin = math.floor((levels[0] + min(0.0, pdb.min())) / LEVEL_STEP) * LEVEL_STEP
        Lmax = math.ceil((levels_x[-1] + max(0.0, pdb.max())) / LEVEL_STEP) * LEVEL_STEP
        Lgrid = np.arange(Lmin, Lmax + 1e-6, LEVEL_STEP)
        Kspec, hp = branch_spectra(levels_x, Hg_x, K, pf, pdb, Lgrid)
        kern, trunc = kernels(Kspec, K)
        ref_norm = np.linalg.norm(kern[0], axis=1)
        branches, modes, errs = [], [], []
        for n in range(1, K + 1):
            bw = basis(kern[n - 1], n, ref_norm)
            if bw is None:
                branches.append(ttm.Branch(n, kern[n - 1].astype(np.float32), None))
                modes.append('per-level')
            else:
                B, W, err = bw
                branches.append(ttm.Branch(n, B.astype(np.float32), W.astype(np.float32)))
                modes.append(f'rank {B.shape[0]}')
                errs.append(err)
        model = ttm.Model(sample_rate=float(FS), ref_peak=float(tc.REF_PEAK), level0_db=float(Lgrid[0]),
                          level_step_db=LEVEL_STEP, num_levels=len(Lgrid), pre_delay=PRE_DELAY, look_ahead_ms=1.0,
                          margin=MARGIN, release_ms=60.0, captured_low_vu=float(levels[0]),
                          captured_high_vu=float(levels[-1]), pemp_freqs=pf.astype(np.float32),
                          pemp_db=pdb.astype(np.float32), branches=branches, flutter_peaks=np.zeros((0, 2)),
                          flutter_cont=np.zeros((0, 2)), flutter_max_delay=0.0, ag_drive0=-18.0, ag_step=1.5,
                          ag_gain_db=np.zeros(25), meta={})
        return model, dict(levels_x=levels_x, Hg_x=Hg_x, xinfo=xinfo, Lgrid=Lgrid, Kspec=Kspec, trunc=trunc,
                           modes=modes, errs=errs)

    # Steady tones fix the kernels, but two things only dense material can decide, because every choice
    # reproduces single tones exactly:
    #  - how much of P(f) is record pre-emphasis into the shared saturation (the highs then push the
    #    whole signal into it) and how much is the tape's own high-frequency saturation (self-erasure:
    #    the highs saturate on their own). P is scaled by a share alpha, from 1 (all pre-emphasis) to 0.
    #  - the detector's release, which decides how a peak envelope above the RMS loads the saturation.
    # Both are chosen where pink noise through the model matches the capture's own pink-noise sections.
    alpha = 1.0
    bbinfo = {'release_ms': 60.0}
    sections = broadband_sections(r0, r6, info) if broadband else []
    if sections:
        quick = broadband_sections(r0, r6, info, BB_SEARCH_WINDOW)
        tried_a, best_m = {}, None

        def try_alpha(al):
            nonlocal best_m
            m_, p_ = assemble(pdb_iso * al)
            tried_a[al] = broadband_error(m_, quick, BB_SEARCH_RELEASE, BB_SEARCH_WINDOW)
            if best_m is None or tried_a[al][0] < tried_a[best_m[0]][0]:
                best_m = (al, m_, p_)

        for al in BB_ALPHAS:
            try_alpha(al)
        a0 = best_m[0]
        for al in (a0 - 0.25, a0 + 0.25):
            if 0.0 <= al <= 1.0 and round(al, 3) not in tried_a:
                try_alpha(round(al, 3))
        alpha, model, parts = best_m
        log(f'  shared pre-emphasis: {", ".join(f"{k:g} -> {v[0]:.2f} dB" for k, v in sorted(tried_a.items()))} '
            f'(pink noise rms); using {alpha:g} of P, P(10k) {alpha * np.interp(10000, pf, pdb_iso):+.1f} dB')
        quick_rel = {rel: broadband_error(model, quick, rel, BB_SEARCH_WINDOW) for rel in BB_RELEASES}
        best = pick_release({r: v[0] for r, v in quick_rel.items()})
        err, per = broadband_error(model, sections, best)
        bbinfo = {'release_ms': best, 'error_db': round(err, 3), 'alpha': alpha,
                  'alpha_tried': {f'{k:g}': round(v[0], 3) for k, v in sorted(tried_a.items())},
                  'sections_vu': [round(s_[0], 1) for s_ in sections],
                  'mean_error_per_section_db': [round(v, 3) for v in per],
                  'release_tried': {str(int(r)): round(v[0], 3) for r, v in quick_rel.items()}}
        log(f'  broadband check: pink noise through the model vs the capture, {len(sections)} sections '
            f'{bbinfo["sections_vu"]} VU: {err:.2f} dB rms (1/3 octaves 63 Hz-10 kHz), release '
            f'{best:.0f} ms; mean per section {", ".join(f"{v:+.2f}" for v in per)} dB')
    else:
        model, parts = assemble(pdb_iso)
    pdb = pdb_iso * alpha
    pinfo['shared_fraction'] = alpha
    model.release_ms = float(bbinfo['release_ms'])
    levels_x, Hg_x, xinfo, Lgrid = parts['levels_x'], parts['Hg_x'], parts['xinfo'], parts['Lgrid']
    Kspec, trunc, modes, errs = parts['Kspec'], parts['trunc'], parts['modes'], parts['errs']
    log(f'  above the top capture: saturation onset {xinfo["onset_vu"]:+.1f} VU, knee {xinfo["knee"]}, fit '
        f'{xinfo["fit_residual_db"]:.1f} dB; continued to {xinfo["top_vu"]:+.0f} VU (1 kHz gain '
        f'{xinfo["gain_1k_db"][1]:+.1f} dB at +6, {xinfo["gain_1k_db"][-1]:+.1f} dB at '
        f'+{xinfo["top_vu"] - levels[-1]:.0f})')
    log(f'  branches: {modes}; truncation {max(trunc):.1f} dB re fundamental')
    peaks, cont, maxd, wf = flutter_model([r for r in (r0, r6) if r is not None])
    log(f'  wow & flutter: {len(peaks)} peaks, max delay {maxd * 1e6:.0f} us, {wf}')
    meta = {'tape': tape, 'speed_ips': speed, 'built': datetime.datetime.now().isoformat(timespec='seconds'),
            'sources': [r.path for r in (r0, r6) if r is not None], 'levels_vu': levels.tolist(), 'orders': K,
            'branch_modes': modes, 'basis_residual_db': errs, 'truncation_db': trunc, 'preemphasis': pinfo,
            'margin': MARGIN, 'extrapolation': xinfo, 'broadband': bbinfo,
            'wow_flutter': wf, 'passes': info,
            'calibration': {'cal_dbfs': r0.cal_dbfs, 'cal_hz': r0.cal_freq, 'cal_thd_db': r0.cal_thd_db,
                            'noise_floor_dbfs': r0.noise_floor_dbfs},
            'warnings': list(r0.warnings) + (list(r6.warnings) if r6 is not None else []) + list(info.get('warnings', []))}
    meta.update(extra_meta or {})
    model.flutter_peaks, model.flutter_cont, model.flutter_max_delay = peaks, cont, maxd
    model.meta = meta
    d0, st, g = auto_gain_table(model, log)
    model.ag_drive0, model.ag_step, model.ag_gain_db = d0, st, g.astype(np.float32)
    log(f'  auto gain: {g[0]:+.1f} dB at -18, {g[12]:+.2f} dB at 0, {g[-1]:+.1f} dB at +18 drive')
    ttm.write(out_path, model)
    return model, dict(levels=levels, H=H, Hg=Hg, noise=noise, presence=presence, K=K, pemp=(pf, pdb, pinfo),
                       Lgrid=Lgrid, info=info, Kspec=Kspec, levels_x=levels_x, Hg_x=Hg_x)
