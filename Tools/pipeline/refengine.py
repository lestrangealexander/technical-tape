"""Reference engine: the plugin's processing in plain numpy, offline. Source/engine must match it
(Tests/compare_engine.py checks the two against each other).

Timeline (base-rate samples, A = look-ahead, D = OS latency, d0 = kernel pre-delay):
  vp(t)     = P * (x(t) * drive / refPeak)                  pre-emphasised, in VU units (0 VU sine peak = 1)
  Lraw(t)   = look-ahead peak detector of vp, covering vp up to t
  Lin(tau)  = mu * Lraw(tau + A - 1)                        envelope for the input sample at time tau
  b1(t)     = vp(t - D)
  bn(t)     = down_M( Lin * T_n(up_M(vp) / Lin) )(t)         n = 2..K, delayed to D in total
  y(t)      = refPeak * sum_n sum_r w_nr(L(t)) (B_nr * bn)(t),  L(t) = Lraw(t - D - d0 + A - 1)
Treble stage (models with one): per band k, Lk = the same detector on W_k * (x * drive / refPeak), W_k
the band's high-pass (and low-pass, below another band); the band's kernels convolve b1 like the
fundamental's, weighted by w_r(Lk(t)) with Lk delayed exactly like L.
"""
from __future__ import annotations

import math

import numpy as np
import scipy.signal as ss

import dsp_shared as ds

LOOK_AHEAD_MAX = 80                      # samples: stays inside the up-sampler's own delay at any rate


def design_pemp(model, fs):
    ntaps = int(round(256 * fs / 48000.0))
    return ds.design_minphase(np.asarray(model.pemp_freqs, float), np.asarray(model.pemp_db, float), fs, ntaps)


def look_ahead(model, fs):
    return max(1, min(LOOK_AHEAD_MAX, int(round(model.look_ahead_ms * 1e-3 * fs))))


def detector(vp, A, fs, release_ms, rms_w=0, peak_mix=0.0, advance=0):
    """Look-ahead peak envelope: max(analytic envelope, |vp|), sliding max over A, moving average over A
    (a smooth attack that completes before the peak arrives), then an exponential release.
    rms_w > 0: the envelope is the energy over the last rms_w samples instead, sqrt(2 * mean(vp^2)) (a sine
    reads its peak), plus peak_mix of the peak envelope: (1 - peak_mix) * energy + peak_mix * peak. The
    caller reads it `advance` samples late (rms_window_advance, as far as the look-ahead allows) and the peak
    envelope is delayed by as much, so it keeps its look-ahead."""
    if rms_w > 0:
        c2 = np.concatenate([[0.0], np.cumsum(np.asarray(vp, float) ** 2)])
        idx = np.arange(1, len(vp) + 1)
        e = np.sqrt(2.0 * np.maximum(c2[idx] - c2[np.maximum(idx - rms_w, 0)], 0.0) / rms_w)
        if peak_mix > 0:
            pk = np.maximum(ds.analytic_envelope(vp), np.abs(vp))
            pk = np.concatenate([np.zeros(advance), pk])[:len(pk)]
            e = (1.0 - peak_mix) * e + peak_mix * pk
    else:
        e = np.maximum(ds.analytic_envelope(vp), np.abs(vp))
    pad = np.concatenate([np.zeros(A - 1), e])
    sm = np.lib.stride_tricks.sliding_window_view(pad, A).max(axis=1)
    c = np.concatenate([[0.0], np.cumsum(np.concatenate([np.zeros(A - 1), sm]))])
    ma = (c[A:] - c[:-A]) / A
    alpha = math.exp(-1.0 / (release_ms * 1e-3 * fs))
    out = np.empty_like(ma)
    prev = 0.0
    for i in range(len(ma)):
        v = ma[i]
        prev = v if v > prev * alpha else prev * alpha
        out[i] = prev
    return np.maximum(out, 1e-6)


def rms_window(tb, fs):
    """A treble band's energy window in samples (0: its detector reads peaks). Same as rmsWindowFor()."""
    ms = getattr(tb, 'rms_ms', 0.0)
    return max(1, int(math.floor(ms * 1e-3 * fs + 0.5))) if ms > 0 else 0


def rms_window_advance(w):
    """Samples by which an energy window of w samples is read early, so that it is centred on the signal."""
    return (w - 1) // 2 if w > 0 else 0


def catmull(W, idx):
    """Rows of W (levels x count) interpolated at fractional row positions idx (clamped)."""
    nl = W.shape[0]
    i0 = np.clip(np.floor(idx).astype(int), 0, nl - 2)
    t = np.clip(idx - i0, 0.0, 1.0)[:, None]
    p0, p1 = W[np.clip(i0 - 1, 0, nl - 1)], W[i0]
    p2, p3 = W[np.clip(i0 + 1, 0, nl - 1)], W[np.clip(i0 + 2, 0, nl - 1)]
    return 0.5 * ((2 * p1) + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t ** 2
                  + (-p0 + 3 * p1 - 3 * p2 + p3) * t ** 3)


def branch_weights(branch, model, L):
    """Per-sample weights (samples, count) for one branch at detector levels L (linear, VU units).
    Inside the level grid: Catmull-Rom across levels. Below: the lowest level's kernel scaled by
    (L/L0)^(n-1) (the small-signal law). Above: the top level's kernel (self-similar)."""
    n = branch.order
    nl = model.num_levels
    W = np.eye(nl) if branch.weights is None else np.asarray(branch.weights, float)
    Ldb = 20 * np.log10(L)
    idx = (Ldb - model.level0_db) / model.level_step_db
    w = catmull(W, idx)
    lo, hi = idx < 0, idx > nl - 1
    L0 = 10 ** (model.level0_db / 20.0)
    w[lo] = W[0][None, :] * (L[lo] / L0)[:, None] ** (n - 1)
    w[hi] = W[-1]
    return w


def branch_rate(n, factor):
    """Oversampling rate a Chebyshev order runs at: the lowest power of two r with n <= 2r - 1 (alias-free),
    capped by the chosen OS factor."""
    if factor <= 1:
        return 1
    r = 2
    while 2 * r - 1 < n:
        r *= 2
    return min(factor, r)


def cheb_inputs(vp, Lraw, A, mu, K, factor):
    """Branch inputs b_n, n = 1..K, each delayed by the OS latency D in total."""
    D = ds.OS_LATENCY
    n_base = len(vp)
    out = {1: np.concatenate([np.zeros(D), vp])[:n_base]}
    by_rate = {}
    for n in range(2, K + 1):
        by_rate.setdefault(branch_rate(n, factor), []).append(n)
    for r, orders in sorted(by_rate.items()):
        if r == 1:
            d_in = A - 1                                   # the detector needs A-1 samples of look-ahead
            v = np.concatenate([np.zeros(d_in), vp])[:n_base]
            Lin = mu * Lraw                                # Lraw(t) covers vp(t - (A-1)) = v(t)
            post = lambda sig: np.concatenate([np.zeros(D - d_in), sig])[:n_base]  # noqa: E731
        else:
            pad, chain = ds.os_pad(r)
            v = ds.upsample(vp, r)
            v = np.concatenate([np.zeros(pad), v])[:len(v)]
            d_up = ds.os_delay(r) / 2 + pad / r           # base samples between vp and v
            tau = np.arange(len(v)) / r - d_up            # input time of every oversampled sample
            Lin = mu * np.interp(tau + A - 1, np.arange(n_base), Lraw)
            post = lambda sig, r=r, chain=chain: np.concatenate([np.zeros(D - chain), ds.downsample(sig, r)])[:n_base]  # noqa: E731
        u = np.clip(v / Lin, -1.0, 1.0)
        Tm, Tc = np.ones_like(u), u.copy()
        for n in range(2, max(orders) + 1):
            Tm, Tc = Tc, 2 * u * Tc - Tm
            if n in orders:
                out[n] = post(Lin * Tc)
    return out


def process(model, x, fs=48000.0, drive_db=0.0, factor=8, max_order=None, auto_gain=False, release_ms=None):
    """Offline render (no flutter). Output latency: OS_LATENCY samples, plus the kernels' pre-delay."""
    assert abs(fs - model.sample_rate) < 1, 'the reference engine runs at the model rate'
    v = np.asarray(x, float) * 10 ** (drive_db / 20.0) / model.ref_peak
    vp = np.convolve(v, design_pemp(model, fs))[:len(v)]
    A = look_ahead(model, fs)
    rel = model.release_ms if release_ms is None else release_ms
    Lraw = detector(vp, A, fs, rel)
    D = ds.OS_LATENCY
    shift = D + model.pre_delay - (A - 1)
    Lw = np.concatenate([np.full(shift, Lraw[0]), Lraw])[:len(Lraw)]
    K = len(model.branches) if max_order is None else min(max_order, len(model.branches))
    b = cheb_inputs(vp, Lraw, A, model.margin, K, factor)
    y = np.zeros(len(v))
    jobs = [(br, Lw) for br in model.branches[:K]]
    for tb in getattr(model, 'treble_bands', []):
        vh, vr = v, vp
        for fb, fa in ds.treble_filters(tb.fc_lo, tb.fc_hi, fs):
            vh, vr = ss.lfilter(fb, fa, vh), ss.lfilter(fb, fa, vr)
        lk = getattr(tb, 'leak', 0.0)
        if lk:
            vh = vh + lk * (vp - vr)          # the rest of what the shared detector hears
        rw = rms_window(tb, fs)
        adv = min(rms_window_advance(rw), shift)                # an energy window is read centred, if it can be
        Lk = detector(vh, A, fs, rel, rw, getattr(tb, 'peak_mix', 0.0), adv)
        sk = shift - adv
        jobs.append((tb.branch, np.concatenate([np.full(sk, Lk[0]), Lk])[:len(Lk)]))
    for br, Lsrc in jobs:
        w = branch_weights(br, model, Lsrc)
        for r in range(br.kernels.shape[0]):
            if not np.any(w[:, r]):
                continue
            c = ss.oaconvolve(b[br.order], br.kernels[r].astype(float))[:len(v)]
            y += w[:, r] * c
    y *= model.ref_peak
    if auto_gain:
        y *= 10 ** (auto_gain_db(model, drive_db) / 20.0)
    return y


def auto_gain_db(model, drive_db):
    g = np.asarray(model.ag_gain_db, float)
    idx = (drive_db - model.ag_drive0) / model.ag_step
    return float(np.interp(idx, np.arange(len(g)), g))


def latency(model):
    return ds.OS_LATENCY + model.pre_delay
