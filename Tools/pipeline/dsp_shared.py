"""DSP definitions shared by the pipeline, the reference engine and the C++ engine.

Everything here is also implemented in Source/engine (or emitted into Source/engine/OsFilters.h by
`python3 dsp_shared.py --emit-cpp`), so the Python reference and the plugin compute the same thing.
"""
from __future__ import annotations

import math
import sys

import numpy as np
import scipy.signal as ss

STOP_DB = 120.0
KAISER_BETA = 0.1102 * (STOP_DB - 8.7)            # 12.265

# ------------------------------------------------------------------ oversampling filters
# Stage A (fs <-> 2fs): passband to 0.4535 fs, stopband from 0.5 fs (nothing folds anywhere below
# Nyquist). Stages B, C, D (further x2): Kaiser half-bands, passband 0.4535 fs, stopband far enough
# out that nothing folds into [0, 0.5 fs].
STAGE_A_TAPS = 339
HALF_BAND_TAPS = (35, 27, 23)                       # rates 4fs, 8fs, 16fs
OS_FACTORS = (1, 2, 4, 8, 16)


def stage_a():
    fc = 0.5 * (0.4535 / 2 + 0.25)                  # middle of the transition, as a fraction of 2 fs
    return ss.firwin(STAGE_A_TAPS, fc * 2.0, window=('kaiser', KAISER_BETA))


def half_band(taps):
    h = ss.firwin(taps, 0.5, window=('kaiser', KAISER_BETA))
    c = taps // 2
    for k in range(taps):                           # exact zeros at even offsets from the centre
        if k != c and (k - c) % 2 == 0:
            h[k] = 0.0
    return h


def os_stages(factor):
    """Filters for each x2 stage of an M-x chain, lowest rate first."""
    st = []
    if factor >= 2:
        st.append(stage_a())
    for i, taps in enumerate(HALF_BAND_TAPS):
        if factor >= 2 ** (i + 2):
            st.append(half_band(taps))
    return st


def os_delay(factor):
    """Group delay of up + down at this factor, in base-rate samples (may be fractional)."""
    d = 0.0
    for s, h in enumerate(os_stages(factor)):
        d += (len(h) - 1) / 2.0 * 2 / 2 ** (s + 1)  # up and down, each (N-1)/2 at rate 2^(s+1) fs
    return d


def os_pad(factor):
    """Extra delay at the M-x rate that makes the chain an integer number of base samples."""
    d = os_delay(factor)
    return int(round((math.ceil(d - 1e-9) - d) * factor)), int(math.ceil(d - 1e-9))


OS_LATENCY = max(os_pad(m)[1] for m in OS_FACTORS)  # every factor is padded to this


def upsample(x, factor):
    y = np.asarray(x, float)
    for h in os_stages(factor):
        z = np.zeros(2 * len(y))
        z[::2] = y
        y = 2.0 * np.convolve(z, h)[:len(z)]
    return y


def downsample(y, factor):
    for h in reversed(os_stages(factor)):
        y = np.convolve(y, h)[:len(y)][::2]
    return y


# ------------------------------------------------------------------ analytic envelope (IIR Hilbert pair)
# Olli Niemitalo's polyphase IIR 90-degree pair: two chains of 4 sections y = a^2 (x + y[n-2]) - x[n-2];
# the first chain's output is taken one sample later. 90 +/- 1.3 degrees from 20 Hz to 22 kHz at 48 kHz.
HILBERT_A = ((0.6923878, 0.9360654322959, 0.9882295226860, 0.9987488452737),
             (0.4021921162426, 0.8561710882420, 0.9722909545651, 0.9952884791278))


def hilbert_pair(x):
    out = []
    for chain in HILBERT_A:
        y = np.asarray(x, float)
        for a in chain:
            a2 = a * a
            y = ss.lfilter([a2, 0.0, -1.0], [1.0, 0.0, -a2], y)   # y = a^2 (x + y[n-2]) - x[n-2]
        out.append(y)
    i, q = out
    i = np.concatenate([[0.0], i[:-1]])
    return i, q


def analytic_envelope(x):
    i, q = hilbert_pair(x)
    return np.sqrt(i * i + q * q)


# ------------------------------------------------------------------ minimum-phase FIR from a magnitude
def design_minphase(freqs_hz, mag_db, fs, ntaps):
    """Minimum-phase FIR whose magnitude follows (freqs_hz, mag_db) (log-f interpolation, held at the
    ends), by the real cepstrum. Same procedure as Source/engine/Design.cpp."""
    nfft = 1 << int(math.ceil(math.log2(32 * ntaps)))
    f = np.arange(nfft // 2 + 1) * fs / nfft
    lf = np.log(np.maximum(f, 1e-3))
    m = np.interp(lf, np.log(freqs_hz), mag_db)
    logmag = m * (math.log(10.0) / 20.0)
    full = np.concatenate([logmag, logmag[-2:0:-1]])
    c = np.fft.ifft(full).real
    fold = np.zeros(nfft)
    fold[0] = c[0]
    fold[1:nfft // 2] = 2.0 * c[1:nfft // 2]
    fold[nfft // 2] = c[nfft // 2]
    spec = np.exp(np.fft.fft(fold))
    h = np.fft.ifft(spec).real[:ntaps]
    nt = max(8, ntaps // 8)
    h[-nt:] *= 0.5 + 0.5 * np.cos(np.linspace(0, np.pi, nt))
    return h


# ------------------------------------------------------------------ K-weighting (ITU-R BS.1770) at 48 kHz
K_SOS_48K = np.array([[1.53512485958697, -2.69169618940638, 1.19839281085285, 1.0, -1.69065929318241, 0.73248077421585],
                      [1.0, -2.0, 1.0, 1.0, -1.99004745483398, 0.99007225036621]])


def treble_hp(fc, fs):
    """The treble detector's weighting: a 2nd-order Butterworth high-pass at fc (RBJ form, Q = 1/sqrt 2),
    (b, a) normalised. Source/engine/Dsp.h computes the same coefficients."""
    w0 = 2.0 * math.pi * fc / fs
    cw, alpha = math.cos(w0), math.sin(w0) / (2.0 * (1.0 / math.sqrt(2.0)))
    a0 = 1.0 + alpha
    b = np.array([(1.0 + cw) / 2.0, -(1.0 + cw), (1.0 + cw) / 2.0]) / a0
    a = np.array([1.0, -2.0 * cw / a0, (1.0 - alpha) / a0])
    return b, a


def treble_lp(fc, fs):
    """2nd-order Butterworth low-pass at fc (RBJ form), the upper edge of a treble band that has another
    band above it. Source/engine/Dsp.h computes the same coefficients."""
    w0 = 2.0 * math.pi * fc / fs
    cw, alpha = math.cos(w0), math.sin(w0) / (2.0 * (1.0 / math.sqrt(2.0)))
    a0 = 1.0 + alpha
    b = np.array([(1.0 - cw) / 2.0, 1.0 - cw, (1.0 - cw) / 2.0]) / a0
    a = np.array([1.0, -2.0 * cw / a0, (1.0 - alpha) / a0])
    return b, a


def treble_filters(fc_lo, fc_hi, fs):
    """The weighting a treble band's detector hears through: a high-pass at fc_lo, and a low-pass at fc_hi
    when the band has another above it (fc_hi > 0)."""
    return [treble_hp(fc_lo, fs)] + ([treble_lp(fc_hi, fs)] if fc_hi > 0 else [])


def treble_split(freqs, fc):
    """Share of the fundamental from fc up: (f/fc)^4 / (1 + (f/fc)^4), the Butterworth high-pass's power
    response (6 % an octave below fc, 94 % an octave above); 0 for fc = 0."""
    if fc <= 0:
        return np.zeros(np.shape(freqs))
    r4 = (np.asarray(freqs, float) / fc) ** 4
    return r4 / (1.0 + r4)


def treble_shares(freqs, bands):
    """How the fundamental is split: the main branch's share and each treble band's, summing to 1.
    bands: (fc_lo, fc_hi) per band, contiguous (each fc_hi the next band's fc_lo, the top one's 0)."""
    if not bands:
        return np.ones(np.shape(freqs)), []
    main = 1.0 - treble_split(freqs, bands[0][0])
    return main, [treble_split(freqs, lo) - treble_split(freqs, hi) for lo, hi in bands]


def program_loop(sec=6.0, fs=48000, seed=3):
    """A music-like loop for calibrating auto gain: kick, bass line, snare, hats, a cymbal wash and chords
    at 100 bpm, equalised to a typical mix's long-term spectrum (per 1/3 octave re 63 Hz: 125 Hz -1, 250 Hz
    -3, 500 Hz -6, 1 kHz -9, 2 kHz -12, 4 kHz -16, 8 kHz -21, 16 kHz -29 dB), -18 dBFS RMS (0 VU).
    Pink noise has far more treble than music; this has a mix's balance and dynamics. Deterministic."""
    n = int(fs * sec)
    t = np.arange(n) / fs
    rng = np.random.default_rng(seed)
    beat = 0.6

    def env(t0, dec, length):
        e = np.zeros(n)
        i0 = int(t0 * fs)
        m = min(n - i0, int(length * fs))
        if m > 0:
            e[i0:i0 + m] = np.exp(-np.arange(m) / fs / dec)
        return e

    def filt(x, kind, f, order=2):
        return ss.sosfilt(ss.butter(order, f, kind, fs=fs, output='sos'), x)

    kick, snare, hats, bass, chords = (np.zeros(n) for _ in range(5))
    for k in range(int(sec / beat) + 1):
        tb = k * beat
        if k % 2 == 0:
            i0 = int(tb * fs)
            m = min(n - i0, int(0.35 * fs))
            if m > 0:
                tt = np.arange(m) / fs
                kick[i0:i0 + m] += (np.sin(2 * np.pi * (50 * tt + 1.8 * (1 - np.exp(-tt / 0.03))))
                                    * np.exp(-tt / 0.12))
        else:
            snare += (filt(rng.standard_normal(n), 'bandpass', [180, 6000]) * env(tb, 0.08, 0.3) * 0.6
                      + np.sin(2 * np.pi * 190 * t) * env(tb, 0.05, 0.2) * 0.4)
        for h in range(2):
            hats += filt(rng.standard_normal(n), 'highpass', 7000, 4) * env(tb + h * beat / 2, 0.03 if h else 0.06, 0.2)
    for k in range(int(sec / (2 * beat)) + 1):
        f0 = (41.2, 41.2, 55.0, 49.0)[k % 4]
        saw = sum(np.sin(2 * np.pi * f0 * m * t) / m for m in range(1, 12))
        bass += filt(saw, 'lowpass', 600) * env(k * 2 * beat, 1.2, 2 * beat)
    for k in range(int(sec / (4 * beat)) + 1):
        for f0 in ((220, 277.2, 329.6), (196, 246.9, 293.7), (174.6, 220, 261.6))[k % 3]:
            chords += sum(np.sin(2 * np.pi * f0 * m * t + rng.uniform(0, 6.28)) * 0.6 ** m
                          for m in range(1, 10)) * env(k * 4 * beat, 3.0, 4 * beat)
    wash = filt(rng.standard_normal(n), 'highpass', 4000) * (0.5 + 0.5 * np.abs(np.sin(np.pi * t / (8 * beat))))
    mix = kick + 0.5 * bass + 0.5 * snare + 0.12 * hats + 0.05 * chords + 0.035 * wash
    # equalise to the target long-term spectrum, keeping the dynamics
    cent = np.array([31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000, 12500, 16000, 20000])
    target = np.array([-4, 0, -1, -3, -6, -9, -12, -16, -21, -25, -29, -35.0])
    f, p = ss.welch(mix, fs, nperseg=8192)
    cur = np.array([10 * np.log10(np.sum(p[(f >= c * 2 ** (-1 / 6)) & (f < c * 2 ** (1 / 6))]) + 1e-30) for c in cent])
    gain_db = np.clip(target - (cur - cur[1]), -40, 40)
    fr = np.concatenate([[0], cent, [fs / 2]]) / (fs / 2)
    h = ss.firwin2(4097, fr, 10 ** (np.concatenate([[gain_db[0]], gain_db, [gain_db[-1]]]) / 20))
    mix = ss.fftconvolve(mix, h, 'same')
    return mix * 10 ** (-18 / 20) / np.sqrt(np.mean(mix ** 2))


def k_loudness_db(x):
    y = ss.sosfilt(K_SOS_48K, x)
    return 10 * np.log10(np.mean(y ** 2) + 1e-30)


# ------------------------------------------------------------------ C++ emitter
def emit_cpp(path):
    def arr(name, h):
        body = ',\n    '.join(', '.join(f'{v:.17g}' for v in h[i:i + 4]) for i in range(0, len(h), 4))
        return f'inline constexpr double {name}[{len(h)}] = {{\n    {body} }};\n'
    lines = ['// Generated by Tools/pipeline/dsp_shared.py --emit-cpp. Do not edit.',
             '#pragma once', '', 'namespace tt::osf', '{',
             f'inline constexpr int kStageATaps = {STAGE_A_TAPS};',
             f'inline constexpr int kHalfBandTaps[3] = {{ {", ".join(str(t) for t in HALF_BAND_TAPS)} }};',
             f'inline constexpr int kLatency = {OS_LATENCY};   // base-rate samples, every factor padded to this', '']
    lines.append(arr('kStageA', stage_a()))
    for i, t in enumerate(HALF_BAND_TAPS):
        lines.append(arr(f'kHalfBand{i}', half_band(t)))
    pads = ', '.join(f'{os_pad(m)[0]}' for m in OS_FACTORS)
    dels = ', '.join(f'{os_pad(m)[1]}' for m in OS_FACTORS)
    lines.append(f'/** Per factor (1, 2, 4, 8, 16): extra delay at the M-x rate, and the chain delay in base samples. */')
    lines.append(f'inline constexpr int kPadAtRate[5] = {{ {pads} }};')
    lines.append(f'inline constexpr int kChainDelay[5] = {{ {dels} }};')
    lines.append('} // namespace tt::osf')
    open(path, 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[1] == '--emit-cpp':
        emit_cpp(sys.argv[2])
    else:
        for m in OS_FACTORS:
            print(m, os_delay(m), os_pad(m))
        print('latency', OS_LATENCY)
