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


def treble_split(freqs, fc):
    """Share of the fundamental handed to the treble branch at each frequency: (f/fc)^4 / (1 + (f/fc)^4),
    the Butterworth high-pass's power response (6 % an octave below fc, 94 % an octave above)."""
    r4 = (np.asarray(freqs, float) / fc) ** 4
    return r4 / (1.0 + r4)


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
