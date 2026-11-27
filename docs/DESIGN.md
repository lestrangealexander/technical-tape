# Technical Tape Special 1500: engine design

Sound Aquarium Recording. This is a tape plugin built from captures of one real machine, the
Technics RS-1500, with two tapes (ATR and Maxell) at three speeds (3¾, 7½ and 15 ips).
It follows the dynamic-convolution approach Acustica took with Nebula and Taupe, and it
changes the parts where that approach is known to be weak.

## 1. What the research says

### Nebula and VVKT

- **Diagonal Volterra kernels.** Nebula's engine is built on "diagonal Volterra kernels". It
  processes many kernels at once and morphs between them like a wavetable synth. Mapping the
  input level to a position in the kernel array is what reproduces compression and distortion
  (Sound On Sound, Nebula 3 review, 2008).
- **Diagonal Volterra in practice.** Farina and Armelloni (AES Italia 2005) describe the form
  Nebula grew out of. The system is a memoryless nonlinearity followed by a linear system with
  memory. Each diagonal kernel is convolved with the input raised to the n-th power, and the
  kernels are derived from an exponential sweep.
- **Origins.** It started as dynamic convolution: impulse responses measured at many levels
  and switched by the input amplitude. This came from Kemp (Sintefex, AES 1999, patent
  WO 98/07141) and from Farina's group. Del Sordo describes the move to Volterra kernels as a
  generalisation of dynamic convolution, and as a way around the existing patents (AES London
  lecture).
- **How many harmonics.** Acustica states "up to nine harmonics" (VVKT page). NAT, the
  sampler, allows at most 10 kernels: the fundamental and nine orders.
- **How it is sampled.** NAT plays sine sweeps at many levels. The default is 30 steps of
  1.5 dB, a 43.5 dB range, with 0 VU = −18 dBFS RMS. "Repeats" average passes, gaining up to
  6/12/18 dB of SNR for 2/4/8 repeats. Hard clipping must be avoided, and Nebula works best on
  soft saturation (Sound On Sound, "Sample your gear").
- **Known weak points.** Kernel switching is crossfaded by a "Liquidity" control. Low values
  are truer but can click or grate; high values add phasing (Sound On Sound, 2008). Envelope
  attack and release also affect preamp and tape programs. The approach assumes the
  nonlinearity has no memory, and both methods were audibly distinguishable from the hardware
  in ABX tests (Farina 2005). Tronchin & Coli (JAES 2015) confirm that the kernels depend on
  the input level, which is why they are sampled per level.
- **A cheaper way to blend levels.** PCA of the level-dependent impulse responses
  (Primavera et al., A3Lab) replaces dozens of level IRs with a few amplitude-weighted FIR
  branches. Perceived quality is kept at a fraction of the cost.

### Taupe (2018) and Taupe Studio (2025)

- **Taupe.** Taupe split tape into three processing "cells": the fundamental response, the
  harmonics, and the signal dynamics (tape compression). Wow, flutter and hiss were left out,
  and the samples were de-noised. Zero-latency versions were provided for tracking
  (Sound On Sound, Taupe review, 2019).
- **Taupe Studio.** Taupe Studio (May 2025) replaced the engine with HTET, "Hyper Tape
  Emulation Technology". It has profiles per machine and per speed. Speed changes headroom,
  frequency response and saturation. It adds wow and flutter, a hysteresis control
  (tight → loose), selectable oversampling, listening modes and clean/saturated response
  views. Later updates improved latency and fixed phase on some models. Acustica does not
  disclose how HTET works.
- **Hyper.** Acustica's Hyper engine computes gain reduction, attack and release at every
  sample. Hyper 3 adds an anti-aliasing algorithm.

### Tape physics that matters here

A tape recorder is record EQ (pre-emphasis), then hysteretic magnetisation with bias, then
playback-head losses and head bump, then playback EQ (de-emphasis), then a transport whose
delay fluctuates (Chowdhury, DAFx 2019; Mikkonen et al., DAFx 2023). The nonlinearity sits
*after* the pre-emphasis. High frequencies therefore saturate earlier, and how hard they
saturate depends on what else is in the signal.

## 2. Design decisions

| Topic | Nebula / Taupe | This plugin |
|---|---|---|
| Kernels per level | Raw harmonic IRs per dynamic step | Generalised Hammerstein branches in **Chebyshev form**, one kernel set per detector level, from synchronised swept-sine separation with the exact harmonic phase (Novak 2015). Branch `n` is fed `μL·T_n(x_p / μL)`: the kernels are then the measured harmonic responses re-indexed, with no ill-conditioned power-series inversion |
| Blending levels | Discrete steps crossfaded ("Liquidity") | Continuous and per sample: Catmull-Rom interpolation across a 1.5 dB level grid. A branch whose kernels are low-rank across levels is stored as an SVD basis (≤ 4 kernels, −70 dB residual re the fundamental); the others keep one kernel per level. No stepping, no crossfade phasing |
| Level detector | Envelope with attack/release | 1 ms look-ahead peak envelope of the *pre-emphasised* signal (analytic envelope, so no ripple on steady tones; inter-sample peaks covered), 60 ms release. The Chebyshev argument never leaves [−1, 1] |
| Frequency-dependent saturation | One broadband level | A pre-emphasis estimate `P(f)`, fitted from the captures by **iso-compression** (how much earlier each band reaches the same compression as 1 kHz), then split: the share `α` that is record pre-emphasis into the shared saturation (the highs push the whole signal into it) feeds the detector and the nonlinearity, `αP·x`; the rest is the tape's own high-frequency saturation (self-erasure), which depends on the highs' own level. A **treble stage** carries it: a second detector hears only the highs (through a high-pass at `f_c`) and indexes the treble share of the fundamental, so loud bass no longer squashes the cymbals. Single tones are exact for any `α` and `f_c`, so both are fitted on the capture's pink noise |
| Beyond the hottest capture | Last step held | A two-parameter saturation curve (onset, knee) is fitted to how the fundamental and 3rd harmonic change over the top 12 dB captured, at every frequency at once; every response is continued to at least +45 VU by that curve's own growth at each frequency's pre-emphasised level (bands already saturated go to a ceiling). Driving harder keeps compressing instead of freezing |
| Broadband behaviour | — | The capture's own pink-noise sections (0 and +3 VU, and the hot pass's from 3 dB under the normal pass's top up) are run through the model and compared in 1/3 octaves; the pre-emphasis share `α` and the detector release are chosen there (tones depend on neither), and the remaining error is reported per model |
| Harmonic orders | Up to 9 | Up to 15; each model keeps the orders that rise above −80 dB re the fundamental |
| Aliasing | Historically none handled; Hyper 3 adds anti-aliasing | Each order runs at the lowest power-of-two rate where it cannot alias (`n ≤ 2r − 1`), capped by the OS setting; 8× covers every order up to 15. Linear-phase −120 dB filters whose stopband starts *at* Nyquist |
| Dynamic range of captures | Repeats | Two passes (normal and hot). The hot pass's offsets are *measured*, not assumed: the input offset (searched over −3…+36 dB) from how the odd harmonics and the compression grow with level, the output offset from the gain at matching tape levels, then time and smooth phase per order. The normal pass stands up to its top level and the hot pass carries on above it; where both cover a level they are not averaged, because the hot pass is played back at another level and the 2nd harmonic, which partly comes from after the tape, follows the playback level |
| Wow and flutter | Left out (Taupe); added in Taupe Studio | Measured: every sweep is de-warped (so kernels stay coherent to 22 kHz) and the warp itself, with the 3150 Hz tone, becomes the flutter model: discrete components plus a continuum, resynthesised as a fluctuating delay |
| Hiss | De-noised | Left out. The noise is measured and used only to gate kernels |
| Auto gain | — | Static loudness match per INPUT setting (K-weighted pink noise through the model): no follower, no pumping |

## 3. Signal flow (per channel, real time)

```
in ─► INPUT ─► ÷ refPeak ─► P(f) (min-phase FIR) ─► x_p
                                     │
                     look-ahead peak detector ─► L(t)            (VU units, pre-emphasised)
                                     │
   in/refPeak ─► W (high-pass at f_c) ─► the same detector ─► L_hf(t)     (treble stage)
                                     │
   x_p, delayed D ─────────────────────────────► Σ_k w_1k(L)·(K_1k ∗ ·) ─┐
   x_p, delayed D ─────────────────────────────► Σ_k w_hk(L_hf)·(K_hk ∗ ·) ┤   (treble share)
   order n = 2..K, each at its rate r_n:                                   │
     up_r(x_p) ─► μL·T_n(x_p / μL) ─► down_r ─► Σ_k w_nk(L)·(K_nk ∗ ·) ───┼─► Σ ─► × refPeak
                                                                           │
   (all branches delayed to the same D = 183 samples)                     ┘
         ─► wow & flutter (fractional delay, measured spectrum, one transport for both channels)
         ─► OUTPUT × auto gain ─► out          (BYPASS: 20 ms crossfade to the latency-matched dry signal)
```

- `K_nk` are kernels (uniformly partitioned FFT convolution, one frequency-domain delay line
  per branch shared by all its kernels, so a kernel that becomes active is exact at once).
- `w_nk(L)` are Catmull-Rom weights across levels (mode "per level") or the interpolated
  basis weights (mode "basis"), per sample. Below the grid a branch scales as
  `(L/L₀)^(n−1)` (the small-signal law); above the (already extended) grid it is held.
- `T_n` by the Chebyshev recurrence. `L` for an oversampled sample is the detector value
  interpolated at that sample's input time, so every rate sees the same envelope.
- Treble stage: `W` is a 2nd-order Butterworth high-pass at `f_c`; `L_hf` is delayed exactly like
  `L`. The fundamental's response is split by `s(f) = (f/f_c)⁴ / (1 + (f/f_c)⁴)`: `1 − s` stays in
  branch 1 (indexed by `L`), `s` goes to the treble branch `K_h` (indexed by `L_hf`). A sine at
  `f` of level `A` reads `A + |W(f)|` dB on the treble detector, so `K_h` at treble level `λ` is
  the measured fundamental at input level `λ − |W(f)|`: single tones are unchanged. The treble
  kernels are 2048 taps (the highs' responses are short), about a quarter of branch 1's cost.

## 4. From captures to a model (`Tools/pipeline`)

1. **Find the test signal.** Each 20 ms chirp of the two sync bursts is matched on its own
   (band-limited, so a dull 3¾ ips top end and wow can't break it) and the 0 / 0.5 / 1.25 s
   pattern picks the right one: bulk delay and any constant speed ratio (resampled out).
2. **Calibrate.** The 1 kHz 0 VU tone of the normal pass sets the gain (unity at 1 kHz, 0 VU).
3. **De-warp each sweep.** Demodulated against its reference phase (60 Hz-wide above 250 Hz,
   12 Hz-wide below). Phase that varies smoothly with level (a weighted quadratic across the
   sweeps) is the machine; the rest is transport jitter, removed by band-limited resampling.
4. **Deconvolve.** Spectral division by the exact sweep, then a window around `−L·ln m` for
   each harmonic order `m`, after rotating the spectrum by `j^(m−1)` (the sine sweep's
   harmonic phase convention). Noise from two empty windows gives the SNR per order,
   frequency and level.
5. **Place the hot pass** (section 2): its input offset from the odd harmonics and compression,
   its output offset from the gain at matching tape levels. The normal pass covers −30…+12 VU,
   the hot pass everything above (the ATR's: up to about +35 VU).
6. **Gate.** Each harmonic is kept where it clears its noise by 12 dB; below the lowest such
   level, `H_m ∝ A^(m−1)`; never valid → zero. A gap between two valid levels (a harmonic that
   dips into the noise mid-range) is bridged in dB and phase, capped at what would have cleared
   the noise, so the valid measurements on both sides are kept. Below each order's excitation (the sweep starts at 10 Hz) responses
   are continued by a minimum-phase roll-off (a cut-off would ring longer than the kernels).
   21–22 kHz taper for every order.
7. **Fit `P(f)`** by iso-compression at −0.25 … −4 dB (monotone curves, 1/6-octave bands to
   18 kHz, peak-preserving smoothing). Falls back to flat if it doesn't explain the data.
8. **Continue above the hottest capture** (section 2): to at least +45 VU, in 3 dB steps.
8b. **Broadband fit.** Pink noise from the capture against the model (section 2): first the
    treble stage's crossover `f_c` (two candidates at fixed recorded wavelengths: 5 and 10 kHz at
    15 ips, halved for each halving of speed) together with the share `α` of `P` that feeds the
    shared saturation (½, ¼ and 0, then ±¼ around the best), then the detector release (the shortest of 20–100 ms within 0.02 dB of the best fit: a longer
    release holds the level up after transients). The error is reported in the build log and the
    QC report.
9. **Re-index, rescale, invert.** For each detector level each harmonic is taken at the input
   level that puts its pre-emphasised amplitude there, divided by `P` (magnitude and `m·∠P`),
   converted to the `μ`-scaled Chebyshev basis and un-mixed (a triangular matrix).
10. **Kernels.** Inverse FFT, 64-sample pre-delay, 8192 taps (fundamental) / 4096 (harmonics).
    Then per branch: SVD basis if ≤ 4 kernels reach −70 dB, else one kernel per level.
11. **Flutter.** Jitter from the de-warp plus the 3150 Hz tone (start/end transients trimmed):
    a spectrum of discrete peaks plus a 1/3-octave continuum.
12. **Auto gain** from the reference engine (K-weighted pink noise, INPUT −18 … +18 dB).

`extract.py` runs all of it for a folder of captures and writes, per tape and speed, the
`.ttm` and an HTML QC report (calibration, warnings, response by level, compression curves
with the continuation, harmonic curves, `P(f)`, flutter, auto gain, build log).

## 5. Zero aliasing

- Order `n` has `n×` the bandwidth. It is computed at the lowest rate `r` (2, 4 or 8) with
  `n ≤ 2r − 1`, so its images fold above Nyquist and are removed. With 15 orders, 8× is
  enough for all of them; 16× is identical in sound and CPU (offered for completeness).
- Filters: the first ×2 up / last ×2 down stage is a 339-tap Kaiser (passband 0.4535·fs,
  −120 dB from 0.5·fs: nothing folds anywhere below Nyquist). Further ×2 stages are
  half-bands (35/27/23 taps, ≈ −117 dB at the fold edge).
- Measured (15 kHz sine at +18 VU, the worst case): worst in-band spur −104 dB re the
  fundamental at 4×, 8× and 16×; at 1× the 3rd harmonic folds to 3 kHz (louder than the
  heavily compressed fundamental), at 2× −9 dB. 1× and 2× are CPU savers for gentle use.
- Latency: every chain is padded to the same 183 samples, so harmonic phase relationships
  are exact and changing OS never changes the plugin's latency.

## 6. Wow and flutter

Captured off the repro head in one pass, so what is measured is the transit-time fluctuation
between the heads, and that is what the plugin reproduces. Both channels share one transport.
A 32-tap windowed-sinc interpolator (256 phases) stays flat to 20 kHz. Amount 0–200 %
(default 100 %) in the menu and as an automatable parameter. The resynthesised spectrum was
checked against the model's: per-band variance within ±10 % from 0.2 to 100 Hz.

## 7. Sample rates, latency, CPU, memory

- **Sample rates.** Models are 48 kHz; kernels are resampled band-limited to the host rate at
  preparation (44.1–192 kHz); block 128 / 256 / 512 samples up to 50 / 100 / 192 kHz.
- **Latency** (constant, reported): block + 183 (OS chain) + 64 (kernel pre-delay) + transport
  headroom (≈ 24–40 samples). 48 kHz: ≈ 400 samples (8 ms).
- **CPU** (the C++ engine, x86-64 SSE2 baseline, one core): a 9-order model runs ≈ 10× real
  time per channel at 8× OS (≈ 20 % of a core for a stereo instance at 48 kHz); 1× ≈ 15–20×.
- **Memory.** Only the model in use is prepared; another is prepared on a background thread
  the first time it is chosen (playing carries on meanwhile), then crossfaded in. Prepared
  models are shared by every instance in the process.
- **Changing tape, speed or OS** starts a second engine core, which runs silently for the
  longest kernel and then crossfades over ≈ 40 ms: no click, no gap (measured: nothing above
  −130 dB re the signal in the 14 kHz+ band around a switch).

## 8. Model file (`.ttm`, version 3)

Little-endian tagged chunks: `TTM1` + version, `HEAD` (sample rate, reference peak, level
grid, pre-delay, look-ahead, release, captured range, margin), `PEMP` (|P| on a log-frequency
grid), `BRCH` × K (order, mode, count, length, float16 kernels, float32 basis weights for
mode 0), `TRBL` (treble crossover) and `BRHF` (the treble branch, laid out like `BRCH`) when
the model has a treble stage, `FLUT` (peaks, continuum, max delay), `AGTB` (auto-gain table)
and `META` (JSON provenance and QC). Version 2 models (no treble stage) still load. `Tools/pipeline/ttm.py` writes it; `Source/engine/TapeModel.cpp` reads it.

## 9. Validation

### On the real ATR captures (November 2026)

The reference is the capture itself. `Tests/validate_capture.py` plays the capture's own steady
tones through the model at the level they went onto the tape; the build's broadband fit plays its
pink-noise sections through it.

| | 15 ips | 7½ ips | 3¾ ips |
|---|---|---|---|
| Hot pass onto the tape / played back | +22.7 / −13.3 dB | +22.4 / −13.1 dB | +22.3 / −12.9 dB |
| Passes agree where they overlap (odd harmonics / compression) | 1.5 / 0.11 dB rms | 1.2 / 0.09 dB rms | 0.7 / 0.23 dB rms |
| Captured range | −30 … +34.7 VU | −30 … +34.4 VU | −30 … +34.2 VU |
| 1 kHz at 0 VU: level; H2, H3 capture → model; waveform match | +0.1 dB; −49.4 → −49.1, −46.5 → −46.5 dB; −69 dB | 0.0 dB; −49.3 → −49.2, −49.6 → −50.3 dB; −66 dB | −0.3 dB; −46.5 → −46.8, −43.2 → −43.6 dB; −58 dB |
| 3.15 kHz at −10 VU: H3 capture → model; waveform match | −65.0 → −64.6 dB; −73 dB | −60.5 → −60.1 dB; −62 dB | −48.2 → −46.2 dB; −56 dB |
| 1 kHz at +22 VU (nearly square): level; THD capture → model; waveform match | +0.1 dB; −8.7 → −8.9 dB; −39 dB | +0.1 dB; −12.8 → −13.0 dB; −36 dB | −0.5 dB; −17.0 → −16.7 dB; −34 dB |
| Treble stage (crossover) and shared pre-emphasis `α`, fitted | 10 kHz, 0 | 5 kHz, ¼ | 2.5 kHz, ½ |
| Pink noise, 1/3-octave rms, all five sections | 0.60 dB | 0.99 dB | 3.12 dB |
| … mean per section at 0 / +3 / +12 / +22 / +25 VU | −0.01 / −0.03 / −0.17 / −0.22 / −0.23 dB | −0.09 / −0.14 / −0.18 / +0.24 / +0.88 dB | +0.02 / −0.62 / +0.02 / +3.39 / +4.48 dB |
| … without the treble stage (best `α`) | 0.76 dB | 1.58 dB | 3.40 dB |
| … with all of `P` shared (`α` = 1), no treble stage | 1.56 dB | 4.24 dB | 6.57 dB |
| Wow & flutter measured (unweighted 2σ / weighted rms) | 0.031 / 0.010 % | 0.055 / 0.013 % | 0.147 / 0.036 % |

What the numbers say. Single tones are reproduced at every level the machine was driven to,
including the hot pass's near-square 1 kHz tone with its droop and overshoot; the treble stage
leaves them exactly as they were. Dense material is very close at 15 ips. At the slower speeds
the tape's own high-frequency saturation (self-erasure) is strong; the treble stage lets the highs
saturate on their own level, which took the 7½ ips error from 1.58 to 0.99 dB: at 0 / +3 VU its
top octave, 1–3 dB over-compressed with one detector, is now within about 1 dB. At 3¾ ips the
bands up to 6.3 kHz at 0 / +3 VU are within about 2 dB (they were up to 8 dB over-compressed),
but the 10 kHz band still comes out 5–7 dB more compressed than the machine's, and at +22 VU
and above the lows and mids 3–6 dB less compressed (section 10). The second harmonic above about 4 kHz is gated
conservatively at the slow speeds, where flutter residue raises the noise estimate: it comes out
up to 9 dB lower than measured at −10 VU (7½ ips: −58 → −67 dB at 6.3 kHz).

The C++ engine matches the Python reference on the real models to −110 dB (3¾ ips, treble
stage included) and runs them at 2.4× real time in stereo at 8× oversampling on one core of the
2.1 GHz build server (the treble stage costs about 12 %).

### On the virtual tape machine

Everything below runs against the virtual tape machine in `Tools/sim` (record EQ → saturating
magnetisation with memory and dynamic HF self-erasure → playback EQ, head bump, transport with
wow and flutter, hiss), whose captures go through exactly the pipeline the real captures will.
"Static" is the same machine without memory or erasure, which isolates the method itself.

| Model vs machine | Static 15 ips | ATR 15 ips | ATR 3¾ ips |
|---|---|---|---|
| Sines −20…+12 VU, 50 Hz–15 kHz, fundamental | ±0.05 dB | ±0.08 dB | ±0.3 dB (−0.9 dB at 10 kHz) |
| Sines at +18 VU (top capture) | ±0.08 dB | ±0.17 dB (0.9 dB at 15 kHz) | ±1 dB (−2.9 dB at 15 kHz) |
| Pink noise, 1/3-octave rms at −10 / 0 / +6 / +12 VU | 0.02 / 0.07 / 0.12 / 0.42 dB | 0.02 / 0.12 / 0.23 / 0.40 dB | 0.10 / 0.62 / 1.54 / 2.87 dB |
| Two-tone IMD sidebands (60 Hz + 7 kHz) | within 0.2 dB | within 0.6–1.2 dB | within 0.2–1.8 dB |
| Transient null at 0 / +6 / +12 VU | −25.0 / −22.6 / −20.1 dB | −23.9 / −21.2 / −18.5 dB | −11.7 / −11.9 / −12.4 dB |

What the numbers say. The method is essentially exact for a machine whose nonlinearity has
no memory (left column), including far past the hottest capture. With memory and HF
self-erasure (the ATR machine) it stays within a fraction of a dB at 15 ips. At 3¾ ips, where
the virtual machine's self-erasure is strongest, dense material at hot levels comes out
over-compressed (≈ 1.5 dB at +6 VU): the single level detector can't separate HF self-erasure
(driven by the HF content itself) from record pre-emphasis. The build measures exactly this on
every real capture (the broadband check in the log and the QC report), so it will be known for
the real machine; if it matters there, a separate self-erasure stage can be fitted from those
same pink-noise sections. The 3¾ ips null is also limited by a 2.4-sample time offset that
heavy wow leaves in the de-warped kernels (inaudible on its own).

The broadband check each build runs (pink noise at 0, +3, +6 and +9 VU from the captures, through
the model, 1/3-octave rms) on the six simulated models: ATR 15 / 7½ / 3¾ ips 0.20 / 0.85 / 1.29 dB,
Maxell 0.94 / 1.01 / 2.16 dB, always on the over-compressed side and growing with level, slowness
and the virtual tape's memory and self-erasure. That is the figure to watch on the real captures.

Before the continuation above the hottest capture, pink noise at +12 VU was off by 2.8–3.5 dB
and the +12 VU null was −1.9 dB (15 ips).

The C++ engine against the Python reference (`Tests/compare_engine.py`): residual −103 dB
re the output (float32 vs float64) at every OS factor, at 0 and +12 dB drive; −111 dB on the
real ATR 3¾ ips model.

## 10. Limits and future work

- Kernels are steady-state per level; hysteresis memory beyond what the detector release
  and the per-level phase capture is approximated (see the null figures above).
- Above the hottest capture the curve is a fit, not a measurement: the hot pass matters (the
  ATR's reach about +35 VU, past any INPUT setting with normal levels).
- Dense material at 3¾ ips (section 9): one treble detector covers everything above 2.5 kHz,
  so the top octave of dense material still over-compresses; a second treble band (a detector
  around 8 kHz) is the next step there. At +22 VU and above, dense material comes out 3–6 dB less
  compressed than the machine's.
- Above the normal pass's top (+12 VU) the 2nd harmonic is the hot pass's, which was played back
  13 dB lower; the part of it that follows the playback level (about −37 dB at +12 VU) is lower
  there than the machine at its normal playback level would give, 15 dB or more under the tape's
  odd harmonics.
- Harmonics above the 15th are not modelled: near square-wave levels (above about +25 VU) the
  17th and up are lost (kernel truncation −18 … −29 dB re the fundamental at the extreme levels).
- Hiss, print-through and dropouts are deliberately absent.
- A zero-latency tracking mode (non-uniform partitions) is a possible later addition.

## Sources

- Sound On Sound: Nebula 3 review (2008); "Sample your gear with Acustica Audio Nebula";
  Acustica Audio Taupe review (2019).
- Acustica Audio: VVKT technology page; Taupe Studio product page and KVR announcement (2025);
  Hyper technology page.
- AES London section: Del Sordo, "Volterra kernel based sampling and the future of
  convolution audio software".
- Farina & Armelloni, "Emulation of not-linear, time-variant devices by the convolution
  technique", AES Italia 2005.
- Kemp, "Analysis and simulation of non-linear audio processes using finite impulse responses
  derived at multiple impulse amplitudes", AES 106th Convention, preprint 4919 (1999).
- Novak et al., "Nonlinear system identification using exponential swept-sine signal",
  IEEE TIM 2010; Novak, Lotton & Simon, "Synchronized swept-sine: theory, application and
  implementation", JAES 2015.
- Tronchin, JAES 60(12) 2012; Tronchin & Coli, JAES 63(9) 2015.
- Primavera et al., "Approximation of dynamic convolution exploiting principal component
  analysis" (A3Lab poster).
- Chowdhury, "Real-time physical modelling for analog tape machines", DAFx 2019;
  Mikkonen et al., "Neural modeling of magnetic tape recorders", DAFx 2023.
