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
| Frequency-dependent saturation | One broadband level | A pre-emphasis estimate `P(f)`, fitted from the captures by **iso-compression** (how much earlier each band reaches the same compression as 1 kHz), then split: the share `α` that is record pre-emphasis into the shared saturation (the highs push the whole signal into it) feeds the detector and the nonlinearity, `αP·x`; the rest is the tape's own high-frequency saturation (self-erasure), which depends on the highs' own level. A **treble stage** carries it: band detectors hear the highs (through a high-pass at `f_c`), optionally plus a *leak* of everything else, and index the treble share of the fundamental, so the highs saturate on their own level. A band detector reads the highs' **energy** over a few milliseconds (scaled so a sine reads its peak) mixed with a share of their peaks: a peak detector reads a mix's highs (cymbals, hats, sibilance) as far louder than a sine of the same energy, so at 3¾ ips it dulled them far sooner than the machine does. Single tones are exact for any `α`, `f_c`, leak and detector, so all of them are fitted on music recorded through the machine (and the capture's pink noise) |
| Beyond the hottest capture | Last step held | A two-parameter saturation curve (onset, knee) is fitted to how the fundamental and 3rd harmonic change over the top 12 dB captured, at every frequency at once; every response is continued to at least +45 VU by that curve's own growth at each frequency's pre-emphasised level (bands already saturated go to a ceiling). Driving harder keeps compressing instead of freezing |
| Broadband behaviour | — | **Music through the machine** at three levels (−10, 0 and about +19 dB onto the tape; `docs/CAPTURES.md`) is the reference: the bounce that was played goes through the model at each pass's measured level and the two are compared in 1/3 octaves from 63 Hz to 16 kHz (the machine's noise and the session's drift taken out), on how far the model falls outside the machine: above its total power, or below the part of its output that follows the input. The treble layout, its detectors and leak, the pre-emphasis share `α` and the detector release are chosen there (tones depend on none of them), with the capture's own pink-noise sections (−10 VU to the hot pass's first) counting 0.3 as much, scored the same way. On tape a growing share of dense highs stops following the input as the level rises (distortion products, the tape's own modulation of the highs: at 3¾ ips about half the power above 6 kHz on pink noise at 0 VU), which the model can't make; in between, the model is as close as it can get. Fitted on pink noise alone (0.9.1–0.9.3), the treble was right for pink noise and wrong for music: too bright when matched on total power, then, matched on the part that follows the input, dulled far too fast as INPUT went up at 3¾ ips |
| Harmonic orders | Up to 9 | Up to 15; each model keeps the orders that rise above −80 dB re the fundamental |
| Aliasing | Historically none handled; Hyper 3 adds anti-aliasing | Each order runs at the lowest power-of-two rate where it cannot alias (`n ≤ 2r − 1`), capped by the OS setting; 8× covers every order up to 15. Linear-phase −120 dB filters whose stopband starts *at* Nyquist |
| Dynamic range of captures | Repeats | Two passes (normal and hot). The hot pass's offsets are *measured*, not assumed: the input offset (searched over −3…+36 dB) from how the odd harmonics and the compression grow with level, the output offset from the gain at matching tape levels, then time and smooth phase per order. The normal pass stands up to its top level and the hot pass carries on above it; where both cover a level they are not averaged, because the hot pass is played back at another level and the 2nd harmonic, which partly comes from after the tape, follows the playback level |
| Wow and flutter | Left out (Taupe); added in Taupe Studio | Measured: every sweep is de-warped (so kernels stay coherent to 22 kHz) and the warp itself, with the 3150 Hz tone, becomes the flutter model: discrete components plus a continuum, resynthesised as a fluctuating delay |
| Hiss | De-noised | Left out. The noise is measured and used only to gate kernels |
| Auto gain | — | Static loudness match per INPUT setting: the mix recorded through the machine (or, without one, a music-like loop at the same level) through the model, K-weighted, held at its loudness at INPUT 0, so AG on and off sound the same at INPUT 0 and only the saturation changes as INPUT turns. No follower, no pumping. (Up to 0.9.3 it was made on pink noise, whose treble collapses with drive at the slow speeds: music got up to 13 dB louder across INPUT at 3¾ ips) |

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
- Treble stage: one or more bands, lowest first. Band `k`'s detector hears the signal through
  `W_k`, a 2nd-order Butterworth high-pass at its `f_lo` and, when another band sits above it, a
  low-pass at its `f_hi`, plus the leak `β` times the rest of what the shared detector hears (`x_p`
  with that band taken out: a loud low end and midrange squash the highs too), so
  `W_k = W_b + β·P·(1 − W_b)` with `W_b` the band's filters. It reads
  `(1 − κ)·√(2·mean(v²)) + κ·peak(v)`: the energy of what it hears over a window of `T` ms (a running
  sum), scaled so a sine reads its peak, plus a share `κ` of the peak envelope; then the same
  look-ahead, smoothing and release as `L`. `T = 0` is the peak detector alone (0.9.3). Its level
  `L_k` is delayed like `L`, less half the window so that the window is centred, as far as the
  look-ahead allows (at 48 kHz, windows up to 8 ms; a longer one lags by the rest); the peak envelope
  is delayed by the same, so its share keeps all of its look-ahead. The fundamental's response
  is split with `s(f; f_c) = (f/f_c)⁴ / (1 + (f/f_c)⁴)`: `1 − s(f; f_lo,1)` stays in branch 1
  (indexed by `L`), band `k` takes `s(f; f_lo) − s(f; f_hi)` into its own kernels (indexed by
  `L_k`). A sine at `f` of level `A` reads `A + |W_k(f)|` dB on band `k`'s detector, so its kernel
  at band level `λ` is the measured fundamental at input level `λ − |W_k(f)|`: single tones are
  unchanged. The treble kernels are 512, 1024 or 2048 taps: the shortest that leaves no more of a
  band's response outside its window than the band-edge pre-ringing ahead of the pre-delay, which
  no causal window keeps (at 3¾ ips in 0.9.2, 1024 instead of 2048 gave the same pink noise to
  0.01 dB and a difference 36 dB down in the highs). A band costs an eighth
  to a quarter of branch 1. Only the selected speed's model runs, so a band costs CPU only at its
  own speed.

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
8b. **Treble fit on music.** The music passes (`Tools/pipeline/music.py`): each recording is lined
    up with the bounce by its two sync bursts, its wow and flutter taken out (the delay tracked in
    1–5 kHz where the music has signal, interpolated through the gaps), its level onto the tape read
    from the −20 VU tone, its noise from the silence before the first burst, and the session's drift
    from the 0 dB pass's 0 VU pink noise against the main capture's. Three 8 s excerpts of each pass
    go through the C++ engine at that pass's level and are compared 1/3 octave by 1/3 octave with the
    machine, on how far the model falls outside it: above the machine's total power (its noise
    subtracted, the drift taken off; bands within 6 dB of the noise left out), or below the part of
    the machine's output that follows the input. At normal levels the two are the same and this is
    plain total power; driven hard, the machine's top octave is mostly its own distortion products
    and modulation (at +19 dB at 7½ ips, 3–6 dB of 10–16 kHz doesn't follow the input), which the
    model can't make, and it isn't asked to make up for them with clean treble (fitted on total power,
    6–10 kHz came out 2 dB too bright there to pay for a dark 16 kHz). Above 2 kHz an error where the
    model is brighter counts twice at 3¾ ips and 1.5 times at 7½ (rather a shade dark than a shade
    bright). The capture's pink noise (−10 VU to the hot pass's first section) counts 0.3 as much,
    scored the same way. The search: the treble layout (one band from 5 or 10 kHz at 15 ips; at 7½
    ips one band from 2.5 or 5 kHz, or two, 5–10 kHz and 10 kHz up; at 3¾ ips two bands, 5–10 kHz and
    10 kHz up or 4–8 kHz and 8 kHz up) at `α` = ¼, ½, ¾ with four band detectors (peaks; energy over 5 ms with ¼, ½
    or ¾ of the peaks), then from the best the peak share in steps of ⅛ and the window along 1, 2, 5,
    10, 20 ms for as long as either improves, `α` ±⅛, the leak (0.1, 0.2, 0.35: only the treble
    kernels change) and the release (20–100 ms). The build log and the QC report give the error per
    pass, the top end and presence against the machine's, the per-band errors, and the pink noise's.
    Without music passes the treble stage is fitted on the pink noise's coherent power alone, as up to
    0.9.3.
9. **Re-index, rescale, invert.** For each detector level each harmonic is taken at the input
   level that puts its pre-emphasised amplitude there, divided by `P` (magnitude and `m·∠P`),
   converted to the `μ`-scaled Chebyshev basis and un-mixed (a triangular matrix).
10. **Kernels.** Inverse FFT, 64-sample pre-delay, 8192 taps (fundamental) / 4096 (harmonics).
    Then per branch: SVD basis if ≤ 4 kernels reach −70 dB, else one kernel per level.
11. **Flutter.** Jitter from the de-warp plus the 3150 Hz tone (start/end transients trimmed):
    a spectrum of discrete peaks plus a 1/3-octave continuum.
12. **Auto gain**: the music passes' mix (without music passes, `ds.program_loop`, a music-like loop
    at the mix's level, −22 dBFS RMS) through the model at INPUT −18 … +18 dB in 1.5 dB steps;
    K-weighted loudness held at INPUT 0's.

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

## 8. Model file (`.ttm`, version 6)

Little-endian tagged chunks: `TTM1` + version, `HEAD` (sample rate, reference peak, level
grid, pre-delay, look-ahead, release, captured range, margin), `PEMP` (|P| on a log-frequency
grid), `BRCH` × K (order, mode, count, length, float16 kernels, float32 basis weights for
mode 0), a `TRBL` (the band's high-pass and low-pass frequencies, its leak, its detector's energy
window in ms and peak share) and `BRHF` (its kernels, laid out like `BRCH`) pair per treble band,
`FLUT` (peaks, continuum, max delay), `AGTB` (auto-gain table) and `META` (JSON provenance and QC,
with the music fit's figures). Version 2 (no treble stage), 3 (one treble band, high-pass only), 4
(no leak) and 5 (peak detectors) models still load. `Tools/pipeline/ttm.py` writes it; `Source/engine/TapeModel.cpp` reads it.

## 9. Validation

### On the real captures (November 2026)

The reference is the machine itself: its own steady tones (`Tests/validate_capture.py` plays the
capture's tones through the model at the level they went onto the tape), a mix recorded through it at
three levels (section 4, step 8b), and the capture's pink noise.

**Music.** All 60 s of the mix through the C++ engine at 8× (no flutter), against the machine's
recording with its flutter, its noise and the session's drift taken out; 1/3 octaves from 63 Hz to
16 kHz, total power, model minus machine. The INPUT rows play the same mix through each model at
every setting from −18 to +18 dB in 1 dB steps.

| | ATR 15 | ATR 7½ | ATR 3¾ | Maxell 15 | Maxell 7½ | Maxell 3¾ |
|---|---|---|---|---|---|---|
| Music onto the tape (dB re the main capture) | −10.3 / −0.3 / +18.8 | −10.2 / −0.1 / +19.0 | −10.0 / −0.2 / +19.1 | −10.3 / −0.6 / +17.8 | −10.0 / +0.2 / +18.1 | −9.8 / −0.2 / +17.9 |
| All of the music, 1/3-octave rms: 0.9.3 → 0.9.4 | 0.49 → 0.43 dB | 0.77 → 0.46 dB | 4.81 → 1.28 dB | 0.47 → 0.45 dB | 0.88 → 0.51 dB | 2.17 → 1.02 dB |
| … per pass, 0.9.4 | 0.33 / 0.28 / 0.60 dB | 0.16 / 0.22 / 0.76 dB | 0.28 / 0.68 / 2.10 dB | 0.31 / 0.28 / 0.65 dB | 0.16 / 0.40 / 0.78 dB | 0.19 / 0.81 / 1.57 dB |
| Top end (6–13 kHz re 250 Hz–2 kHz) at 0 dB: machine / 0.9.3 / 0.9.4 | +1.5 / +1.2 / +1.2 | +1.9 / +1.8 / +1.8 | +4.3 / +1.5 / +3.8 | +0.7 / +0.5 / +0.5 | −0.6 / −0.9 / −1.0 | +0.6 / −1.2 / −0.1 |
| … driven hard | +1.6 / −0.0 / +0.5 | +1.4 / +1.6 / +2.4 | +0.4 / −14.3 / +0.0 | +0.5 / +0.2 / +0.6 | −1.2 / −1.1 / +0.0 | −4.1 / −10.2 / −3.6 |
| Presence (2.5–5 kHz) driven hard: machine / 0.9.3 / 0.9.4 | +1.0 / +1.0 / +0.9 | +0.4 / +0.1 / −0.2 | +0.4 / −2.5 / −1.3 | +0.8 / +1.0 / +0.9 | −0.4 / −0.6 / −0.6 | −1.0 / −0.6 / −0.6 |
| Top end from INPUT −18 to +18 dB: 0.9.3 → 0.9.4 | −1.5 → −1.0 dB | −0.3 → +0.5 dB | −18.0 → −4.0 dB | −0.1 → +0.2 dB | +0.1 → +1.1 dB | −10.7 → −4.1 dB |
| … largest change per 1 dB of INPUT | 0.10 → 0.08 dB | 0.04 → 0.09 dB | 0.94 → 0.25 dB | 0.09 → 0.09 dB | 0.10 → 0.12 dB | 0.56 → 0.34 dB |
| Loudness of the mix across INPUT with AG on (spread) | 3.5 dB → under 0.1 dB | 8.7 dB → under 0.1 dB | 13.5 dB → under 0.1 dB | 2.5 dB → under 0.1 dB | 6.5 dB → under 0.1 dB | 13.0 dB → under 0.1 dB |
| Treble stage (fitted): bands; detectors; shared pre-emphasis; leak; release | 10 kHz up; peak; α 0.125; leak 0.2; release 20 ms | 5–10 kHz and 10 kHz up; energy 1 ms + 0.375 peak; α 0.375; leak 0; release 60 ms | 5–10 kHz and 10 kHz up; energy 20 ms + 0.5 peak; α 0.5; leak 0; release 20 ms | 10 kHz up; peak; α 0.875; leak 0.35; release 20 ms | 5–10 kHz and 10 kHz up; energy 1 ms + 0.75 peak; α 0.375; leak 0; release 60 ms | 5–10 kHz and 10 kHz up; energy 20 ms + 0.875 peak; α 0.625; leak 0; release 20 ms |
| Pink noise to +13 VU: how far outside the machine (above its total power, or below the part that follows the input), rms | 0.54 dB | 0.85 dB | 1.25 dB | 0.52 dB | 1.10 dB | 1.13 dB |
| … total power above 1.8 kHz per section | −0.0 / −0.3 / −0.4 / −0.3 / −0.2 / +0.1 dB | −0.1 / −0.2 / −0.2 / −0.1 / −0.1 / +1.2 dB | +0.4 / +0.1 / −0.4 / −0.3 / +2.4 / +4.3 dB | −0.1 / −0.1 / −0.2 / −0.0 / −0.2 / +0.2 dB | −0.2 / −0.6 / −1.2 / −1.3 / −0.1 / +1.7 dB | +0.4 / +0.4 / +0.4 / −0.7 / +5.9 / +7.9 dB |
| … sections (VU) | −10 / 0 / +3 / +13 / +23 / +26 | −10 / 0 / +3 / +12 / +22 / +25 | −10 / 0 / +3 / +12 / +22 / +25 | −10 / 0 / +3 / +13 / +23 / +26 | −10 / 0 / +3 / +13 / +23 / +26 | −10 / 0 / +3 / +13 / +23 / +26 |

What the music says. 0.9.3's treble stage read the highs by their peaks: a mix's highs (cymbals,
hats, sibilance) are far spikier than a sine of the same energy, so they read as much louder and
were dulled far too soon. At 3¾ ips the machine's top end comes down about 4 dB from 0 to +19 dB onto
the tape; 0.9.3's came down 16 dB (ATR), and moved up to 0.9 dB for every dB of INPUT. Read by
their energy with a share of their peaks, the highs now come down as the machine's do, and every
model moves smoothly (at most 0.34 dB of top end per dB of INPUT). At the 15 ips speeds the fit
chose peak detectors again, with a leak, and the music error is a little lower than 0.9.3's.
Driven hard at 7½ ips the top end still comes up about a dB against the mids where the machine's
eases off; at 3¾ ips it is a shade darker than the machine's at 0 dB, as asked. Auto gain, made on
the mix, holds its loudness to within 0.1 dB across the whole INPUT range (0.9.3, made on pink
noise: up to 13.5 dB).

**ATR**

| | 15 ips | 7½ ips | 3¾ ips |
|---|---|---|---|
| Hot pass onto the tape / played back | +22.7 / −13.3 dB | +22.4 / −13.1 dB | +22.3 / −12.9 dB |
| Passes agree where they overlap (odd harmonics / compression) | 1.5 / 0.11 dB rms | 1.2 / 0.09 dB rms | 0.7 / 0.23 dB rms |
| Captured range | −30 … +34.7 VU | −30 … +34.4 VU | −30 … +34.2 VU |
| 1 kHz at 0 VU: level; H2, H3 capture → model; waveform match | +0.1 dB; −49.4 → −49.1, −46.5 → −46.5 dB; −70 dB | 0.0 dB; −49.3 → −49.2, −49.6 → −50.3 dB; −66 dB | −0.3 dB; −46.5 → −46.7, −43.2 → −43.6 dB; −58 dB |
| 3.15 kHz at −10 VU: level; H3 capture → model; waveform match | +0.1 dB; −65.0 → −64.6 dB; −73 dB | −0.2 dB; −60.5 → −60.1 dB; −62 dB | +0.2 dB; −48.2 → −46.2 dB; −56 dB |
| 1 kHz at +22 VU (nearly square): level; THD capture → model; waveform match | +0.1 dB; −8.7 → −8.8 dB; −39 dB | +0.1 dB; −12.8 → −13.0 dB; −35 dB | −0.5 dB; −17.0 → −16.7 dB; −34 dB |
| Pink noise, the part that follows the input, 1/3-octave rms (63 Hz–16 kHz, six sections) | 0.66 dB | 1.32 dB | 4.15 dB |
| … mean per section at −10 / +0 / +3 / +13 / +23 / +26 VU | −0.07 / −0.06 / −0.10 / −0.25 / −0.02 / −0.10 dB | +0.02 / −0.01 / −0.01 / +0.04 / +0.61 / +1.42 dB | +0.34 / +1.31 / +1.28 / +3.93 / +4.36 / +4.37 dB |
| … above 1.8 kHz | −0.02 / −0.26 / −0.39 / −0.08 / +0.68 / +0.70 dB | +0.06 / −0.07 / +0.04 / +1.04 / +1.68 / +2.40 dB | +0.87 / +3.38 / +4.00 / +6.91 / +3.96 / +3.37 dB |
| … total power instead, rms; above 1.8 kHz | 0.62 dB; −0.0 / −0.3 / −0.4 / −0.3 / −0.2 / +0.1 | 1.06 dB; −0.1 / −0.2 / −0.2 / −0.1 / −0.1 / +1.2 | 3.01 dB; +0.4 / +0.1 / −0.4 / −0.3 / +2.4 / +4.3 |
| Wow & flutter measured (unweighted 2σ / weighted rms) | 0.031 / 0.010 % | 0.055 / 0.013 % | 0.147 / 0.036 % |

**Maxell** (15 ips from the second normal pass, `V2`: see `docs/CAPTURES.md`)

| | 15 ips | 7½ ips | 3¾ ips |
|---|---|---|---|
| Hot pass onto the tape / played back | +22.9 / −12.7 dB | +22.9 / −12.5 dB | +22.9 / −12.6 dB |
| Passes agree where they overlap (odd harmonics / compression) | 1.9 / 0.06 dB rms | 1.3 / 0.09 dB rms | 1.0 / 0.18 dB rms |
| Captured range | −30 … +34.9 VU | −30 … +34.8 VU | −30 … +34.9 VU |
| 1 kHz at 0 VU: level; H2, H3 capture → model; waveform match | 0.0 dB; −50.0 → −49.9, −53.5 → −54.4 dB; −68 dB | +0.2 dB; −49.8 → −49.8, −54.1 → −54.6 dB; −72 dB | +0.2 dB; −47.1 → −47.8, −45.1 → −44.4 dB; −61 dB |
| 3.15 kHz at −10 VU: level; H3 capture → model; waveform match | +0.2 dB; −69.4 → −68.7 dB; −70 dB | 0.0 dB; −63.0 → −64.9 dB; −60 dB | +2.5 dB; −52.7 → −48.7 dB; −56 dB |
| 1 kHz at +22 VU (nearly square): level; THD capture → model; waveform match | +0.2 dB; −8.1 → −8.0 dB; −37 dB | −0.2 dB; −11.7 → −11.7 dB; −36 dB | −0.4 dB; −16.4 → −16.2 dB; −37 dB |
| Pink noise, the part that follows the input, 1/3-octave rms (63 Hz–16 kHz, six sections) | 0.75 dB | 1.22 dB | 3.90 dB |
| … mean per section at −10 / +0 / +3 / +13 / +23 / +26 VU | −0.06 / −0.07 / −0.15 / −0.46 / −0.44 / −0.48 dB | −0.01 / −0.15 / −0.44 / −0.25 / +0.66 / +1.44 dB | +0.46 / +1.01 / +1.04 / +2.10 / +4.12 / +4.31 dB |
| … above 1.8 kHz | −0.03 / −0.07 / −0.11 / +0.12 / +0.53 / +0.58 dB | −0.08 / −0.49 / −1.06 / −0.35 / +0.90 / +1.99 dB | +0.98 / +3.69 / +4.46 / +4.16 / +5.33 / +5.41 dB |
| … total power instead, rms; above 1.8 kHz | 0.72 dB; −0.1 / −0.1 / −0.2 / −0.0 / −0.2 / +0.2 | 1.29 dB; −0.2 / −0.6 / −1.2 / −1.3 / −0.1 / +1.7 | 3.50 dB; +0.4 / +0.4 / +0.4 / −0.7 / +5.9 / +7.9 |
| Wow & flutter measured (unweighted 2σ / weighted rms) | 0.026 / 0.008 % | 0.055 / 0.012 % | 0.152 / 0.040 % |

What the numbers say. Single tones are reproduced at every level the machine was driven to,
including the hot pass's near-square 1 kHz tone with its droop and overshoot; the treble stage
leaves them exactly as they were, whatever its detectors (`validate_capture.py` gives 0.9.3's
figures to 0.1 dB). (Maxell 3¾ ips is the exception at 3.15 kHz, −10 VU: its sweeps wander ±1 dB
from one level to the next around 2–4 kHz, the tape's contact at that speed, and the model follows
them, so the steady tone comes out 2.5 dB above the capture's.) Pink noise is now a check, not the
fit: up to +13 VU its total treble power is within about 1 dB of the machine's at every speed (at 3¾
ips 0.9.3's was 1–4 dB under it, the part that follows the input matched instead); the hot pass's
pink noise at 3¾ ips (+22 VU and up, past any INPUT setting) is still 2–8 dB brighter. A growing
share of the machine's highs stops following the input as the level goes up (at 3¾ ips about half
the power above 6 kHz on pink noise at 0 VU; on the mix at +19 dB at 7½ ips, 3–6 dB of 10–16 kHz);
the model can't make that part and isn't asked to make up for it. The second harmonic above about
4 kHz is gated conservatively at the slow speeds, where flutter residue raises the noise estimate:
it comes out up to 9 dB lower than measured at −10 VU (7½ ips: −58 → −67 dB at 6.3 kHz).

The C++ engine matches the Python reference on the real models to −110 dB (ATR 3¾ ips −111.2,
Maxell 3¾ ips −109.6, Maxell 15 ips −112.6 dB), the energy detectors included, and runs them at
about 2× real time in stereo at 8× oversampling on one core of the build server (0.9.3's models: the
same within 7 %; the second treble band at 7½ ips costs about that).

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
- The part of the machine's highs that stops following the input as the level rises (its
  distortion products and its own modulation of the highs) is not made: the model gives the part
  that follows, and the fit doesn't let it stand in for the rest with clean treble.
- Driven very hard (the music at +18–19 dB onto the tape), at the slow speeds a mix's 6–10 kHz comes
  out up to 2 dB brighter than the machine's and its top octave (12.5–16 kHz) 2–8 dB darker: there the
  machine's top octave is mostly its own distortion products and modulation, and the upper harmonics
  of the highs are gated conservatively at the slow speeds. The lowest octave (50–60 Hz) comes out 1–2
  dB under, and at ATR 3¾ ips the presence (2–4 kHz) 1.5–3 dB under: the shared detector compresses a
  whole mix by its loudest part, the low end. At 7½ ips this shows as a top end that comes up about a
  dB against the mids at the top of the INPUT range, where the machine's eases off half a dB.
- The treble detectors are a compromise between two kinds of material: a mix's highs saturate on the
  machine about as an energy detector says, pink noise's (continuous, dense highs) more as a peak
  detector says. The fit takes a mix of the two, weighted to music. Pink noise driven far past any INPUT
  setting (the hot pass's, +22 VU and up) still comes out brighter than the machine's.
- Auto gain is calibrated on one mix (the music recordings'): material with a very different balance or
  level drifts by a dB or two at the ends of the INPUT range, as with any static auto gain.
- The leak (a treble band's detector hearing the rest of the signal) fitted at 15 ips (0.2 and 0.35,
  with peak detectors), never at the slower speeds, where it dulled a mix's highs far too fast.
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
