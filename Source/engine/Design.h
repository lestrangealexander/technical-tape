// Technical Tape engine — filter design done at prepare time (never on the audio thread)
#pragma once

#include <vector>

namespace tt::design
{
/** Minimum-phase FIR following a magnitude curve given at (freqsHz, magDb) (log-frequency interpolation,
    held flat beyond the ends), by the real cepstrum. The same procedure as design_minphase() in
    Tools/pipeline/dsp_shared.py, so the plugin's pre-emphasis matches the one the kernels were built with. */
std::vector<float> minimumPhase (const std::vector<float>& freqsHz, const std::vector<float>& magDb, double fs, int numTaps);

/** Band-limited resampling of an impulse response from fsIn to fsOut, keeping its frequency response
    (and its time structure: a peak at t seconds stays at t seconds). Kaiser-windowed sinc, 64 taps. */
std::vector<float> resampleKernel (const float* in, int length, double fsIn, double fsOut);

/** In-place complex FFT (double), n a power of two; inverse is unscaled. */
void fft (std::vector<double>& re, std::vector<double>& im, bool inverse);

} // namespace tt::design
