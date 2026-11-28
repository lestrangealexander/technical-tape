// Technical Tape engine — model files (.ttm) and their preparation for a host sample rate
//
// A model is one tape at one speed, captured from the machine (Tools/pipeline). See docs/DESIGN.md and
// Tools/pipeline/ttm.py for the format; this reader follows it exactly.
#pragma once

#include "Fft.h"

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace tt
{
struct ModelData
{
    double sampleRate = 48000.0;
    float refPeak = 0.0f;                // full-scale amplitude of a 0 VU sine's peak
    int numLevels = 0;
    float level0Db = 0.0f, levelStepDb = 1.5f;
    int preDelay = 64;                   // kernel pre-delay at sampleRate
    float lookAheadMs = 1.0f, releaseMs = 60.0f;
    float capturedLowVu = -30.0f, capturedHighVu = 18.0f;
    int flags = 0;
    float margin = 1.02f;

    std::vector<float> pempFreqs, pempDb;          // record pre-emphasis estimate |P|

    struct Branch
    {
        int order = 1;
        int mode = 1;                    // 0: basis mixed by `weights`; 1: one kernel per level
        int count = 0, length = 0;
        std::vector<float> kernels;      // count * length
        std::vector<float> weights;      // numLevels * count (mode 0)
    };
    std::vector<Branch> branches;

    // treble stage (version 3+ models that have one): bands of the fundamental from fcLo up, each indexed
    // by its own detector, which hears the signal through a high-pass at fcLo (and a low-pass at fcHi
    // when another band sits above it), lowest first
    static constexpr int kMaxTrebleBands = 4;
    struct TrebleBand
    {
        float fcLo = 0.0f, fcHi = 0.0f;
        Branch branch;
        float leak = 0.0f;               // how much of the shared detector's pre-emphasised signal it also hears
    };
    std::vector<TrebleBand> treble;
    bool hasTreble() const noexcept { return ! treble.empty(); }

    std::vector<std::pair<float, float>> flutterPeaks;   // Hz, RMS delay (s)
    std::vector<std::pair<float, float>> flutterCont;    // Hz, one-sided PSD of the delay (s^2/Hz)
    float flutterMaxDelay = 0.0f;                         // s

    float agDrive0 = -18.0f, agStep = 1.5f;
    std::vector<float> agGainDb;

    std::string metaJson;
    uint64_t serial = 0;                 // unique per parse: the prepared-model cache's key

    /** Auto-gain compensation (dB) for an input drive (dB). */
    float autoGainDb (float driveDb) const noexcept;
};

/** Parses a .ttm image. Returns an empty string on success, otherwise what's wrong. */
std::string parseModel (const void* data, size_t size, ModelData& out);

//==============================================================================
/** A model made ready for one host sample rate and convolution block size: kernels resampled,
    the pre-emphasis designed, every kernel partitioned and transformed. Immutable once built,
    so any number of engine instances can share it. */
struct PreparedModel
{
    const ModelData* source = nullptr;
    double fs = 48000.0;
    int block = 128;
    int fftSize = 256;

    std::vector<float> pemp;             // minimum-phase FIR of the pre-emphasis estimate
    int lookAhead = 48;
    int preDelay = 64;
    float releaseCoef = 0.0f;
    float margin = 1.02f, refPeak = 0.0f;
    float level0Db = 0.0f, levelStepDb = 1.5f;
    int numLevels = 0;

    struct Branch
    {
        int order = 1, mode = 1, count = 0, partitions = 0;
        AlignedBuffer spectra;           // count * partitions * fftSize
        std::vector<float> weights;      // numLevels * count (mode 0)
        const float* kernelSpectrum (int k, int p) const noexcept
        {
            return spectra.data() + ((size_t) k * (size_t) partitions + (size_t) p) * (size_t) fftSizeOf;
        }
        int fftSizeOf = 256;
    };
    std::vector<Branch> branches;
    int maxOrder = 1;
    int maxPartitions = 0;
    int maxCount = 0;

    struct TrebleBand
    {
        Branch branch;                   // order 1, indexed by this band's detector
        double hpB[3] = { 1.0, 0.0, 0.0 }, hpA[3] = { 1.0, 0.0, 0.0 };
        double lpB[3] = { 1.0, 0.0, 0.0 }, lpA[3] = { 1.0, 0.0, 0.0 };
        bool lowPass = false;
        float leak = 0.0f;
    };
    std::vector<TrebleBand> treble;

    static std::shared_ptr<const PreparedModel> build (const ModelData& m, double fs, int block);
};

/** Convolution block size the engine uses at a sample rate. */
int blockSizeFor (double fs) noexcept;

} // namespace tt
