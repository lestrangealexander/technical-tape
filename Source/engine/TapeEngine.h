// Technical Tape engine — the whole machine: drive, the captured tape model, wow & flutter, output.
//
// Signal flow (docs/DESIGN.md, section 3):
//   in -> INPUT -> pre-emphasis estimate -> look-ahead peak detector
//                                       -> Chebyshev branches (order n at the lowest alias-free rate)
//                                       -> per-level kernels / basis (partitioned FFT convolution)
//              -> per treble band: band filter -> its own detector -> that band's share of the fundamental
//      -> wow & flutter (fractional delay, measured spectrum) -> OUTPUT x auto gain -> out
// Nothing on the audio thread allocates. Changing tape, speed or oversampling starts a second core that
// warms up and crossfades in, so there is never a click or a gap.
#pragma once

#include "Convolver.h"
#include "Dsp.h"
#include "TapeModel.h"

#include <array>
#include <atomic>
#include <memory>
#include <random>
#include <thread>
#include <vector>

namespace tt
{
//==============================================================================
/** One model at one oversampling factor, for all channels. */
class Core
{
public:
    struct Limits { int maxOrder = 1, maxPartitions = 1, maxCount = 1, maxPemp = 256, block = 128, maxTreble = 0; };

    void allocate (int numChannels, const Limits& lim);
    /** Points the core at a prepared model and an OS factor (1, 2, 4, 8, 16) and resets it. No allocation. */
    void configure (const PreparedModel* pm, int osFactor);
    void reset();
    const PreparedModel* model() const noexcept { return pm; }
    int osFactor() const noexcept { return factor; }

    /** One block (the model's block size) of one channel: x in full-scale units in, y added out. */
    void processBlock (int ch, const float* x, float* y) noexcept;

    /** Latency of the core alone, in samples: OS chain + kernel pre-delay + one block. */
    int latency() const noexcept;

private:
    struct DownChain
    {
        DownStage st[4];
        int nst = 0;
        float pend[4] = {};
        bool has[4] = {};
        void reset() { for (int i = 0; i < 4; ++i) { st[i].reset(); has[i] = false; } }
        /** Pushes one sample at the branch's rate; returns true and sets out when a base sample emerges.
            st[0..3] decimate 16x->8x, 8x->4x, 4x->2x, 2x->1x; a chain at rate 2^nst uses the last nst. */
        bool push (float v, float& out) noexcept
        {
            for (int i = 4 - nst; i < 4; ++i)
            {
                if (! has[i]) { pend[i] = v; has[i] = true; return false; }
                has[i] = false;
                v = st[i].process (pend[i], v);
            }
            out = v;
            return true;
        }
    };

    struct Channel
    {
        Fir pemp;
        Detector det;
        History lraw;
        Delay linDelay, rate1In;
        UpStage up[4];
        Delay pad[5];
        int64_t osCount[5] = {};
        DownChain down[32];
        Delay extra[32];
        BranchConvolver conv[32];
        std::vector<float> bIn[32];
        std::vector<float> lw;
        // treble stage: one of these per band
        struct Treble
        {
            Biquad hp, lp;
            Detector det;
            History lraw;
            BranchConvolver conv;
            std::vector<float> lw;
        };
        Treble treble[ModelData::kMaxTrebleBands];
        int64_t n = 0;
        void reset();
    };

    void handleRate (Channel& c, int ri, float s, int i) noexcept;
    void computeWeights (const PreparedModel::Branch& b, const float* L, const float* Ldb, int& kmin, int& kmax) noexcept;
    void runBranch (BranchConvolver& cv, const PreparedModel::Branch& b, const float* in, const float* L, const float* Ldb) noexcept;

    const PreparedModel* pm = nullptr;
    int factor = 8;
    int B = 128;
    Limits lim;
    std::vector<Channel> chans;
    std::unique_ptr<RealFft> fft;

    // configuration derived in configure()
    int rateOfOrder[32] = {};
    std::vector<int> ordersAtRate[5];
    double dUp[5] = {};
    int maxRateIndex = 0;
    int shiftLw = 0;
    float invRefPeak = 1.0f;
    int numOrders = 1;

    // per-block scratch
    std::vector<float> weights;          // maxCount * B
    std::vector<float> ldb, ldbHf;
    AlignedBuffer yacc;
};

//==============================================================================
/** Transport jitter with the measured spectrum: discrete components (capstan, pinch roller, reels) with
    a slow random wander, plus shaped noise for the continuum. Output in seconds per sample. */
class Flutter
{
public:
    struct Spec
    {
        std::vector<float> freq, amp;    // sinusoids: Hz, peak seconds
        std::vector<float> fir;          // noise shaping at the control rate
        float maxDelay = 0.0f;
    };
    static Spec makeSpec (const ModelData& m);

    void prepare (double fs);
    void setSpec (const Spec* s, uint32_t seed);
    void reset();
    float next() noexcept;

private:
    const Spec* spec = nullptr;
    double fs = 48000.0;
    static constexpr double kControlRate = 1000.0;
    int decim = 48;
    int phaseInControl = 0;
    std::vector<double> phase, wander, wanderVel;
    size_t noisePos = 0;
    float ctl[4] = {};                   // last four control samples (Catmull-Rom)
    uint64_t rng = 0x9E3779B97F4A7C15ull;
    float gauss() noexcept;
    float controlSample() noexcept;
};

/** Fractional delay line (32-tap windowed sinc, 256 phases). */
class FractionalDelay
{
public:
    void prepare (int maxDelaySamples);
    void reset();
    float process (float x, double delaySamples) noexcept;

private:
    std::vector<float> buf;
    size_t mask = 0, w = 0;
    static constexpr int kTaps = 32, kPhases = 256;
    static const std::vector<float>& table();
};

//==============================================================================
struct EngineParams
{
    int model = 0;              // index into the model list given to setModels()
    int osIndex = 3;            // 0..4 = 1x, 2x, 4x, 8x, 16x
    float inputDb = 0.0f;
    float outputDb = 0.0f;
    float wowFlutter = 1.0f;    // 0..2
    bool autoGain = true;
    bool bypass = false;
};

class TapeEngine
{
public:
    TapeEngine() = default;
    ~TapeEngine();
    TapeEngine (const TapeEngine&) = delete;
    TapeEngine& operator= (const TapeEngine&) = delete;

    /** The models this engine can switch between (owned by the caller, must outlive the engine). */
    void setModels (std::vector<const ModelData*> models);
    /** Prepares the model selected in the last setParams() at once; the others are prepared in the
        background the first time they are chosen (playing carries on with the current one meanwhile). */
    void prepare (double fs, int maxBlock, int numChannels);
    void reset();
    int latencySamples() const noexcept { return latency; }

    /** Call from the audio thread before process(). */
    void setParams (const EngineParams& p) noexcept { pending = p; }
    void process (float* const* io, int numChannels, int numSamples) noexcept;

    /** Mean-square levels (full-scale units) of what goes onto the tape (after INPUT) and of the output,
        averaged with VU-meter timing (99 % of a step in 300 ms). Any thread. */
    float inputMeanSquare() const noexcept { return meterIn.load(); }
    float outputMeanSquare() const noexcept { return meterOut.load(); }
    bool isPrepared() const noexcept { return prepared; }

    static constexpr int kOsFactors[5] = { 1, 2, 4, 8, 16 };

private:
    void processBlock() noexcept;
    void startTransition (int model, int osIndex) noexcept;
    void stopWorker();
    void workerLoop();
    const PreparedModel* readyModel (int m) const noexcept { return ready[(size_t) m].load (std::memory_order_acquire); }

    std::vector<const ModelData*> models;
    std::vector<std::shared_ptr<const PreparedModel>> owned;            // written by prepare() and the worker
    std::unique_ptr<std::atomic<const PreparedModel*>[]> readyStore;
    std::atomic<const PreparedModel*>* ready = nullptr;                  // what the audio thread reads
    std::atomic<uint32_t> requests { 0 };                               // models the audio thread wants
    std::atomic<bool> workerRun { false };
    std::thread worker;
    std::vector<Flutter::Spec> flutterSpecs;
    double fs = 48000.0;
    int B = 128, numCh = 2, latency = 0, flutterBase = 0;
    bool prepared = false;

    Core cores[2];
    int active = 0;                     // index of the core that is playing
    enum class Phase { Steady, Warming, Fading } phase = Phase::Steady;
    int phaseCount = 0, warmLen = 0, fadeLen = 0;
    int curModel = -1, curOs = -1, wantModel = 0, wantOs = 3;

    Flutter flutter[2];
    std::vector<FractionalDelay> fdelay;
    std::vector<Delay> dryDelay;

    EngineParams pending, cur;
    float driveGain = 1.0f, outGain = 1.0f, wfAmount = 1.0f, bypassMix = 0.0f;
    float smoothCoef = 0.0f;

    std::vector<std::vector<float>> inFifo, outFifo, blockX, blockY, blockY2;
    std::vector<float> deltaBuf, ramp;
    int fifoPos = 0;
    std::atomic<float> meterIn { 0.0f }, meterOut { 0.0f };
    float vuIn = 0.0f, vuOut = 0.0f, vuCoef = 0.0f;
};

} // namespace tt
