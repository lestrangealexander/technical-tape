// Technical Tape engine — the whole machine
#include "TapeEngine.h"
#include "Design.h"
#include "OsFilters.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <map>
#include <mutex>
#include <tuple>

namespace tt
{
namespace
{
    constexpr double kPi = 3.14159265358979323846;

    /** Rate (as an index into 1, 2, 4, 8, 16) a Chebyshev order runs at: the lowest alias-free one, capped
        by the chosen factor. Same as refengine.branch_rate(). */
    int rateIndexFor (int order, int factor) noexcept
    {
        if (factor <= 1) return 0;
        int r = 2;
        while (2 * r - 1 < order) r *= 2;
        r = std::min (r, factor);
        int ri = 0;
        while ((1 << ri) < r) ++ri;
        return ri;
    }

    /** Up-sampling delay (base samples) of the chain up to rate index ri, without padding. */
    double upDelay (int ri) noexcept
    {
        const int taps[4] = { osf::kStageATaps, osf::kHalfBandTaps[0], osf::kHalfBandTaps[1], osf::kHalfBandTaps[2] };
        double d = 0.0;
        for (int s = 0; s < ri; ++s)
            d += (double) (taps[s] - 1) / 2.0 / (double) (1 << (s + 1));
        return d;
    }

    /** Shared prepared models, so many plugin instances at one rate keep one copy. */
    std::shared_ptr<const PreparedModel> preparedFor (const ModelData* m, double fs, int block)
    {
        static std::mutex mtx;
        static std::map<std::tuple<uint64_t, long long, int>, std::weak_ptr<const PreparedModel>> cache;
        const auto key = std::make_tuple (m->serial, (long long) std::llround (fs * 1000.0), block);
        std::lock_guard<std::mutex> lock (mtx);
        if (auto it = cache.find (key); it != cache.end())
            if (auto sp = it->second.lock())
                return sp;
        auto sp = PreparedModel::build (*m, fs, block);
        cache[key] = sp;
        return sp;
    }
} // namespace

//==============================================================================
void Core::Channel::reset()
{
    pemp.reset();
    det.reset();
    lraw.reset();
    linDelay.reset();
    rate1In.reset();
    for (auto& u : up) u.reset();
    for (auto& p : pad) p.reset();
    for (auto& c : osCount) c = 0;
    for (auto& d : down) d.reset();
    for (auto& e : extra) e.reset();
    for (auto& cv : conv) cv.reset();
    for (auto& b : bIn) std::fill (b.begin(), b.end(), 0.0f);
    std::fill (lw.begin(), lw.end(), 0.0f);
    for (auto& t : treble)
    {
        t.hp.reset();
        t.lp.reset();
        t.hpRest.reset();
        t.lpRest.reset();
        t.det.reset();
        t.lraw.reset();
        t.conv.reset();
        std::fill (t.lw.begin(), t.lw.end(), 0.0f);
    }
    n = 0;
}

void Core::allocate (int numChannels, const Limits& l)
{
    lim = l;
    B = l.block;
    fft = std::make_unique<RealFft> (2 * B);
    chans.clear();
    chans.resize ((size_t) numChannels);
    for (auto& c : chans)
    {
        std::vector<float> reserveTaps ((size_t) l.maxPemp, 0.0f);
        c.pemp.setTaps (reserveTaps);
        c.det.prepare (80, 0.9f);
        c.lraw.prepare (2048);
        c.linDelay.setLength (osf::kLatency);
        c.rate1In.setLength (80);
        c.up[0].setTaps (osf::kStageA, osf::kStageATaps);
        c.up[1].setTaps (osf::kHalfBand0, osf::kHalfBandTaps[0]);
        c.up[2].setTaps (osf::kHalfBand1, osf::kHalfBandTaps[1]);
        c.up[3].setTaps (osf::kHalfBand2, osf::kHalfBandTaps[2]);
        for (auto& p : c.pad) p.setLength (16);
        for (int o = 0; o < 32; ++o)
        {
            // stages from 16x down to 1x; a chain at rate 2^k uses the last k of them
            c.down[o].st[0].setTaps (osf::kHalfBand2, osf::kHalfBandTaps[2]);
            c.down[o].st[1].setTaps (osf::kHalfBand1, osf::kHalfBandTaps[1]);
            c.down[o].st[2].setTaps (osf::kHalfBand0, osf::kHalfBandTaps[0]);
            c.down[o].st[3].setTaps (osf::kStageA, osf::kStageATaps);
            c.extra[o].setLength (osf::kLatency);
            if (o >= 1 && o <= l.maxOrder)
            {
                c.conv[o].prepare (B, l.maxPartitions);
                c.bIn[o].assign ((size_t) B, 0.0f);
            }
        }
        c.lw.assign ((size_t) B, 0.0f);
        for (int k = 0; k < l.maxTreble; ++k)
        {
            auto& t = c.treble[k];
            t.det.prepare (80, 0.9f);
            t.lraw.prepare (2048);
            t.conv.prepare (B, l.maxPartitions);
            t.lw.assign ((size_t) B, 0.0f);
        }
    }
    for (auto& o : ordersAtRate) o.reserve (32);
    weights.assign ((size_t) l.maxCount * (size_t) B, 0.0f);
    ldb.assign ((size_t) B, 0.0f);
    ldbHf.assign ((size_t) B, 0.0f);
    yacc.allocate ((size_t) B);
}

void Core::configure (const PreparedModel* model, int osFactor)
{
    pm = model;
    factor = osFactor;
    for (auto& o : ordersAtRate) o.clear();
    maxRateIndex = 0;
    numOrders = 0;
    for (const auto& b : pm->branches) numOrders = std::max (numOrders, b.order);
    for (int n = 2; n <= numOrders; ++n)
    {
        const int ri = rateIndexFor (n, factor);
        rateOfOrder[n] = ri;
        ordersAtRate[ri].push_back (n);          // capacity grows once, on the first configure
        maxRateIndex = std::max (maxRateIndex, ri);
    }
    for (int ri = 1; ri < 5; ++ri)
        dUp[ri] = upDelay (ri) + (double) osf::kPadAtRate[ri] / (double) (1 << ri);
    shiftLw = osf::kLatency + pm->preDelay - (pm->lookAhead - 1);
    invRefPeak = 1.0f / pm->refPeak;
    for (auto& c : chans)
    {
        c.pemp.setTaps (pm->pemp);
        c.det.prepare (pm->lookAhead, pm->releaseCoef);
        c.linDelay.setLength (osf::kLatency);
        c.rate1In.setLength (pm->lookAhead - 1);
        for (int ri = 1; ri < 5; ++ri) c.pad[ri].setLength (osf::kPadAtRate[ri]);
        for (int n = 2; n <= numOrders; ++n)
        {
            const int ri = rateOfOrder[n];
            const int chain = ri == 0 ? pm->lookAhead - 1 : osf::kChainDelay[ri];
            c.extra[n].setLength (osf::kLatency - chain);
            c.down[n].nst = ri;
        }
        for (const auto& b : pm->branches)
            c.conv[b.order].setBranch (&b);
        for (size_t k = 0; k < pm->treble.size(); ++k)
        {
            const auto& band = pm->treble[k];
            auto& t = c.treble[k];
            t.hp.setCoefficients (band.hpB, band.hpA);
            t.lp.setCoefficients (band.lpB, band.lpA);
            t.hpRest.setCoefficients (band.hpB, band.hpA);
            t.lpRest.setCoefficients (band.lpB, band.lpA);
            t.det.prepare (pm->lookAhead, pm->releaseCoef);
            t.conv.setBranch (&band.branch);
        }
    }
    reset();
}

void Core::reset()
{
    for (auto& c : chans) c.reset();
}

int Core::latency() const noexcept
{
    return osf::kLatency + (pm != nullptr ? pm->preDelay : 0) + B;
}

void Core::handleRate (Channel& c, int ri, float s, int i) noexcept
{
    const float sp = c.pad[ri].process (s);
    const double tau = (double) c.osCount[ri] / (double) (1 << ri) - dUp[ri];
    ++c.osCount[ri];
    const auto& orders = ordersAtRate[ri];
    if (orders.empty()) return;
    const float Lin = pm->margin * c.lraw.at (tau + (double) (pm->lookAhead - 1));
    const float u = std::clamp (sp / Lin, -1.0f, 1.0f);
    float tm = 1.0f, tc = u;
    const int top = orders.back();
    for (int n = 2; n <= top; ++n)
    {
        const float tn = 2.0f * u * tc - tm;
        tm = tc;
        tc = tn;
        if (rateOfOrder[n] == ri)
        {
            float out;
            if (c.down[n].push (Lin * tn, out))
                c.bIn[n][(size_t) i] = c.extra[n].process (out);
        }
    }
}

void Core::computeWeights (const PreparedModel::Branch& b, const float* L, const float* ldbIn, int& kmin, int& kmax) noexcept
{
    const float* lev = ldbIn;            // the levels in dB this branch is indexed by
    const int nl = pm->numLevels;
    const int n = b.order;
    const float L0 = std::pow (10.0f, pm->level0Db / 20.0f);
    float* w = weights.data();
    if (b.mode == 0)
    {
        const int R = b.count;
        kmin = 0;
        kmax = R - 1;
        for (int i = 0; i < B; ++i)
        {
            const float idx = (lev[(size_t) i] - pm->level0Db) / pm->levelStepDb;
            if (idx < 0.0f)
            {
                const float s = n > 1 ? std::pow (L[i] / L0, (float) (n - 1)) : 1.0f;
                for (int r = 0; r < R; ++r) w[(size_t) r * (size_t) B + (size_t) i] = b.weights[(size_t) r] * s;
                continue;
            }
            if (idx > (float) (nl - 1))
            {
                for (int r = 0; r < R; ++r) w[(size_t) r * (size_t) B + (size_t) i] = b.weights[(size_t) (nl - 1) * (size_t) R + (size_t) r];
                continue;
            }
            const int i0 = std::clamp ((int) std::floor (idx), 0, nl - 2);
            const float t = std::clamp (idx - (float) i0, 0.0f, 1.0f);
            const float t2 = t * t, t3 = t2 * t;
            const float c0 = 0.5f * (-t + 2.0f * t2 - t3), c1 = 0.5f * (2.0f - 5.0f * t2 + 3.0f * t3);
            const float c2 = 0.5f * (t + 4.0f * t2 - 3.0f * t3), c3 = 0.5f * (-t2 + t3);
            const int j0 = std::max (i0 - 1, 0), j2 = std::min (i0 + 1, nl - 1), j3 = std::min (i0 + 2, nl - 1);
            for (int r = 0; r < R; ++r)
                w[(size_t) r * (size_t) B + (size_t) i] = c0 * b.weights[(size_t) j0 * (size_t) R + (size_t) r]
                                                        + c1 * b.weights[(size_t) i0 * (size_t) R + (size_t) r]
                                                        + c2 * b.weights[(size_t) j2 * (size_t) R + (size_t) r]
                                                        + c3 * b.weights[(size_t) j3 * (size_t) R + (size_t) r];
        }
        return;
    }
    // mode 1: one kernel per level; only the neighbours of each sample's level are non-zero
    kmin = nl;
    kmax = -1;
    for (int i = 0; i < B; ++i)
    {
        const float idx = (lev[(size_t) i] - pm->level0Db) / pm->levelStepDb;
        int lo, hi;
        if (idx < 0.0f) { lo = hi = 0; }
        else if (idx > (float) (nl - 1)) { lo = hi = nl - 1; }
        else
        {
            const int i0 = std::clamp ((int) std::floor (idx), 0, nl - 2);
            lo = std::max (i0 - 1, 0);
            hi = std::min (i0 + 2, nl - 1);
        }
        kmin = std::min (kmin, lo);
        kmax = std::max (kmax, hi);
    }
    for (int k = kmin; k <= kmax; ++k)
        std::fill_n (w + (size_t) k * (size_t) B, B, 0.0f);
    for (int i = 0; i < B; ++i)
    {
        const float idx = (lev[(size_t) i] - pm->level0Db) / pm->levelStepDb;
        if (idx < 0.0f)
        {
            w[(size_t) i] += n > 1 ? std::pow (L[i] / L0, (float) (n - 1)) : 1.0f;
            continue;
        }
        if (idx > (float) (nl - 1))
        {
            w[(size_t) (nl - 1) * (size_t) B + (size_t) i] += 1.0f;
            continue;
        }
        const int i0 = std::clamp ((int) std::floor (idx), 0, nl - 2);
        const float t = std::clamp (idx - (float) i0, 0.0f, 1.0f);
        const float t2 = t * t, t3 = t2 * t;
        const float cc[4] = { 0.5f * (-t + 2.0f * t2 - t3), 0.5f * (2.0f - 5.0f * t2 + 3.0f * t3),
                              0.5f * (t + 4.0f * t2 - 3.0f * t3), 0.5f * (-t2 + t3) };
        const int js[4] = { std::max (i0 - 1, 0), i0, std::min (i0 + 1, nl - 1), std::min (i0 + 2, nl - 1) };
        for (int j = 0; j < 4; ++j)
            w[(size_t) js[j] * (size_t) B + (size_t) i] += cc[j];
    }
}

void Core::processBlock (int ch, const float* x, float* y) noexcept
{
    if (pm == nullptr) return;
    auto& c = chans[(size_t) ch];
    const float mu = pm->margin;
    for (int i = 0; i < B; ++i)
    {
        const float vp = c.pemp.process (x[i] * invRefPeak);
        const float L = c.det.process (vp);
        c.lraw.push (L);
        c.bIn[1][(size_t) i] = c.linDelay.process (vp);
        if (! ordersAtRate[0].empty())
        {
            const float vin = c.rate1In.process (vp);
            const float Lin = mu * L;
            const float u = std::clamp (vin / Lin, -1.0f, 1.0f);
            float tm = 1.0f, tc = u;
            for (int n = 2; n <= ordersAtRate[0].back(); ++n)
            {
                const float tn = 2.0f * u * tc - tm;
                tm = tc;
                tc = tn;
                if (rateOfOrder[n] == 0) c.bIn[n][(size_t) i] = c.extra[n].process (Lin * tn);
            }
        }
        if (maxRateIndex >= 1)
        {
            float s2[2];
            c.up[0].process (vp, s2);
            for (float a : s2)
            {
                handleRate (c, 1, a, i);
                if (maxRateIndex < 2) continue;
                float s4[2];
                c.up[1].process (a, s4);
                for (float b4 : s4)
                {
                    handleRate (c, 2, b4, i);
                    if (maxRateIndex < 3) continue;
                    float s8[2];
                    c.up[2].process (b4, s8);
                    for (float b8 : s8)
                    {
                        handleRate (c, 3, b8, i);
                        if (maxRateIndex < 4) continue;
                        float s16[2];
                        c.up[3].process (b8, s16);
                        handleRate (c, 4, s16[0], i);
                        handleRate (c, 4, s16[1], i);
                    }
                }
            }
        }
        c.lw[(size_t) i] = c.lraw.at ((double) (c.n - shiftLw));
        for (size_t k = 0; k < pm->treble.size(); ++k)
        {
            // each treble band's detector: the same look-ahead detector, hearing that band, plus some of the
            // rest of what the shared detector hears (a loud low end and midrange squash the highs too)
            auto& t = c.treble[k];
            const auto& band = pm->treble[k];
            float v = t.hp.process (x[i] * invRefPeak);
            if (band.lowPass) v = t.lp.process (v);
            if (band.leak > 0.0f)
            {
                float r = t.hpRest.process (vp);
                if (band.lowPass) r = t.lpRest.process (r);
                v += band.leak * (vp - r);
            }
            t.lraw.push (t.det.process (v));
            t.lw[(size_t) i] = t.lraw.at ((double) (c.n - shiftLw));
        }
        ++c.n;
    }
    for (int i = 0; i < B; ++i)
        ldb[(size_t) i] = 20.0f * std::log10 (c.lw[(size_t) i]);
    yacc.clear();
    for (const auto& b : pm->branches)
        runBranch (c.conv[b.order], b, c.bIn[b.order].data(), c.lw.data(), ldb.data());
    for (size_t k = 0; k < pm->treble.size(); ++k)
    {
        auto& t = c.treble[k];
        for (int i = 0; i < B; ++i)
            ldbHf[(size_t) i] = 20.0f * std::log10 (t.lw[(size_t) i]);
        runBranch (t.conv, pm->treble[k].branch, c.bIn[1].data(), t.lw.data(), ldbHf.data());
    }
    const float g = pm->refPeak;
    for (int i = 0; i < B; ++i) y[i] += yacc[(size_t) i] * g;
}

void Core::runBranch (BranchConvolver& cv, const PreparedModel::Branch& b, const float* in, const float* L,
                      const float* Ldb) noexcept
{
    cv.pushBlock (in, *fft);
    int kmin = 0, kmax = -1;
    computeWeights (b, L, Ldb, kmin, kmax);
    for (int k = kmin; k <= kmax; ++k)
    {
        const float* w = weights.data() + (size_t) k * (size_t) B;
        bool any = false;
        for (int i = 0; i < B && ! any; ++i) any = w[i] != 0.0f;
        if (any) cv.addKernel (k, w, yacc.data(), *fft);
    }
}

//==============================================================================
Flutter::Spec Flutter::makeSpec (const ModelData& m)
{
    Spec s;
    for (const auto& [f, rms] : m.flutterPeaks)
        if (f > 0.0f && rms > 0.0f)
        {
            s.freq.push_back (f);
            s.amp.push_back (rms * std::sqrt (2.0f));
        }
    // continuum: a 65.5 s loop of noise at the control rate with the measured spectrum (random phases)
    const size_t N = 1 << 16;
    std::vector<double> re (N, 0.0), im (N, 0.0);
    uint64_t st = 0x243F6A8885A308D3ull;
    auto uni = [&st]
    {
        st ^= st << 13; st ^= st >> 7; st ^= st << 17;
        return (double) (st >> 11) * (1.0 / 9007199254740992.0);
    };
    double var = 0.0;
    for (size_t k = 1; k < N / 2; ++k)
    {
        const double f = (double) k * kControlRate / (double) N;
        double psd = 0.0;
        if (! m.flutterCont.empty() && f >= 0.05 && f <= 300.0)
        {
            const auto& c = m.flutterCont;
            const double lf = std::log (std::max (f, (double) c.front().first));
            size_t j = 1;
            while (j < c.size() - 1 && std::log ((double) c[j].first) < lf) ++j;
            const double x0 = std::log ((double) c[j - 1].first), x1 = std::log ((double) c[j].first);
            const double y0 = std::log (std::max (1e-30, (double) c[j - 1].second)), y1 = std::log (std::max (1e-30, (double) c[j].second));
            const double t = std::clamp ((lf - x0) / (x1 - x0), 0.0, 1.0);
            psd = std::exp (y0 + t * (y1 - y0));
            if (f < 0.2) psd *= (f - 0.05) / 0.15;                       // fade out below the measured band
            if (f > (double) c.back().first) psd *= std::max (0.0, (300.0 - f) / (300.0 - (double) c.back().first));
        }
        var += psd * kControlRate / (double) N;
        const double a = std::sqrt (psd);
        const double ph = 2.0 * kPi * uni();
        re[k] = a * std::cos (ph);
        im[k] = a * std::sin (ph);
        re[N - k] = re[k];
        im[N - k] = -im[k];
    }
    design::fft (re, im, true);
    double ms = 0.0;
    for (size_t i = 0; i < N; ++i) ms += re[i] * re[i];
    ms /= (double) N;
    const double scale = ms > 0.0 ? std::sqrt (var / ms) : 0.0;
    s.fir.resize (N);
    for (size_t i = 0; i < N; ++i) s.fir[i] = (float) (re[i] * scale);
    s.maxDelay = m.flutterMaxDelay;
    return s;
}

void Flutter::prepare (double sampleRate)
{
    fs = sampleRate;
    decim = std::max (1, (int) std::lround (fs / kControlRate));
    phase.assign (16, 0.0);
    wander.assign (16, 0.0);
    wanderVel.assign (16, 0.0);
}

void Flutter::setSpec (const Spec* s, uint32_t seed)
{
    spec = s;
    rng = 0x9E3779B97F4A7C15ull ^ ((uint64_t) seed * 0xBF58476D1CE4E5B9ull);
    reset();
}

void Flutter::reset()
{
    for (size_t k = 0; k < phase.size(); ++k)
    {
        phase[k] = 2.0 * kPi * (double) ((rng >> (k * 3)) & 1023) / 1024.0;
        wander[k] = 0.0;
        wanderVel[k] = 0.0;
    }
    noisePos = (size_t) (rng & 0xFFFF);
    phaseInControl = 0;
    for (float& c : ctl) c = 0.0f;
    if (spec != nullptr)
        for (int k = 0; k < 4; ++k) ctl[k] = controlSample();
}

float Flutter::gauss() noexcept
{
    rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
    const double u1 = ((double) (rng >> 11) + 0.5) * (1.0 / 9007199254740992.0);
    rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
    const double u2 = ((double) (rng >> 11) + 0.5) * (1.0 / 9007199254740992.0);
    return (float) (std::sqrt (-2.0 * std::log (u1)) * std::cos (2.0 * kPi * u2));
}

float Flutter::controlSample() noexcept
{
    if (spec == nullptr || spec->fir.empty()) return 0.0f;
    const float v = spec->fir[noisePos];
    noisePos = noisePos + 1 == spec->fir.size() ? 0 : noisePos + 1;
    // slow random wander of each component's frequency (+/- ~0.5 %), so it is never perfectly periodic
    for (size_t k = 0; k < spec->freq.size() && k < wander.size(); ++k)
    {
        wanderVel[k] = 0.995 * wanderVel[k] + 0.0005 * gauss();
        wander[k] = std::clamp (0.999 * wander[k] + wanderVel[k], -1.0, 1.0);
    }
    return v;
}

float Flutter::next() noexcept
{
    if (spec == nullptr) return 0.0f;
    // continuum: Catmull-Rom between control samples
    const float t = (float) phaseInControl / (float) decim;
    const float p0 = ctl[0], p1 = ctl[1], p2 = ctl[2], p3 = ctl[3];
    const float noise = 0.5f * ((2.0f * p1) + (-p0 + p2) * t + (2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3) * t * t
                                + (-p0 + 3.0f * p1 - 3.0f * p2 + p3) * t * t * t);
    if (++phaseInControl >= decim)
    {
        phaseInControl = 0;
        ctl[0] = ctl[1]; ctl[1] = ctl[2]; ctl[2] = ctl[3];
        ctl[3] = controlSample();
    }
    double tone = 0.0;
    for (size_t k = 0; k < spec->freq.size() && k < phase.size(); ++k)
    {
        tone += (double) spec->amp[k] * std::sin (phase[k]);
        phase[k] += 2.0 * kPi * (double) spec->freq[k] * (1.0 + 0.005 * wander[k]) / fs;
        if (phase[k] > 2.0 * kPi) phase[k] -= 2.0 * kPi;
    }
    return noise + (float) tone;
}

//==============================================================================
const std::vector<float>& FractionalDelay::table()
{
    static const std::vector<float> tab = []
    {
        std::vector<float> t ((size_t) (kPhases + 1) * kTaps);
        auto i0 = [] (double x)
        {
            double s = 1.0, term = 1.0;
            for (int k = 1; k < 50; ++k) { term *= (x * x / 4.0) / ((double) k * k); s += term; }
            return s;
        };
        const double beta = 8.0, half = kTaps / 2 + 0.5;
        for (int p = 0; p <= kPhases; ++p)
        {
            const double mu = (double) p / (double) kPhases;
            double sum = 0.0;
            std::vector<double> c ((size_t) kTaps);
            for (int j = 0; j < kTaps; ++j)
            {
                const double x = (double) (j - kTaps / 2 + 1) - mu;       // taps k = -15..16
                const double s = std::abs (x) < 1e-12 ? 1.0 : std::sin (kPi * x) / (kPi * x);
                const double r = x / half;
                const double w = std::abs (r) < 1.0 ? i0 (beta * std::sqrt (1.0 - r * r)) / i0 (beta) : 0.0;
                c[(size_t) j] = s * w;
                sum += c[(size_t) j];
            }
            for (int j = 0; j < kTaps; ++j) t[(size_t) p * kTaps + (size_t) j] = (float) (c[(size_t) j] / sum);
        }
        return t;
    }();
    return tab;
}

void FractionalDelay::prepare (int maxDelaySamples)
{
    size_t n = 1;
    while (n < (size_t) maxDelaySamples + kTaps + 8) n <<= 1;
    buf.assign (n, 0.0f);
    mask = n - 1;
    w = 0;
    (void) table();
}

void FractionalDelay::reset()
{
    std::fill (buf.begin(), buf.end(), 0.0f);
    w = 0;
}

float FractionalDelay::process (float x, double d) noexcept
{
    buf[w] = x;
    d = std::clamp (d, (double) (kTaps / 2), (double) (mask + 1 - kTaps / 2 - 2));
    const double di = std::floor (d);
    const int n = (int) di;
    const double mu = d - di;
    const double pp = mu * kPhases;
    const int p = std::min ((int) pp, kPhases - 1);
    const float f = (float) (pp - p);
    const float* c0 = table().data() + (size_t) p * kTaps;
    const float* c1 = c0 + kTaps;
    float acc = 0.0f;
    for (int j = 0; j < kTaps; ++j)
    {
        const int k = j - kTaps / 2 + 1;                                  // -15..16
        const float s = buf[(w - (size_t) (n + k)) & mask];
        acc += s * (c0[j] + f * (c1[j] - c0[j]));
    }
    w = (w + 1) & mask;
    return acc;
}

//==============================================================================
TapeEngine::~TapeEngine()
{
    stopWorker();
}

void TapeEngine::setModels (std::vector<const ModelData*> m)
{
    stopWorker();
    models = std::move (m);
}

void TapeEngine::stopWorker()
{
    workerRun.store (false);
    if (worker.joinable()) worker.join();
}

void TapeEngine::workerLoop()
{
    while (workerRun.load())
    {
        const uint32_t want = requests.load();
        for (size_t m = 0; m < models.size() && workerRun.load(); ++m)
            if ((want & (1u << m)) != 0 && ready[m].load() == nullptr)
            {
                owned[m] = preparedFor (models[m], fs, B);
                ready[m].store (owned[m].get(), std::memory_order_release);
            }
        std::this_thread::sleep_for (std::chrono::milliseconds (5));
    }
}

void TapeEngine::prepare (double sampleRate, int, int channels)
{
    stopWorker();
    prepared = false;
    fs = sampleRate;
    numCh = std::max (1, channels);
    B = blockSizeFor (fs);
    flutterSpecs.clear();
    owned.assign (models.size(), nullptr);
    readyStore = std::make_unique<std::atomic<const PreparedModel*>[]> (std::max<size_t> (1, models.size()));
    ready = readyStore.get();
    for (size_t m = 0; m < models.size(); ++m) ready[m].store (nullptr);
    requests.store (0);
    if (models.empty()) return;

    // room for the largest model, worked out without preparing them
    Core::Limits lim;
    lim.block = B;
    float maxDelay = 20e-6f;
    for (const auto* m : models)
    {
        auto room = [&] (const ModelData::Branch& b)
        {
            const int len = (int) std::ceil ((double) b.length * fs / m->sampleRate);
            lim.maxOrder = std::max (lim.maxOrder, b.order);
            lim.maxPartitions = std::max (lim.maxPartitions, (len + B - 1) / B);
            lim.maxCount = std::max (lim.maxCount, b.count);
        };
        for (const auto& b : m->branches) room (b);
        for (const auto& band : m->treble) room (band.branch);
        lim.maxTreble = std::max (lim.maxTreble, (int) m->treble.size());
        lim.maxPemp = std::max (lim.maxPemp, (int) std::lround (256.0 * fs / 48000.0));
        flutterSpecs.push_back (Flutter::makeSpec (*m));
        maxDelay = std::max (maxDelay, m->flutterMaxDelay);
    }
    curModel = std::clamp (pending.model, 0, (int) models.size() - 1);
    owned[(size_t) curModel] = preparedFor (models[(size_t) curModel], fs, B);
    ready[curModel].store (owned[(size_t) curModel].get());

    for (auto& c : cores) c.allocate (numCh, lim);
    // mean delay of the transport: room for the largest excursion (99.99th percentile of the measured delay,
    // x2 for the W&F amount, x2.5 margin) on both sides, plus half the interpolator's length
    flutterBase = std::clamp ((int) std::ceil (fs * 5.0 * (double) maxDelay) + 20, 24, (int) (0.02 * fs));
    for (auto& f : flutter) f.prepare (fs);
    fdelay.assign ((size_t) numCh, {});
    for (auto& d : fdelay) d.prepare (2 * flutterBase + 64);
    const int coreLat = osf::kLatency + owned[(size_t) curModel]->preDelay;
    dryDelay.assign ((size_t) numCh, {});
    for (auto& d : dryDelay) d.setLength (coreLat + flutterBase);
    latency = B + coreLat + flutterBase;
    smoothCoef = (float) (1.0 - std::exp (-1.0 / (0.02 * fs)));
    vuCoef = (float) (1.0 - std::exp (-(double) B / (0.065 * fs)));       // tau 65 ms: 99 % in 300 ms
    vuIn = vuOut = 0.0f;
    // the incoming core runs silently for its longest kernel (its convolution state is then complete)
    warmLen = std::max (lim.maxPartitions * B, ((int) std::lround (4096.0 * fs / 48000.0) + B - 1) / B * B);
    fadeLen = ((int) std::lround (2048.0 * fs / 48000.0) + B - 1) / B * B;
    inFifo.assign ((size_t) numCh, std::vector<float> ((size_t) B, 0.0f));
    outFifo.assign ((size_t) numCh, std::vector<float> ((size_t) B, 0.0f));
    blockX.assign ((size_t) numCh, std::vector<float> ((size_t) B, 0.0f));
    blockY.assign ((size_t) numCh, std::vector<float> ((size_t) B, 0.0f));
    blockY2.assign ((size_t) numCh, std::vector<float> ((size_t) B, 0.0f));
    deltaBuf.assign ((size_t) B, 0.0f);
    ramp.assign ((size_t) B, 0.0f);
    cur = pending;
    curOs = std::clamp (pending.osIndex, 0, 4);
    wantModel = curModel;
    wantOs = curOs;
    active = 0;
    cores[0].configure (owned[(size_t) curModel].get(), kOsFactors[curOs]);
    flutter[0].setSpec (&flutterSpecs[(size_t) curModel], 1);
    phase = Phase::Steady;
    driveGain = std::pow (10.0f, pending.inputDb / 20.0f);
    const float ag = pending.autoGain ? models[(size_t) curModel]->autoGainDb (pending.inputDb) : 0.0f;
    outGain = std::pow (10.0f, (pending.outputDb + ag) / 20.0f);
    wfAmount = pending.wowFlutter;
    bypassMix = pending.bypass ? 1.0f : 0.0f;
    fifoPos = 0;
    prepared = true;
    reset();
    workerRun.store (true);
    worker = std::thread ([this] { workerLoop(); });
}

void TapeEngine::reset()
{
    if (! prepared) return;
    cores[active].reset();
    for (auto& d : fdelay) d.reset();
    for (auto& d : dryDelay) d.reset();
    for (auto& f : inFifo) std::fill (f.begin(), f.end(), 0.0f);
    for (auto& f : outFifo) std::fill (f.begin(), f.end(), 0.0f);
    flutter[active].reset();
    fifoPos = 0;
}

void TapeEngine::startTransition (int model, int osIndex) noexcept
{
    const int other = 1 - active;
    cores[other].configure (readyModel (model), kOsFactors[osIndex]);
    flutter[other].setSpec (&flutterSpecs[(size_t) model], (uint32_t) (model * 7 + osIndex + 3));
    curModel = model;
    curOs = osIndex;
    phase = Phase::Warming;
    phaseCount = 0;
}

void TapeEngine::process (float* const* io, int channels, int numSamples) noexcept
{
    if (! prepared) return;
    const int nc = std::min (channels, numCh);
    int pos = 0;
    while (pos < numSamples)
    {
        const int m = std::min (numSamples - pos, B - fifoPos);
        for (int c = 0; c < nc; ++c)
        {
            std::memcpy (inFifo[(size_t) c].data() + fifoPos, io[c] + pos, (size_t) m * sizeof (float));
            std::memcpy (io[c] + pos, outFifo[(size_t) c].data() + fifoPos, (size_t) m * sizeof (float));
        }
        for (int c = nc; c < numCh; ++c)          // mono host buffer on a stereo-prepared engine
            std::memcpy (inFifo[(size_t) c].data() + fifoPos, io[0] + pos, (size_t) m * sizeof (float));
        fifoPos += m;
        pos += m;
        if (fifoPos == B)
        {
            processBlock();
            fifoPos = 0;
        }
    }
}

void TapeEngine::processBlock() noexcept
{
    cur = pending;
    const int wantM = std::clamp (cur.model, 0, (int) models.size() - 1);
    const int wantO = std::clamp (cur.osIndex, 0, 4);
    if (phase == Phase::Steady && (wantM != curModel || wantO != curOs))
    {
        if (readyModel (wantM) != nullptr)
            startTransition (wantM, wantO);
        else
            requests.fetch_or (1u << wantM);           // the worker prepares it; carry on meanwhile
    }

    const ModelData* modelNow = models[(size_t) curModel];
    const float driveTarget = std::pow (10.0f, cur.inputDb / 20.0f);
    const float ag = cur.autoGain ? modelNow->autoGainDb (cur.inputDb) : 0.0f;
    const float outTarget = std::pow (10.0f, (cur.outputDb + ag) / 20.0f);
    const float wfTarget = std::clamp (cur.wowFlutter, 0.0f, 2.0f);
    const float bypassTarget = cur.bypass ? 1.0f : 0.0f;

    // per-sample gain ramps, shared by all channels
    float dg = driveGain;
    for (int i = 0; i < B; ++i)
    {
        dg += smoothCoef * (driveTarget - dg);
        ramp[(size_t) i] = dg;
    }
    double msIn = 0.0;
    for (int c = 0; c < numCh; ++c)
        for (int i = 0; i < B; ++i)
        {
            const float v = inFifo[(size_t) c][(size_t) i] * ramp[(size_t) i];
            blockX[(size_t) c][(size_t) i] = v;
            msIn += (double) v * (double) v;
        }
    driveGain = dg;
    vuIn += vuCoef * ((float) (msIn / (double) (B * numCh)) - vuIn);
    meterIn.store (vuIn);

    const int other = 1 - active;
    const bool two = phase != Phase::Steady;
    for (int c = 0; c < numCh; ++c)
    {
        std::fill (blockY[(size_t) c].begin(), blockY[(size_t) c].end(), 0.0f);
        cores[active].processBlock (c, blockX[(size_t) c].data(), blockY[(size_t) c].data());
        if (two)
        {
            std::fill (blockY2[(size_t) c].begin(), blockY2[(size_t) c].end(), 0.0f);
            cores[other].processBlock (c, blockX[(size_t) c].data(), blockY2[(size_t) c].data());
        }
    }
    // wow & flutter (one transport for all channels), crossfading the generators with the cores
    for (int i = 0; i < B; ++i)
    {
        float d = flutter[active].next();
        if (two)
        {
            const float dn = flutter[other].next();
            if (phase == Phase::Fading)
            {
                const float a = (float) (phaseCount + i) / (float) fadeLen;
                d += a * (dn - d);
            }
        }
        wfAmount += smoothCoef * (wfTarget - wfAmount);
        deltaBuf[(size_t) i] = (float) (flutterBase + (double) d * fs * wfAmount);
    }
    if (two && phase == Phase::Fading)
        for (int c = 0; c < numCh; ++c)
            for (int i = 0; i < B; ++i)
            {
                const float a = (float) (phaseCount + i) / (float) fadeLen;
                blockY[(size_t) c][(size_t) i] += a * (blockY2[(size_t) c][(size_t) i] - blockY[(size_t) c][(size_t) i]);
            }

    double msOut = 0.0;
    float og = outGain, bm = bypassMix;
    for (int c = 0; c < numCh; ++c)
    {
        og = outGain;
        bm = bypassMix;
        auto& fd = fdelay[(size_t) c];
        auto& dd = dryDelay[(size_t) c];
        for (int i = 0; i < B; ++i)
        {
            og += smoothCoef * (outTarget - og);
            bm += smoothCoef * (bypassTarget - bm);
            const float wet = fd.process (blockY[(size_t) c][(size_t) i], (double) deltaBuf[(size_t) i]) * og;
            const float dry = dd.process (inFifo[(size_t) c][(size_t) i]);
            const float o = wet + bm * (dry - wet);
            outFifo[(size_t) c][(size_t) i] = o;
            msOut += (double) o * (double) o;
        }
    }
    outGain = og;
    bypassMix = bm;
    vuOut += vuCoef * ((float) (msOut / (double) (B * numCh)) - vuOut);
    meterOut.store (vuOut);

    if (two)
    {
        phaseCount += B;
        if (phase == Phase::Warming && phaseCount >= warmLen) { phase = Phase::Fading; phaseCount = 0; }
        else if (phase == Phase::Fading && phaseCount >= fadeLen)
        {
            phase = Phase::Steady;
            active = other;
            phaseCount = 0;
        }
    }
}

} // namespace tt
