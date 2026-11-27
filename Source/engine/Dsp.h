// Technical Tape engine — per-sample building blocks
//
// Every block here has a twin in Tools/pipeline (refengine.py / dsp_shared.py); Tests/compare_engine.py
// checks that the C++ engine and the Python reference produce the same output.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace tt
{
//==============================================================================
/** Dot product with eight partial sums, so the compiler vectorises it without -ffast-math. */
inline float dot (const float* a, const float* b, int n) noexcept
{
    float acc[8] = {};
    int i = 0;
    for (; i + 8 <= n; i += 8)
        for (int k = 0; k < 8; ++k)
            acc[k] += a[i + k] * b[i + k];
    float s = ((acc[0] + acc[4]) + (acc[1] + acc[5])) + ((acc[2] + acc[6]) + (acc[3] + acc[7]));
    for (; i < n; ++i) s += a[i] * b[i];
    return s;
}

//==============================================================================
/** Direct-form FIR (the pre-emphasis estimate). */
class Fir
{
public:
    void setTaps (const std::vector<float>& taps)
    {
        h = taps;
        hist.assign (2 * h.size(), 0.0f);
        pos = 0;
    }
    void reset() { std::fill (hist.begin(), hist.end(), 0.0f); pos = 0; }

    float process (float x) noexcept
    {
        const size_t n = h.size();
        if (n == 0) return x;
        // history kept twice so the dot product never wraps
        pos = (pos == 0 ? n : pos) - 1;
        hist[pos] = hist[pos + n] = x;
        return dot (h.data(), hist.data() + pos, (int) n);
    }

private:
    std::vector<float> h, hist;
    size_t pos = 0;
};

//==============================================================================
/** Integer delay line with fixed capacity. */
class Delay
{
public:
    void setLength (int samples) { len = std::max (0, samples); buf.assign ((size_t) len + 1, 0.0f); pos = 0; }
    void reset() { std::fill (buf.begin(), buf.end(), 0.0f); pos = 0; }
    float process (float x) noexcept
    {
        if (len == 0) return x;
        buf[pos] = x;
        pos = pos + 1 == buf.size() ? 0 : pos + 1;
        return buf[pos];
    }
    int length() const noexcept { return len; }

private:
    std::vector<float> buf;
    size_t pos = 0;
    int len = 0;
};

//==============================================================================
/** The look-ahead peak detector of the pre-emphasised signal:
    max(analytic envelope, |x|) -> sliding max over A -> moving average over A -> exponential release.
    Same as refengine.detector(). */
class Detector
{
public:
    void prepare (int lookAhead, float releaseCoef)
    {
        A = std::max (1, lookAhead);
        rel = releaseCoef;
        envRing.assign ((size_t) A, 0.0f);
        maxRing.assign ((size_t) A, 0.0f);
        reset();
    }

    void reset()
    {
        for (auto& c : hs) c = {};
        x1 = 0.0f;
        std::fill (envRing.begin(), envRing.end(), 0.0f);
        std::fill (maxRing.begin(), maxRing.end(), 0.0f);
        pos = 0;
        sum = 0.0;
        prev = 0.0f;
    }

    float process (float x) noexcept
    {
        // Hilbert pair: chain 0 taken one sample later
        float i = x, q = x;
        for (int s = 0; s < 4; ++s) i = hs[s].process (i, kA2[0][s]);
        for (int s = 0; s < 4; ++s) q = hs[4 + s].process (q, kA2[1][s]);
        const float iDel = x1;
        x1 = i;
        const float env = std::max (std::sqrt (iDel * iDel + q * q), std::abs (x));
        envRing[pos] = env;
        float sm = 0.0f;                                   // sliding max over the last A (A <= 80)
        for (float v : envRing) sm = v > sm ? v : sm;
        sum += (double) sm - (double) maxRing[pos];        // moving average over the last A maxima
        maxRing[pos] = sm;
        pos = pos + 1 == envRing.size() ? 0 : pos + 1;
        const float ma = (float) (sum / (double) A);
        prev = ma > prev * rel ? ma : prev * rel;
        return std::max (prev, 1e-6f);
    }

private:
    struct Section
    {
        float x1 = 0, x2 = 0, y1 = 0, y2 = 0;
        float process (float x, float a2) noexcept
        {
            const float y = a2 * (x + y2) - x2;
            x2 = x1; x1 = x;
            y2 = y1; y1 = y;
            return y;
        }
    };
    static constexpr float kA2[2][4] = {
        { 0.6923878f * 0.6923878f, 0.9360654322959f * 0.9360654322959f, 0.9882295226860f * 0.9882295226860f, 0.9987488452737f * 0.9987488452737f },
        { 0.4021921162426f * 0.4021921162426f, 0.8561710882420f * 0.8561710882420f, 0.9722909545651f * 0.9722909545651f, 0.9952884791278f * 0.9952884791278f } };

    Section hs[8];
    float x1 = 0.0f;
    int A = 1;
    float rel = 0.0f;
    std::vector<float> envRing, maxRing;
    size_t pos = 0;
    double sum = 0.0;
    float prev = 0.0f;
};

//==============================================================================
/** x2 interpolation stage: zero-stuff then 2 * FIR (as dsp_shared.upsample). Polyphase. */
class UpStage
{
public:
    void setTaps (const double* taps, int numTaps)
    {
        ph[0].clear(); ph[1].clear();
        for (int i = 0; i < numTaps; ++i)
            ph[i & 1].push_back ((float) (2.0 * taps[i]));
        len = (int) std::max (ph[0].size(), ph[1].size());
        ph[1].resize ((size_t) len, 0.0f);
        ph[0].resize ((size_t) len, 0.0f);
        hist.assign (2 * (size_t) len, 0.0f);
        pos = 0;
    }
    void reset() { std::fill (hist.begin(), hist.end(), 0.0f); pos = 0; }

    /** One input sample in, two output samples out. */
    void process (float x, float* out2) noexcept
    {
        pos = (pos == 0 ? (size_t) len : pos) - 1;
        hist[pos] = hist[pos + (size_t) len] = x;
        const float* s = hist.data() + pos;
        out2[0] = dot (ph[0].data(), s, len);
        out2[1] = dot (ph[1].data(), s, len);
    }

private:
    std::vector<float> ph[2], hist;
    size_t pos = 0;
    int len = 0;
};

/** /2 decimation stage: FIR then keep even samples (as dsp_shared.downsample), in polyphase form:
    y[k] = sum_j h[2j] x[2k-2j] + sum_j h[2j+1] x[2k-2j-1]. Each phase keeps only the span of its
    non-zero taps (a half-band's odd phase is its centre tap alone). */
class DownStage
{
public:
    void setTaps (const double* taps, int numTaps)
    {
        for (int q = 0; q < 2; ++q)
        {
            std::vector<float> c;
            for (int i = q; i < numTaps; i += 2) c.push_back ((float) taps[i]);
            int a = 0, b = (int) c.size();
            // a half-band's zero taps are exact zeros by design
            while (a < b && std::fpclassify (c[(size_t) a]) == FP_ZERO) ++a;
            while (b > a && std::fpclassify (c[(size_t) b - 1]) == FP_ZERO) --b;
            off[q] = a;
            ph[q].assign (c.begin() + a, c.begin() + b);
            len[q] = (int) ph[q].size();
            cap[q] = std::max (1, off[q] + len[q]);
            hist[q].assign (2 * (size_t) cap[q], 0.0f);
            pos[q] = 0;
        }
    }
    void reset() { for (int q = 0; q < 2; ++q) { std::fill (hist[q].begin(), hist[q].end(), 0.0f); pos[q] = 0; } }

    /** Two input samples in (the even one first), one output sample out. */
    float process (float even, float odd) noexcept
    {
        push (0, even);                                  // E[j] = even sample k-j
        const float y = dot (ph[0].data(), hist[0].data() + pos[0] + off[0], len[0])
                      + dot (ph[1].data(), hist[1].data() + pos[1] + off[1], len[1]);   // O[j] = odd sample k-1-j
        push (1, odd);
        return y;
    }

private:
    void push (int q, float x) noexcept
    {
        pos[q] = (pos[q] == 0 ? (size_t) cap[q] : pos[q]) - 1;
        hist[q][pos[q]] = hist[q][pos[q] + (size_t) cap[q]] = x;
    }
    std::vector<float> ph[2], hist[2];
    size_t pos[2] = {};
    int off[2] = {}, len[2] = {}, cap[2] = { 1, 1 };
};

//==============================================================================
/** Ring of recent values of a base-rate signal, read back with linear interpolation (the detector
    level, read at the fractional input time of each oversampled sample). */
class History
{
public:
    void prepare (int capacity) { buf.assign ((size_t) capacity, 0.0f); n = 0; }
    void reset() { std::fill (buf.begin(), buf.end(), 0.0f); n = 0; first = 0.0f; }
    void push (float v) noexcept
    {
        if (n == 0) first = v;
        buf[(size_t) (n % (int64_t) buf.size())] = v;
        ++n;
    }
    /** Value at base time t (samples since the start; fractional), clamped to what exists. */
    float at (double t) const noexcept
    {
        if (t <= 0.0) return n > 0 ? first : 0.0f;
        const double last = (double) (n - 1);
        if (t >= last) return buf[(size_t) ((n - 1) % (int64_t) buf.size())];
        const int64_t i = (int64_t) t;
        const float fr = (float) (t - (double) i);
        const float a = buf[(size_t) (i % (int64_t) buf.size())], b = buf[(size_t) ((i + 1) % (int64_t) buf.size())];
        return a + fr * (b - a);
    }
    int64_t count() const noexcept { return n; }

private:
    std::vector<float> buf;
    int64_t n = 0;
    float first = 0.0f;
};

} // namespace tt
