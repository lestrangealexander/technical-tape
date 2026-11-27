// Technical Tape engine — filter design done at prepare time
#include "Design.h"

#include <algorithm>
#include <cmath>

namespace tt::design
{
namespace
{
    constexpr double kPi = 3.14159265358979323846;

    double besselI0 (double x)
    {
        double sum = 1.0, term = 1.0;
        const double q = x * x * 0.25;
        for (int k = 1; k < 60; ++k)
        {
            term *= q / ((double) k * (double) k);
            sum += term;
            if (term < 1e-17 * sum) break;
        }
        return sum;
    }

    /** np.interp: linear, held flat outside [xp.front(), xp.back()]. */
    double interp (double x, const std::vector<double>& xp, const std::vector<double>& fp)
    {
        if (x <= xp.front()) return fp.front();
        if (x >= xp.back()) return fp.back();
        const auto it = std::upper_bound (xp.begin(), xp.end(), x);
        const size_t i = (size_t) (it - xp.begin());
        const double t = (x - xp[i - 1]) / (xp[i] - xp[i - 1]);
        return fp[i - 1] + t * (fp[i] - fp[i - 1]);
    }
} // namespace

void fft (std::vector<double>& re, std::vector<double>& im, bool inverse)
{
    const size_t n = re.size();
    for (size_t i = 1, j = 0; i < n; ++i)
    {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { std::swap (re[i], re[j]); std::swap (im[i], im[j]); }
    }
    for (size_t len = 2; len <= n; len <<= 1)
    {
        const double ang = 2.0 * kPi / (double) len * (inverse ? 1.0 : -1.0);
        const double wr = std::cos (ang), wi = std::sin (ang);
        for (size_t i = 0; i < n; i += len)
        {
            double cr = 1.0, ci = 0.0;
            for (size_t k = 0; k < len / 2; ++k)
            {
                const size_t a = i + k, b = i + k + len / 2;
                const double xr = re[b] * cr - im[b] * ci, xi = re[b] * ci + im[b] * cr;
                re[b] = re[a] - xr; im[b] = im[a] - xi;
                re[a] += xr; im[a] += xi;
                const double t = cr * wr - ci * wi;
                ci = cr * wi + ci * wr;
                cr = t;
            }
        }
    }
}

std::vector<float> minimumPhase (const std::vector<float>& freqsHz, const std::vector<float>& magDb, double fs, int numTaps)
{
    size_t nfft = 1;
    while (nfft < (size_t) (32 * numTaps)) nfft <<= 1;
    std::vector<double> xp (freqsHz.size()), fp (magDb.size());
    for (size_t i = 0; i < xp.size(); ++i)
    {
        xp[i] = std::log ((double) freqsHz[i]);
        fp[i] = (double) magDb[i];
    }
    std::vector<double> re (nfft, 0.0), im (nfft, 0.0);
    const size_t half = nfft / 2;
    for (size_t k = 0; k <= half; ++k)
    {
        const double f = (double) k * fs / (double) nfft;
        const double m = interp (std::log (std::max (f, 1e-3)), xp, fp);
        re[k] = m * std::log (10.0) / 20.0;
    }
    for (size_t k = 1; k < half; ++k)
        re[nfft - k] = re[k];
    // real cepstrum
    fft (re, im, true);
    for (size_t k = 0; k < nfft; ++k) { re[k] /= (double) nfft; im[k] = 0.0; }
    std::vector<double> fr (nfft, 0.0), fi (nfft, 0.0);
    fr[0] = re[0];
    for (size_t k = 1; k < half; ++k) fr[k] = 2.0 * re[k];
    fr[half] = re[half];
    fft (fr, fi, false);
    for (size_t k = 0; k < nfft; ++k)
    {
        const double mag = std::exp (fr[k]);
        const double ph = fi[k];
        fr[k] = mag * std::cos (ph);
        fi[k] = mag * std::sin (ph);
    }
    fft (fr, fi, true);
    std::vector<float> h ((size_t) numTaps);
    for (int i = 0; i < numTaps; ++i) h[(size_t) i] = (float) (fr[(size_t) i] / (double) nfft);
    const int nt = std::max (8, numTaps / 8);
    for (int k = 0; k < nt; ++k)
    {
        const double a = kPi * (double) k / (double) (nt - 1);
        h[(size_t) (numTaps - nt + k)] *= (float) (0.5 + 0.5 * std::cos (a));
    }
    return h;
}

std::vector<float> resampleKernel (const float* in, int length, double fsIn, double fsOut)
{
    if (std::abs (fsIn - fsOut) < 1e-6)
        return std::vector<float> (in, in + length);
    const double ratio = fsIn / fsOut;                          // input samples per output sample
    const double cutoff = 0.5 * std::min (1.0, fsOut / fsIn);   // cycles per input sample
    const int half = 32;
    const double beta = 9.0;
    const double i0b = besselI0 (beta);
    constexpr int tableSize = 8192;
    std::vector<double> win (tableSize + 2);
    for (int i = 0; i <= tableSize + 1; ++i)
    {
        const double r = std::min (1.0, (double) i / (double) tableSize);
        win[(size_t) i] = besselI0 (beta * std::sqrt (std::max (0.0, 1.0 - r * r))) / i0b;
    }
    const int outLen = (int) std::ceil ((double) length / ratio);
    std::vector<float> out ((size_t) outLen, 0.0f);
    for (int n = 0; n < outLen; ++n)
    {
        const double t = (double) n * ratio;
        const int c = (int) std::floor (t);
        double acc = 0.0;
        for (int m = c - half + 1; m <= c + half; ++m)
        {
            if (m < 0 || m >= length) continue;
            const double d = t - (double) m;
            const double x = 2.0 * cutoff * d;
            const double s = std::abs (x) < 1e-12 ? 1.0 : std::sin (kPi * x) / (kPi * x);
            const double r = std::abs (d) / (double) half;
            if (r >= 1.0) continue;
            const double pos = r * (double) tableSize;
            const int ti = (int) pos;
            const double w = win[(size_t) ti] + (pos - (double) ti) * (win[(size_t) ti + 1] - win[(size_t) ti]);
            acc += (double) in[m] * 2.0 * cutoff * s * w;
        }
        out[(size_t) n] = (float) (acc * ratio);
    }
    return out;
}

} // namespace tt::design
