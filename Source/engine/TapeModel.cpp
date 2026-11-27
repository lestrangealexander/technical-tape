// Technical Tape engine — model files (.ttm) and their preparation for a host sample rate
#include "TapeModel.h"

#include <atomic>
#include "Design.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>

namespace tt
{
namespace
{
    struct Reader
    {
        const uint8_t* p;
        size_t n, pos = 0;
        bool ok = true;

        template <typename T> T get()
        {
            T v {};
            if (pos + sizeof (T) > n) { ok = false; return v; }
            std::memcpy (&v, p + pos, sizeof (T));
            pos += sizeof (T);
            return v;
        }
        bool floats (std::vector<float>& out, size_t count)
        {
            if (pos + 4 * count > n) { ok = false; return false; }
            out.resize (count);
            std::memcpy (out.data(), p + pos, 4 * count);
            pos += 4 * count;
            return true;
        }
        bool halves (std::vector<float>& out, size_t count)
        {
            if (pos + 2 * count > n) { ok = false; return false; }
            out.resize (count);
            for (size_t i = 0; i < count; ++i)
            {
                uint16_t h;
                std::memcpy (&h, p + pos + 2 * i, 2);
                out[i] = halfToFloat (h);
            }
            pos += 2 * count;
            return true;
        }
        static float halfToFloat (uint16_t h)
        {
            const uint32_t sign = (uint32_t) (h & 0x8000u) << 16;
            const uint32_t exp = (h >> 10) & 0x1Fu;
            const uint32_t man = h & 0x3FFu;
            uint32_t bits;
            if (exp == 0)
            {
                if (man == 0) bits = sign;
                else
                {
                    const float v = std::ldexp ((float) man, -24);
                    return sign ? -v : v;
                }
            }
            else if (exp == 31) bits = sign | 0x7F800000u | (man << 13);
            else bits = sign | ((exp + 112u) << 23) | (man << 13);
            float f;
            std::memcpy (&f, &bits, 4);
            return f;
        }
    };
} // namespace

float ModelData::autoGainDb (float driveDb) const noexcept
{
    if (agGainDb.empty()) return -driveDb;
    const float idx = (driveDb - agDrive0) / agStep;
    if (idx <= 0.0f) return agGainDb.front() + (agDrive0 - driveDb);          // keep compensating the level
    const float last = (float) (agGainDb.size() - 1);
    if (idx >= last) return agGainDb.back() - (driveDb - (agDrive0 + last * agStep));
    const int i = (int) idx;
    const float t = idx - (float) i;
    return agGainDb[(size_t) i] + t * (agGainDb[(size_t) i + 1] - agGainDb[(size_t) i]);
}

std::string parseModel (const void* data, size_t size, ModelData& out)
{
    const auto* b = static_cast<const uint8_t*> (data);
    if (size < 8 || std::memcmp (b, "TTM1", 4) != 0) return "not a Technical Tape model";
    Reader top { b, size, 4 };
    const auto version = top.get<uint32_t>();
    if (version != 2) return "unsupported model version " + std::to_string (version);
    static std::atomic<uint64_t> serials { 0 };
    bool haveHead = false;
    out = ModelData();
    out.serial = ++serials;
    while (top.pos + 8 <= size)
    {
        char tag[5] = {};
        std::memcpy (tag, b + top.pos, 4);
        uint32_t len;
        std::memcpy (&len, b + top.pos + 4, 4);
        const size_t start = top.pos + 8;
        if (start + len > size) return std::string ("truncated chunk ") + tag;
        Reader r { b + start, len };
        const std::string t (tag);
        if (t == "HEAD")
        {
            out.sampleRate = r.get<double>();
            out.refPeak = (float) r.get<double>();
            out.numLevels = r.get<int32_t>();
            out.level0Db = r.get<float>();
            out.levelStepDb = r.get<float>();
            (void) r.get<int32_t>();                    // number of branches (counted from BRCH chunks)
            out.preDelay = r.get<int32_t>();
            out.lookAheadMs = r.get<float>();
            out.releaseMs = r.get<float>();
            out.capturedLowVu = r.get<float>();
            out.capturedHighVu = r.get<float>();
            out.flags = r.get<int32_t>();
            out.margin = r.get<float>();
            haveHead = r.ok;
        }
        else if (t == "PEMP")
        {
            const int n = r.get<int32_t>();
            if (n <= 0 || n > 4096) return "bad pre-emphasis table";
            r.floats (out.pempFreqs, (size_t) n);
            r.floats (out.pempDb, (size_t) n);
        }
        else if (t == "BRCH")
        {
            if (! haveHead) return "branch before header";
            ModelData::Branch br;
            br.order = r.get<int32_t>();
            br.mode = r.get<int32_t>();
            br.count = r.get<int32_t>();
            br.length = r.get<int32_t>();
            if (br.order < 1 || br.order > 31 || br.count <= 0 || br.length <= 0 || br.length > (1 << 20)
                || (br.mode == 1 && br.count != out.numLevels) || (br.mode != 0 && br.mode != 1))
                return "bad branch";
            r.halves (br.kernels, (size_t) br.count * (size_t) br.length);
            if (br.mode == 0)
                r.floats (br.weights, (size_t) out.numLevels * (size_t) br.count);
            out.branches.push_back (std::move (br));
        }
        else if (t == "FLUT")
        {
            const int np = r.get<int32_t>();
            for (int i = 0; i < np && r.ok; ++i)
            {
                const float f = r.get<float>(), a = r.get<float>();
                out.flutterPeaks.emplace_back (f, a);
            }
            const int nc = r.get<int32_t>();
            for (int i = 0; i < nc && r.ok; ++i)
            {
                const float f = r.get<float>(), p = r.get<float>();
                out.flutterCont.emplace_back (f, p);
            }
            out.flutterMaxDelay = r.get<float>();
        }
        else if (t == "AGTB")
        {
            const int n = r.get<int32_t>();
            out.agDrive0 = r.get<float>();
            out.agStep = r.get<float>();
            if (n > 0 && n < 10000) r.floats (out.agGainDb, (size_t) n);
        }
        else if (t == "META")
        {
            out.metaJson.assign (reinterpret_cast<const char*> (b + start), len);
        }
        if (! r.ok) return std::string ("bad chunk ") + tag;
        top.pos = start + len;
    }
    if (! haveHead) return "no header";
    if (out.branches.empty() || out.branches.front().order != 1) return "no linear branch";
    if (out.pempFreqs.empty()) return "no pre-emphasis";
    std::sort (out.branches.begin(), out.branches.end(), [] (const auto& a, const auto& c) { return a.order < c.order; });
    return {};
}

int blockSizeFor (double fs) noexcept
{
    if (fs <= 50000.0) return 128;
    if (fs <= 100000.0) return 256;
    return 512;
}

std::shared_ptr<const PreparedModel> PreparedModel::build (const ModelData& m, double fs, int block)
{
    auto pm = std::make_shared<PreparedModel>();
    pm->source = &m;
    pm->fs = fs;
    pm->block = block;
    pm->fftSize = 2 * block;
    const double ratio = fs / m.sampleRate;
    pm->pemp = design::minimumPhase (m.pempFreqs, m.pempDb, fs, (int) std::lround (256.0 * fs / 48000.0));
    pm->lookAhead = std::clamp ((int) std::lround (m.lookAheadMs * 1e-3 * fs), 1, 80);
    pm->preDelay = (int) std::lround (m.preDelay * ratio);
    pm->releaseCoef = (float) std::exp (-1.0 / (m.releaseMs * 1e-3 * fs));
    pm->margin = m.margin;
    pm->refPeak = m.refPeak;
    pm->level0Db = m.level0Db;
    pm->levelStepDb = m.levelStepDb;
    pm->numLevels = m.numLevels;

    RealFft fft (pm->fftSize);
    AlignedBuffer tmp ((size_t) pm->fftSize);
    for (const auto& b : m.branches)
    {
        Branch pb;
        pb.order = b.order;
        pb.mode = b.mode;
        pb.count = b.count;
        pb.weights = b.weights;
        pb.fftSizeOf = pm->fftSize;
        std::vector<std::vector<float>> ks ((size_t) b.count);
        int len = 0;
        for (int k = 0; k < b.count; ++k)
        {
            ks[(size_t) k] = design::resampleKernel (b.kernels.data() + (size_t) k * (size_t) b.length, b.length, m.sampleRate, fs);
            len = std::max (len, (int) ks[(size_t) k].size());
        }
        pb.partitions = (len + block - 1) / block;
        pb.spectra.allocate ((size_t) b.count * (size_t) pb.partitions * (size_t) pm->fftSize);
        for (int k = 0; k < b.count; ++k)
        {
            const auto& kv = ks[(size_t) k];
            for (int p = 0; p < pb.partitions; ++p)
            {
                tmp.clear();
                for (int i = 0; i < block; ++i)
                {
                    const size_t src = (size_t) p * (size_t) block + (size_t) i;
                    if (src < kv.size()) tmp[(size_t) i] = kv[src];
                }
                fft.forward (tmp.data(), const_cast<float*> (pb.kernelSpectrum (k, p)));
            }
        }
        pm->maxOrder = std::max (pm->maxOrder, b.order);
        pm->maxPartitions = std::max (pm->maxPartitions, pb.partitions);
        pm->maxCount = std::max (pm->maxCount, pb.count);
        pm->branches.push_back (std::move (pb));
    }
    return pm;
}

} // namespace tt
