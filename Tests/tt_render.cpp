// tt_render — runs the plugin's engine offline on a WAV file (tests, validation, benchmarks).
//
//   tt_render <model.ttm>[,<model2.ttm>...] <in.wav> <out.wav> [options]
//     --os N            oversampling factor 1, 2, 4, 8 or 16 (default 8)
//     --drive dB        INPUT (default 0)
//     --output dB       OUTPUT (default 0)
//     --wf X            wow & flutter amount 0..2 (default 1); --no-flutter = 0
//     --no-autogain     auto gain off (default on)
//     --block N         host block size (default 512); --random-blocks varies it per call (1..N)
//     --switch-model I  switch to model I (of the comma-separated list) at --at seconds
//     --switch-os N     switch to OS factor N at --at seconds
//     --at s            time of the switch (default 1.0)
//     --bypass-at s     engage BYPASS at s
//   Prints "latency N" (samples) and the real-time factor.
#include "../Source/engine/TapeEngine.h"

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace
{
struct Wav
{
    double fs = 48000.0;
    int channels = 1;
    std::vector<std::vector<float>> data;
};

bool readWav (const std::string& path, Wav& w)
{
    std::ifstream f (path, std::ios::binary);
    if (! f) return false;
    std::vector<uint8_t> b ((std::istreambuf_iterator<char> (f)), std::istreambuf_iterator<char>());
    if (b.size() < 12 || std::memcmp (b.data(), "RIFF", 4) != 0 || std::memcmp (b.data() + 8, "WAVE", 4) != 0) return false;
    auto u16 = [&] (size_t p) { return (uint32_t) b[p] | ((uint32_t) b[p + 1] << 8); };
    auto u32 = [&] (size_t p) { return u16 (p) | (u16 (p + 2) << 16); };
    int fmt = 0, bits = 0;
    size_t p = 12;
    while (p + 8 <= b.size())
    {
        const uint32_t len = u32 (p + 4);
        if (std::memcmp (b.data() + p, "fmt ", 4) == 0)
        {
            fmt = (int) u16 (p + 8);
            w.channels = (int) u16 (p + 10);
            w.fs = (double) u32 (p + 12);
            bits = (int) u16 (p + 22);
            if (fmt == 0xFFFE) fmt = (int) u16 (p + 32);
        }
        else if (std::memcmp (b.data() + p, "data", 4) == 0)
        {
            const size_t bytes = std::min<size_t> (len, b.size() - p - 8);
            const int bps = bits / 8;
            const size_t frames = bytes / (size_t) (bps * w.channels);
            w.data.assign ((size_t) w.channels, std::vector<float> (frames));
            const uint8_t* d = b.data() + p + 8;
            for (size_t i = 0; i < frames; ++i)
                for (int c = 0; c < w.channels; ++c)
                {
                    const uint8_t* s = d + (i * (size_t) w.channels + (size_t) c) * (size_t) bps;
                    float v = 0.0f;
                    if (fmt == 3 && bits == 32) std::memcpy (&v, s, 4);
                    else if (fmt == 3 && bits == 64) { double dv; std::memcpy (&dv, s, 8); v = (float) dv; }
                    else if (bits == 16) v = (float) (int16_t) (s[0] | (s[1] << 8)) / 32768.0f;
                    else if (bits == 24) v = (float) ((int32_t) ((uint32_t) s[0] << 8 | (uint32_t) s[1] << 16 | (uint32_t) s[2] << 24) >> 8) / 8388608.0f;
                    else if (bits == 32) { int32_t iv; std::memcpy (&iv, s, 4); v = (float) ((double) iv / 2147483648.0); }
                    w.data[(size_t) c][i] = v;
                }
            return true;
        }
        p += 8 + len + (len & 1);
    }
    return false;
}

bool writeWav (const std::string& path, const Wav& w)
{
    std::ofstream f (path, std::ios::binary);
    if (! f) return false;
    const size_t frames = w.data.empty() ? 0 : w.data[0].size();
    const uint32_t bytes = (uint32_t) (frames * (size_t) w.channels * 4);
    auto put32 = [&] (uint32_t v) { f.write ((const char*) &v, 4); };
    auto put16 = [&] (uint16_t v) { f.write ((const char*) &v, 2); };
    f.write ("RIFF", 4); put32 (36 + bytes); f.write ("WAVE", 4);
    f.write ("fmt ", 4); put32 (16); put16 (3); put16 ((uint16_t) w.channels); put32 ((uint32_t) w.fs);
    put32 ((uint32_t) w.fs * (uint32_t) w.channels * 4); put16 ((uint16_t) (w.channels * 4)); put16 (32);
    f.write ("data", 4); put32 (bytes);
    for (size_t i = 0; i < frames; ++i)
        for (int c = 0; c < w.channels; ++c)
            f.write ((const char*) &w.data[(size_t) c][i], 4);
    return true;
}

bool loadModel (const std::string& path, tt::ModelData& m)
{
    std::ifstream f (path, std::ios::binary);
    if (! f) { std::cerr << "cannot open " << path << "\n"; return false; }
    std::vector<char> b ((std::istreambuf_iterator<char> (f)), std::istreambuf_iterator<char>());
    const auto err = tt::parseModel (b.data(), b.size(), m);
    if (! err.empty()) { std::cerr << path << ": " << err << "\n"; return false; }
    return true;
}
} // namespace

int main (int argc, char** argv)
{
    if (argc < 4)
    {
        std::cerr << "usage: tt_render <model.ttm>[,more] <in.wav> <out.wav> [--os N] [--drive dB] [--output dB]"
                     " [--wf X | --no-flutter] [--no-autogain] [--block N] [--random-blocks]"
                     " [--switch-model I] [--switch-os N] [--at s] [--bypass-at s]\n";
        return 2;
    }
    std::vector<std::string> modelPaths;
    {
        std::string s = argv[1];
        size_t a = 0;
        while (a <= s.size())
        {
            const size_t e = s.find (',', a);
            modelPaths.push_back (s.substr (a, e == std::string::npos ? std::string::npos : e - a));
            if (e == std::string::npos) break;
            a = e + 1;
        }
    }
    tt::EngineParams p;
    int block = 512, switchModel = -1, switchOs = -1;
    bool randomBlocks = false;
    double at = 1.0, bypassAt = -1.0;
    for (int i = 4; i < argc; ++i)
    {
        const std::string a = argv[i];
        auto next = [&] { return i + 1 < argc ? std::string (argv[++i]) : std::string ("0"); };
        auto osIndex = [] (int f) { int k = 0; while (k < 4 && tt::TapeEngine::kOsFactors[k] < f) ++k; return k; };
        if (a == "--os") p.osIndex = osIndex (std::atoi (next().c_str()));
        else if (a == "--drive") p.inputDb = (float) std::atof (next().c_str());
        else if (a == "--output") p.outputDb = (float) std::atof (next().c_str());
        else if (a == "--wf") p.wowFlutter = (float) std::atof (next().c_str());
        else if (a == "--no-flutter") p.wowFlutter = 0.0f;
        else if (a == "--no-autogain") p.autoGain = false;
        else if (a == "--block") block = std::max (1, std::atoi (next().c_str()));
        else if (a == "--random-blocks") randomBlocks = true;
        else if (a == "--switch-model") switchModel = std::atoi (next().c_str());
        else if (a == "--switch-os") switchOs = osIndex (std::atoi (next().c_str()));
        else if (a == "--at") at = std::atof (next().c_str());
        else if (a == "--bypass-at") bypassAt = std::atof (next().c_str());
        else { std::cerr << "unknown option " << a << "\n"; return 2; }
    }

    std::vector<tt::ModelData> models (modelPaths.size());
    for (size_t i = 0; i < modelPaths.size(); ++i)
        if (! loadModel (modelPaths[i], models[i])) return 1;
    Wav in;
    if (! readWav (argv[2], in)) { std::cerr << "cannot read " << argv[2] << "\n"; return 1; }

    tt::TapeEngine eng;
    std::vector<const tt::ModelData*> ptrs;
    for (auto& m : models) ptrs.push_back (&m);
    eng.setModels (ptrs);
    eng.setParams (p);
    const auto t0 = std::chrono::steady_clock::now();
    eng.prepare (in.fs, block, in.channels);
    const auto t1 = std::chrono::steady_clock::now();
    if (! eng.isPrepared()) { std::cerr << "engine not prepared\n"; return 1; }
    std::cout << "latency " << eng.latencySamples() << "\n";
    std::cout << "prepare " << std::chrono::duration<double> (t1 - t0).count() << " s\n";

    const size_t frames = in.data[0].size();
    Wav out = in;
    std::vector<float*> ptr ((size_t) in.channels);
    std::mt19937 rng (7);
    size_t pos = 0;
    const auto t2 = std::chrono::steady_clock::now();
    while (pos < frames)
    {
        int n = randomBlocks ? (int) (1 + rng() % (uint32_t) block) : block;
        n = (int) std::min<size_t> ((size_t) n, frames - pos);
        const double t = (double) pos / in.fs;
        tt::EngineParams q = p;
        if (t >= at && switchModel >= 0) q.model = switchModel;
        if (t >= at && switchOs >= 0) q.osIndex = switchOs;
        if (bypassAt >= 0.0 && t >= bypassAt) q.bypass = true;
        eng.setParams (q);
        for (int c = 0; c < in.channels; ++c) ptr[(size_t) c] = out.data[(size_t) c].data() + pos;
        eng.process (ptr.data(), in.channels, n);
        pos += (size_t) n;
    }
    const double secs = std::chrono::duration<double> (std::chrono::steady_clock::now() - t2).count();
    std::cout << "realtime " << ((double) frames / in.fs) / secs << " x (" << in.channels << " ch, "
              << tt::TapeEngine::kOsFactors[p.osIndex] << "x)\n";
    if (! writeWav (argv[3], out)) { std::cerr << "cannot write " << argv[3] << "\n"; return 1; }
    return 0;
}
