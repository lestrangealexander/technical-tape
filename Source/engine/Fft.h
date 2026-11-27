// Technical Tape engine — FFT and aligned storage (PFFFT underneath)
#pragma once

#include <cstddef>
#include <utility>

struct PFFFT_Setup;

namespace tt
{
/** SIMD-aligned float array (fixed size once allocated). */
class AlignedBuffer
{
public:
    AlignedBuffer() = default;
    explicit AlignedBuffer (size_t n) { allocate (n); }
    ~AlignedBuffer() { release(); }
    AlignedBuffer (const AlignedBuffer&) = delete;
    AlignedBuffer& operator= (const AlignedBuffer&) = delete;
    AlignedBuffer (AlignedBuffer&& o) noexcept : ptr (std::exchange (o.ptr, nullptr)), count (std::exchange (o.count, 0)) {}
    AlignedBuffer& operator= (AlignedBuffer&& o) noexcept
    {
        if (this != &o) { release(); ptr = std::exchange (o.ptr, nullptr); count = std::exchange (o.count, 0); }
        return *this;
    }

    void allocate (size_t n);       // zeroed
    void clear();
    float* data() noexcept { return ptr; }
    const float* data() const noexcept { return ptr; }
    size_t size() const noexcept { return count; }
    float& operator[] (size_t i) noexcept { return ptr[i]; }
    float operator[] (size_t i) const noexcept { return ptr[i]; }

private:
    void release();
    float* ptr = nullptr;
    size_t count = 0;
};

/** Real FFT of size n (a multiple of 32) in PFFFT's fast internal ("unordered") layout, which is all
    a convolver needs: forward, pointwise multiply-accumulate, inverse. The inverse is unscaled. */
class RealFft
{
public:
    explicit RealFft (int n);
    ~RealFft();
    RealFft (const RealFft&) = delete;
    RealFft& operator= (const RealFft&) = delete;

    int size() const noexcept { return n; }
    void forward (const float* in, float* out) const;
    void inverse (const float* in, float* out) const;
    /** acc += a * b * scale (spectra from forward()). */
    void multiplyAccumulate (const float* a, const float* b, float* acc, float scale) const;

private:
    int n;
    PFFFT_Setup* setup = nullptr;
    mutable AlignedBuffer work;
};

} // namespace tt
