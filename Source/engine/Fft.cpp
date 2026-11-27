// Technical Tape engine — FFT and aligned storage
#include "Fft.h"
#include "pffft/pffft.h"

#include <cstring>
#include <new>
#include <stdexcept>

namespace tt
{
void AlignedBuffer::allocate (size_t n)
{
    release();
    if (n == 0) return;
    ptr = static_cast<float*> (pffft_aligned_malloc (n * sizeof (float)));
    if (ptr == nullptr) throw std::bad_alloc();
    count = n;
    clear();
}

void AlignedBuffer::clear()
{
    if (ptr != nullptr) std::memset (ptr, 0, count * sizeof (float));
}

void AlignedBuffer::release()
{
    if (ptr != nullptr) pffft_aligned_free (ptr);
    ptr = nullptr;
    count = 0;
}

RealFft::RealFft (int size) : n (size)
{
    setup = pffft_new_setup (n, PFFFT_REAL);
    if (setup == nullptr) throw std::runtime_error ("PFFFT: unsupported size");
    work.allocate ((size_t) n);
}

RealFft::~RealFft()
{
    if (setup != nullptr) pffft_destroy_setup (setup);
}

void RealFft::forward (const float* in, float* out) const
{
    pffft_transform (setup, in, out, work.data(), PFFFT_FORWARD);
}

void RealFft::inverse (const float* in, float* out) const
{
    pffft_transform (setup, in, out, work.data(), PFFFT_BACKWARD);
}

void RealFft::multiplyAccumulate (const float* a, const float* b, float* acc, float scale) const
{
    pffft_zconvolve_accumulate (setup, a, b, acc, scale);
}

} // namespace tt
