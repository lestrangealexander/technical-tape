// Technical Tape engine — uniformly partitioned convolution
#include "Convolver.h"

#include <cstring>

namespace tt
{
void BranchConvolver::prepare (int block, int maxPartitions)
{
    B = block;
    N = 2 * block;
    maxP = maxPartitions;
    fdl.allocate ((size_t) maxPartitions * (size_t) N);
    input.allocate ((size_t) N);
    acc.allocate ((size_t) N);
    time.allocate ((size_t) N);
    head = 0;
}

void BranchConvolver::reset() noexcept
{
    fdl.clear();
    input.clear();
    head = 0;
}

void BranchConvolver::pushBlock (const float* in, const RealFft& fft) noexcept
{
    // input = [previous block | this block]
    std::memmove (input.data(), input.data() + B, (size_t) B * sizeof (float));
    std::memcpy (input.data() + B, in, (size_t) B * sizeof (float));
    const int P = branch != nullptr ? branch->partitions : 1;
    head = head + 1 >= P ? 0 : head + 1;
    fft.forward (input.data(), fdl.data() + (size_t) head * (size_t) N);
}

void BranchConvolver::addKernel (int k, const float* weights, float* out, const RealFft& fft) noexcept
{
    if (branch == nullptr) return;
    const int P = branch->partitions;
    acc.clear();
    const float scale = 1.0f / (float) N;
    int slot = head;
    for (int p = 0; p < P; ++p)
    {
        fft.multiplyAccumulate (fdl.data() + (size_t) slot * (size_t) N, branch->kernelSpectrum (k, p), acc.data(), scale);
        slot = slot == 0 ? P - 1 : slot - 1;
    }
    fft.inverse (acc.data(), time.data());
    const float* t = time.data() + B;
    for (int i = 0; i < B; ++i) out[i] += weights[i] * t[i];
}

} // namespace tt
