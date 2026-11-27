// Technical Tape engine — uniformly partitioned convolution of one branch signal with a set of kernels
//
// One frequency-domain delay line (the input's recent spectra) is shared by every kernel of the branch,
// so a kernel that becomes active (its weight turns non-zero) is immediately exact: no warm-up.
#pragma once

#include "Fft.h"
#include "TapeModel.h"

namespace tt
{
class BranchConvolver
{
public:
    /** Allocates for the largest branch this convolver will ever run. */
    void prepare (int block, int maxPartitions);
    void setBranch (const PreparedModel::Branch* b) noexcept { branch = b; reset(); }
    void reset() noexcept;

    /** Feeds one block of the branch signal. */
    void pushBlock (const float* in, const RealFft& fft) noexcept;

    /** Adds weights[i] * (kernel k convolved with the input)[i] to out[0..block). */
    void addKernel (int k, const float* weights, float* out, const RealFft& fft) noexcept;

private:
    const PreparedModel::Branch* branch = nullptr;
    int B = 0, N = 0, maxP = 0, head = 0;
    AlignedBuffer fdl, input, acc, time;
};

} // namespace tt
