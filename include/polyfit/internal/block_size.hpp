#pragma once

#include <cstddef>

#include "feature_macros.hpp"
#include "tags.hpp"

namespace poly_eval::detail {

// optimal_block_size<>: a consteval picker for the Hybrid block size K, from a
// two-regime cost model. For `NCOEFFS` coefficients split into
// B = ceil(NCOEFFS/K) blocks, the single-stream critical path is roughly
//
//     cycles_chain ~= L_fma * (K - 1 + ceil(log2(B))) + K             (Latency)
//
// with L_fma the target FMA latency. A caller that saturates `SIMD_W`
// independent lanes per FMA pipe runs one long chain per lane instead:
//
//     cycles_simd  ~= (NCOEFFS - 1) / SIMD_W                        (Throughput)
//
// The Balanced policy clamps K so that B accumulators plus ceil(log2(K)) Estrin
// powers fit the NREG vector-register budget. The models are qualitative; the
// K-sweep benchmark validates the choice per target.
//
// The production Hybrid path uses `hybrid_block_size` in poly_eval.hpp instead.

PF_CXX20_CONSTEVAL std::size_t ceil_log2(std::size_t n) noexcept {
    std::size_t l = 0;
    std::size_t v = 1;
    while (v < n) { v <<= 1; ++l; }
    return l;
}

PF_CXX20_CONSTEVAL std::size_t ceil_div(std::size_t a, std::size_t b) noexcept {
    return (a + b - 1) / b;
}

template<std::size_t NCOEFFS, EvalPolicy P>
PF_CXX20_CONSTEVAL std::size_t optimal_block_size_latency() noexcept {
    // Return the integer K in [2, NCOEFFS-1] that minimises K + ceil(log2(B))
    // with B = ceil(NCOEFFS / K). Ties resolve to the smaller K, which costs
    // less register pressure.
    std::size_t bestK = 2;
    std::size_t bestCost = NCOEFFS; // sentinel upper bound
    for (std::size_t k = 2; k < NCOEFFS; ++k) {
        const std::size_t b = ceil_div(NCOEFFS, k);
        const std::size_t cost = k + ceil_log2(b);
        if (cost < bestCost) {
            bestCost = cost;
            bestK = k;
        }
    }
    (void)P;
    return bestK;
}

/// Consteval block-size picker. Inputs:
///   - NCOEFFS : compile-time polynomial coefficient count.
///   - SIMD_W  : SIMD lane width on the target (1 for scalar).
///   - NREG    : available vector register budget (32 on AVX-512/NEON, 16 on
///               AVX2/SSE).
///   - P       : `EvalPolicy::{Latency,Throughput,Balanced}`.
///
/// Returns K in [2, NCOEFFS-1]. For `NCOEFFS <= 4` returns `NCOEFFS`, because a
/// single Horner chain beats an Estrin split there.
template<std::size_t NCOEFFS, std::size_t SIMD_W, std::size_t NREG, EvalPolicy P>
PF_CXX20_CONSTEVAL std::size_t optimal_block_size() noexcept {
    if constexpr (NCOEFFS <= 4) {
        return NCOEFFS;
    } else {
        constexpr std::size_t kMax = NCOEFFS - 1;
        if constexpr (P == EvalPolicy::Throughput) {
            // SIMD_W >= 2: each lane is an independent evaluation, so one long
            // chain per lane keeps max ILP across lanes (degenerate Horner).
            // SIMD_W == 1: no lanes to exploit, so use the Latency-optimal K.
            if constexpr (SIMD_W >= 2) {
                return kMax;
            } else {
                return optimal_block_size_latency<NCOEFFS, P>();
            }
        } else if constexpr (P == EvalPolicy::Latency) {
            const std::size_t k = optimal_block_size_latency<NCOEFFS, P>();
            return k < 2 ? std::size_t{2} : (k > kMax ? kMax : k);
        } else {
            // Balanced: start from the Latency choice, then clamp by register
            // pressure. Live registers at the combine point are roughly
            //   B + ceil(log2(K)) + scratch (+1 for the running power).
            // The clamp requires B + ceil(log2(K)) + 1 <= NREG - 2.
            std::size_t k = optimal_block_size_latency<NCOEFFS, P>();
            for (;;) {
                const std::size_t b = ceil_div(NCOEFFS, k);
                const std::size_t live = b + ceil_log2(k) + 1;
                if (live + 2 <= NREG) break;
                if (k >= kMax) break;
                ++k;
            }
            if (k < 2) k = 2;
            if (k > kMax) k = kMax;
            return k;
        }
    }
}

// Self-tests: sanity checks on the heuristic for representative microarchs.
// AVX-512: SIMD_W = 8, NREG = 32.
static_assert(optimal_block_size<3, 8, 32, EvalPolicy::Latency>() == 3,
              "NCOEFFS<=4 returns NCOEFFS");
static_assert(optimal_block_size<4, 8, 32, EvalPolicy::Balanced>() == 4,
              "NCOEFFS<=4 returns NCOEFFS");

static_assert(optimal_block_size<8, 8, 32, EvalPolicy::Throughput>() == 7,
              "Throughput with wide SIMD picks degenerate K = NCOEFFS-1");
static_assert(optimal_block_size<16, 8, 32, EvalPolicy::Throughput>() == 15,
              "Throughput with wide SIMD picks degenerate K = NCOEFFS-1");

static_assert(optimal_block_size<16, 8, 32, EvalPolicy::Latency>() >= 2
                  && optimal_block_size<16, 8, 32, EvalPolicy::Latency>() <= 15,
              "Latency K is in [2, NCOEFFS-1]");
static_assert(optimal_block_size<32, 8, 32, EvalPolicy::Latency>() >= 2
                  && optimal_block_size<32, 8, 32, EvalPolicy::Latency>() <= 31,
              "Latency K is in [2, NCOEFFS-1]");

// Throughput K >= Latency K when SIMD_W >= 2.
static_assert(optimal_block_size<16, 8, 32, EvalPolicy::Throughput>()
                  >= optimal_block_size<16, 8, 32, EvalPolicy::Latency>(),
              "Throughput K >= Latency K when SIMD_W >= 2");

// AVX2: SIMD_W = 4, NREG = 16.
static_assert(optimal_block_size<8, 4, 16, EvalPolicy::Balanced>() >= 2
                  && optimal_block_size<8, 4, 16, EvalPolicy::Balanced>() <= 7,
              "AVX2 Balanced in range");
static_assert(optimal_block_size<16, 4, 16, EvalPolicy::Balanced>() >= 2
                  && optimal_block_size<16, 4, 16, EvalPolicy::Balanced>() <= 15,
              "AVX2 Balanced in range");

// NEON-fp64: SIMD_W = 2, NREG = 32.
static_assert(optimal_block_size<16, 2, 32, EvalPolicy::Throughput>() == 15,
              "NEON Throughput picks degenerate K");

// Scalar: SIMD_W = 1, Throughput falls back to Latency.
static_assert(optimal_block_size<16, 1, 16, EvalPolicy::Throughput>()
                  == optimal_block_size<16, 1, 16, EvalPolicy::Latency>(),
              "Scalar Throughput == Latency");

} // namespace poly_eval::detail
