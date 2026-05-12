#pragma once

#include <cstddef>

#include "feature_macros.hpp"
#include "tags.hpp"

namespace poly_eval::detail {

// optimal_block_size<>: closed-form, consteval picker for the Hybrid block
// size K. Modelled after the Initial Model from the K-sweep plan:
//
//   For NCOEFFS coefficients split into B = ceil(NCOEFFS/K) blocks, the
//   single-stream Estrin-Horner critical path is roughly
//
//       cycles_chain ~= L_fma * (K - 1 + ceil(log2(B))) + K           (Latency)
//
//   where L_fma is the FMA latency (cycles) on the target microarch. When
//   the caller already saturates `SIMD_W` independent lanes per FMA pipe,
//   the throughput-limited cost is
//
//       cycles_simd  ~= (NCOEFFS - 1) / SIMD_W                        (Throughput)
//
//   (one long dependent chain per lane, all lanes in flight). NREG is the
//   available vector-register budget; the Balanced choice clamps K so that
//   K accumulators + ceil(log2(K)) Estrin powers + a few scratch slots fit
//   below NREG-2.
//
//   These are qualitative models; the actual K is validated empirically per
//   target via the K-sweep benchmark.

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
    // Pick the integer K in [2, NCOEFFS-1] minimising K + ceil(log2(B))
    // where B = ceil(NCOEFFS / K). Ties resolve to the smaller K (lower
    // register pressure).
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
///   - NREG    : available vector register budget (e.g. 32 on AVX-512/NEON,
///               16 on AVX2/SSE).
///   - P       : `EvalPolicy::{Latency,Throughput,Balanced}`.
///
/// Returns K in [2, NCOEFFS-1]. Special cases: `NCOEFFS <= 4` returns
/// `NCOEFFS` (single Horner chain, no Estrin split worthwhile).
template<std::size_t NCOEFFS, std::size_t SIMD_W, std::size_t NREG, EvalPolicy P>
PF_CXX20_CONSTEVAL std::size_t optimal_block_size() noexcept {
    if constexpr (NCOEFFS <= 4) {
        return NCOEFFS;
    } else {
        constexpr std::size_t kMax = NCOEFFS - 1;
        if constexpr (P == EvalPolicy::Throughput) {
            // SIMD_W >= 2: each lane is an independent evaluation, so we want
            // one long chain per lane (degenerate Horner, max ILP across lanes).
            // SIMD_W == 1: fall back to the Latency-optimal K — there are no
            // lanes to exploit, and a balanced split is best.
            if constexpr (SIMD_W >= 2) {
                return kMax;
            } else {
                return optimal_block_size_latency<NCOEFFS, P>();
            }
        } else if constexpr (P == EvalPolicy::Latency) {
            const std::size_t k = optimal_block_size_latency<NCOEFFS, P>();
            return k < 2 ? std::size_t{2} : (k > kMax ? kMax : k);
        } else {
            // Balanced: start from Latency choice, then clamp by register
            // pressure. Live regs at the combine point are roughly
            //   B + ceil(log2(K)) + a few scratch (+1 for the running power).
            // We require B + ceil(log2(K)) + 1 <= NREG - 2.
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
// AVX-512 — SIMD_W = 8, NREG = 32
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

// AVX2 — SIMD_W = 4, NREG = 16
static_assert(optimal_block_size<8, 4, 16, EvalPolicy::Balanced>() >= 2
                  && optimal_block_size<8, 4, 16, EvalPolicy::Balanced>() <= 7,
              "AVX2 Balanced in range");
static_assert(optimal_block_size<16, 4, 16, EvalPolicy::Balanced>() >= 2
                  && optimal_block_size<16, 4, 16, EvalPolicy::Balanced>() <= 15,
              "AVX2 Balanced in range");

// NEON-fp64 — SIMD_W = 2, NREG = 32
static_assert(optimal_block_size<16, 2, 32, EvalPolicy::Throughput>() == 15,
              "NEON Throughput picks degenerate K");

// Scalar — SIMD_W = 1: Throughput must fall back to Latency.
static_assert(optimal_block_size<16, 1, 16, EvalPolicy::Throughput>()
                  == optimal_block_size<16, 1, 16, EvalPolicy::Latency>(),
              "Scalar Throughput == Latency");

} // namespace poly_eval::detail
