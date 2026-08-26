// Across-points batch evaluation of ND polynomial coefficients that the
// caller owns. Eval-only: no fitting machinery. The `FuncEvalND` AoS,
// scatter-AoS and SoA batch overloads forward here; external coefficient
// owners pass an mdspan plus a `domain_nd_view` to the same entry points.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>

#include <poet/poet.hpp>
#include <xsimd/xsimd.hpp>

#include "helpers.hpp"
#include "mdspan.hpp"
#include "poly_eval.hpp"
#include "tags.hpp"

namespace poly_eval {

/**
 * @brief Domain map for coefficients fitted on `[a, b]`.
 *
 * The map to the canonical domain is `fma(2, x, -sum_endpoints) * inv_span`
 * per axis. `identity == true` skips the map (coefficients fitted on
 * `[-1, 1]`, or the map fused into the coefficients).
 */
template<class Scalar, std::size_t DIM>
struct domain_nd_view {
    std::array<Scalar, DIM> inv_span{};
    std::array<Scalar, DIM> sum_endpoints{};
    bool identity = true;
};

/**
 * @brief mdspan over an external Horner-order ND coefficient buffer.
 *
 * The extents match `FuncEvalND` storage: `NCOEFFS^DIM x OUT_DIM`,
 * `layout_right`. `NCOEFFS == 0` selects runtime extents built from
 * @p nCoeffsRt.
 */
namespace detail {

template<std::size_t NCOEFFS, std::size_t OUT_DIM, std::size_t... Is>
PF_CXX20_CONSTEVAL auto makeStaticCoeffExtents(std::index_sequence<Is...>) {
    return stdex::extents<std::size_t, ((void)Is, NCOEFFS)..., OUT_DIM>{};
}

template<std::size_t DIM, std::size_t... Is>
constexpr auto makeDynamicCoeffExtents(int nCoeffsRt, std::size_t outDim, std::index_sequence<Is...>) noexcept {
    return stdex::dextents<std::size_t, DIM + 1>{((void)Is, static_cast<std::size_t>(nCoeffsRt))..., outDim};
}

} // namespace detail

template<std::size_t NCOEFFS, std::size_t DIM, std::size_t OUT_DIM, class Scalar>
[[nodiscard]] PF_ALWAYS_INLINE constexpr auto make_coeffs_mdspan(Scalar *coeffs, [[maybe_unused]] int nCoeffsRt) noexcept {
    if constexpr (NCOEFFS > 0) {
        using extents_t = decltype(detail::makeStaticCoeffExtents<NCOEFFS, OUT_DIM>(std::make_index_sequence<DIM>{}));
        return stdex::mdspan<Scalar, extents_t, stdex::layout_right>{coeffs, extents_t{}};
    } else {
        const auto ext = detail::makeDynamicCoeffExtents<DIM>(nCoeffsRt, OUT_DIM, std::make_index_sequence<DIM>{});
        return stdex::mdspan<Scalar, decltype(ext), stdex::layout_right>{coeffs, ext};
    }
}

namespace detail {

// Scalar single-point evaluation over a coefficient mdspan and domain view.
template<std::size_t NCOEFFS, std::size_t OUT_DIM, ScalarKernel SK, std::size_t HK, class Mdspan, class InScalar,
         std::size_t DIM>
PF_ALWAYS_INLINE constexpr auto horner_nd_view_scalar(const Mdspan &coeffs, int nCoeffsRt,
                                                      const domain_nd_view<InScalar, DIM> &dom,
                                                      const std::array<InScalar, DIM> &x) noexcept {
    using Scalar = std::remove_const_t<typename Mdspan::value_type>;
    using CanonicalOutput = std::array<Scalar, OUT_DIM>;
    const auto xm = dom.identity
                        ? x
                        : polyfit::internal::helpers::mapFromDomainArray(x, dom.inv_span, dom.sum_endpoints);
    if constexpr (SK == ScalarKernel::Horner) {
        return poly_eval::horner<NCOEFFS, true, CanonicalOutput>(xm, coeffs, nCoeffsRt);
    } else {
        return poly_eval::hybrid_nd<NCOEFFS, CanonicalOutput, std::array<InScalar, DIM>, Mdspan, HK>(xm, coeffs,
                                                                                                     nCoeffsRt);
    }
}

// Across-points SIMD tile driver shared by the AoS, scatter-AoS, and SoA
// fronts. Each lane of `batch_t` carries one point's intermediate Horner
// accumulator, so FMAs become packed. StoreBatch receives
// (base, std::array<batch_t, OUT_DIM>); StoreScalar receives
// (i, std::array<Scalar, OUT_DIM>). PF_FLATTEN forces inlining so the store
// callback fuses into the tile loop and codegen does not outline.
template<std::size_t NCOEFFS, std::size_t OUT_DIM, ScalarKernel SK, std::size_t HK, class Mdspan, class InScalar,
         std::size_t DIM, class StoreBatch, class StoreScalar>
PF_FLATTEN constexpr void horner_nd_batch_impl(const Mdspan &coeffs, const domain_nd_view<InScalar, DIM> &dom,
                                               const std::array<InScalar, DIM> *pts, std::size_t count,
                                               StoreBatch storeBatch, StoreScalar storeScalar) noexcept {
    using Scalar = std::remove_const_t<typename Mdspan::value_type>;
    using batch_t = xsimd::batch<Scalar>;
    constexpr std::size_t B = batch_t::size;
    const int nCoeffsRt = (NCOEFFS ? static_cast<int>(NCOEFFS) : static_cast<int>(coeffs.extent(0)));
    if constexpr (B <= 1) {
        for (std::size_t i = 0; i < count; ++i)
            storeScalar(i, horner_nd_view_scalar<NCOEFFS, OUT_DIM, SK, HK>(coeffs, nCoeffsRt, dom, pts[i]));
        return;
    } else {
        constexpr std::size_t kAlign = batch_t::arch_type::alignment();
        // Independent Horner accumulator chains per outer iteration give the
        // FMA-throughput-bound inner loop more ILP. Higher DIM/OUT_DIM grows
        // accumulator register footprint per chain, so unrolling spills —
        // hand-tuned: hank103 DIM=1 OUT_DIM=4 prefers U=1 over U=2.
        constexpr std::size_t U = (DIM <= 2 && OUT_DIM <= 2) ? 2 : 1;

        auto loadPointsAt = [&](std::size_t base) PF_ALWAYS_INLINE_LAMBDA {
            alignas(kAlign) Scalar buf[DIM][B];
            for (std::size_t b = 0; b < B; ++b)
                for (std::size_t d = 0; d < DIM; ++d)
                    buf[d][b] = pts[base + b][d];

            std::array<batch_t, DIM> x_v;
            for (std::size_t d = 0; d < DIM; ++d) {
                batch_t v = batch_t::load_aligned(buf[d]);
                if (!dom.identity) {
                    v = polyfit::internal::helpers::mapFromDomainScalar<batch_t, Scalar>(
                        v, static_cast<Scalar>(dom.inv_span[d]), static_cast<Scalar>(dom.sum_endpoints[d]));
                }
                x_v[d] = v;
            }
            return x_v;
        };

        std::size_t i = 0;
        for (; i + U * B <= count; i += U * B) {
            std::array<std::array<batch_t, DIM>, U> x_vU;
            for (std::size_t u = 0; u < U; ++u) x_vU[u] = loadPointsAt(i + u * B);

            std::array<std::array<batch_t, OUT_DIM>, U> resU;
            for (std::size_t u = 0; u < U; ++u) {
                resU[u] = detail::horner_nd_acrossPts_multi<DIM, NCOEFFS, OUT_DIM, batch_t>(x_vU[u], coeffs,
                                                                                            nCoeffsRt);
            }
            for (std::size_t u = 0; u < U; ++u) storeBatch(i + u * B, resU[u]);
        }
        for (; i + B <= count; i += B) {
            auto x_v = loadPointsAt(i);
            auto res = detail::horner_nd_acrossPts_multi<DIM, NCOEFFS, OUT_DIM, batch_t>(x_v, coeffs, nCoeffsRt);
            storeBatch(i, res);
        }
        // Tail: scalar fallback for the leftover < B points.
        for (; i < count; ++i)
            storeScalar(i, horner_nd_view_scalar<NCOEFFS, OUT_DIM, SK, HK>(coeffs, nCoeffsRt, dom, pts[i]));
    }
}

} // namespace detail

/**
 * @brief Contiguous AoS batch evaluation over external ND coefficients.
 *
 * @p coeffs is an mdspan in `FuncEvalND` layout (see `make_coeffs_mdspan`);
 * result for @p pts[k] lands in @p out[k].
 *
 * @tparam Tag linkage anchor (here and on the two overloads below). An
 * internal-linkage Tag makes the specialization and everything it
 * instantiates internal, so a per-ISA TU exports no weak symbols that could
 * collide across TUs compiled at different `-march` levels. `void` keeps the
 * ordinary external-linkage behavior.
 */
template<std::size_t NCOEFFS, std::size_t OUT_DIM, ScalarKernel SK = ScalarKernel::Horner, std::size_t HK = 0,
         class Tag = void, class Mdspan, class InScalar, std::size_t DIM>
constexpr void horner_nd_batch(const Mdspan &coeffs, const domain_nd_view<InScalar, DIM> &dom,
                               const std::array<InScalar, DIM> *pts,
                               std::array<std::remove_const_t<typename Mdspan::value_type>, OUT_DIM> *out,
                               std::size_t count) noexcept {
    using Scalar = std::remove_const_t<typename Mdspan::value_type>;
    using batch_t = xsimd::batch<Scalar>;
    constexpr std::size_t B = batch_t::size;
    detail::horner_nd_batch_impl<NCOEFFS, OUT_DIM, SK, HK>(
        coeffs, dom, pts, count,
        [&](std::size_t base, const std::array<batch_t, OUT_DIM> &res) PF_ALWAYS_INLINE_LAMBDA {
            // SoA -> AoS shuffle on the output side: spill each per-component
            // batch to aligned scratch, then re-interleave into the caller's
            // AoS buffer. Scalar stores beat xsimd transpose helpers at low
            // OUT_DIM (transpose adds cross-lane shuffles for a tiny re-pack).
            alignas(batch_t::arch_type::alignment()) Scalar outbuf[OUT_DIM][B];
            poet::static_for<OUT_DIM>([&](auto d) { res[d].store_aligned(outbuf[d]); });
            for (std::size_t b = 0; b < B; ++b)
                poet::static_for<OUT_DIM>([&](auto d) { out[base + b][d] = outbuf[d][b]; });
        },
        [&](std::size_t i, const std::array<Scalar, OUT_DIM> &r) PF_ALWAYS_INLINE_LAMBDA { out[i] = r; });
}

/**
 * @brief Scatter-write AoS batch evaluation — result for @p pts[k] lands in @p out[perm[k]].
 */
template<std::size_t NCOEFFS, std::size_t OUT_DIM, ScalarKernel SK = ScalarKernel::Horner, std::size_t HK = 0,
         class Tag = void, class Mdspan, class InScalar, std::size_t DIM>
constexpr void horner_nd_batch(const Mdspan &coeffs, const domain_nd_view<InScalar, DIM> &dom,
                               const std::array<InScalar, DIM> *pts,
                               std::array<std::remove_const_t<typename Mdspan::value_type>, OUT_DIM> *out,
                               const std::uint32_t *perm, std::size_t count) noexcept {
    using Scalar = std::remove_const_t<typename Mdspan::value_type>;
    using batch_t = xsimd::batch<Scalar>;
    constexpr std::size_t B = batch_t::size;
    detail::horner_nd_batch_impl<NCOEFFS, OUT_DIM, SK, HK>(
        coeffs, dom, pts, count,
        [&](std::size_t base, const std::array<batch_t, OUT_DIM> &res) PF_ALWAYS_INLINE_LAMBDA {
            alignas(batch_t::arch_type::alignment()) Scalar outbuf[OUT_DIM][B];
            poet::static_for<OUT_DIM>([&](auto d) { res[d].store_aligned(outbuf[d]); });
            for (std::size_t b = 0; b < B; ++b) {
                const std::uint32_t dst = perm[base + b];
                poet::static_for<OUT_DIM>([&](auto d) { out[dst][d] = outbuf[d][b]; });
            }
        },
        [&](std::size_t i, const std::array<Scalar, OUT_DIM> &r) PF_ALWAYS_INLINE_LAMBDA { out[perm[i]] = r; });
}

/**
 * @brief SoA-output batch evaluation — component @p d of point @p k lands in @p soa_out[d][k].
 */
template<std::size_t NCOEFFS, std::size_t OUT_DIM, ScalarKernel SK = ScalarKernel::Horner, std::size_t HK = 0,
         class Tag = void, class Mdspan, class InScalar, std::size_t DIM>
constexpr void horner_nd_batch_soa(const Mdspan &coeffs, const domain_nd_view<InScalar, DIM> &dom,
                                   const std::array<InScalar, DIM> *pts,
                                   std::array<std::remove_const_t<typename Mdspan::value_type> *, OUT_DIM> soa_out,
                                   std::size_t count) noexcept {
    using Scalar = std::remove_const_t<typename Mdspan::value_type>;
    using batch_t = xsimd::batch<Scalar>;
    detail::horner_nd_batch_impl<NCOEFFS, OUT_DIM, SK, HK>(
        coeffs, dom, pts, count,
        [&](std::size_t base, const std::array<batch_t, OUT_DIM> &res) PF_ALWAYS_INLINE_LAMBDA {
            // Per-component stride-1 store — direct, no AoS shuffle.
            poet::static_for<OUT_DIM>([&](auto d) { res[d].store_unaligned(soa_out[d] + base); });
        },
        [&](std::size_t i, const std::array<Scalar, OUT_DIM> &r) PF_ALWAYS_INLINE_LAMBDA {
            for (std::size_t d = 0; d < OUT_DIM; ++d) soa_out[d][i] = r[d];
        });
}

} // namespace poly_eval
