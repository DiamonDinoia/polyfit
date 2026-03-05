#pragma once

#include "macros.h"
#include "simd_utils.h"
#include "utils.h"

#include <poet/poet.hpp>

#include <cassert>
#include <cstddef>

//------------------------------------------------------------------------------
// Forward Declarations (Public API)
//------------------------------------------------------------------------------

namespace poly_eval {

// Scalar Horner (one-point)
/**
 * @brief Evaluate a polynomial at a single point using Horner's method (scalar version).
 *
 * @tparam N_total Compile-time number of coefficients (0 for runtime).
 * @tparam OutputType Output value type.
 * @tparam InputType Input value type.
 * @param x The input value at which to evaluate the polynomial.
 * @param c_ptr Pointer to the coefficients array (reversed order: highest degree first).
 * @param c_size Number of coefficients (used if N_total == 0).
 * @return The evaluated polynomial value at x.
 */
template<std::size_t N_total = 0, typename OutputType, typename InputType>
PF_ALWAYS_INLINE constexpr OutputType horner(InputType xin, const OutputType *c_ptr, std::size_t c_size = 0) noexcept;

/**
 * @brief Evaluate a polynomial at multiple points using SIMD-accelerated Horner's method.
 *
 * Coefficients are in reversed order (highest degree first).
 *
 * @tparam N_monomials Compile-time number of monomials (0 for runtime).
 * @tparam pts_aligned Whether input points are SIMD-aligned.
 * @tparam out_aligned Whether output is SIMD-aligned.
 * @tparam UNROLL Unroll factor for SIMD loop (0 for default).
 * @tparam InputType Input value type.
 * @tparam OutputType Output value type.
 * @tparam MapFunc Mapping function applied to each input point used to center the coefficients.
 * @param pts Pointer to input points array.
 * @param out Pointer to output array.
 * @param num_points Number of points to evaluate.
 * @param monomials Pointer to coefficients array (reversed order).
 * @param monomials_size Number of coefficients.
 * @param map_func Function to map/transform each input point (default: identity).
 */
template<std::size_t N_monomials = 0, bool pts_aligned = false, bool out_aligned = false, int UNROLL = 0,
         typename InputType, typename OutputType, typename MapFunc>
PF_ALWAYS_INLINE constexpr void horner(
    const InputType *pts, OutputType *out, std::size_t num_points, const OutputType *monomials,
    std::size_t monomials_size, MapFunc map_func = [](auto v) { return v; }) noexcept;

/**
 * @brief Evaluate multiple polynomials at a single point using Horner's method.
 *
 * Optionally supports scaling of the input.
 *
 * @tparam M_total Compile-time number of polynomials (0 for runtime).
 * @tparam N_total Compile-time number of coefficients per polynomial (0 for runtime).
 * @tparam scaling Whether to apply scaling to the input.
 * @tparam OutputType Output value type.
 * @tparam InputType Input value type.
 * @param x The input value at which to evaluate.
 * @param coeffs Pointer to the coefficients array (row-major: M polynomials × N coefficients).
 * @param out Pointer to output array (size M).
 * @param M Number of polynomials (used if M_total == 0).
 * @param N Number of coefficients per polynomial (used if N_total == 0).
 * @param low Optional pointer to scaling lower bounds (used if scaling).
 * @param hi Optional pointer to scaling upper bounds (used if scaling).
 */
template<std::size_t M_total = 0, std::size_t N_total = 0, bool scaling = false, typename OutputType,
         typename InputType>
PF_ALWAYS_INLINE constexpr void horner_many(InputType xin, const OutputType *coeffs, OutputType *out, std::size_t M = 0,
                                            std::size_t N = 0, const InputType *low = nullptr,
                                            const InputType *high = nullptr) noexcept;

/**
 * @brief Evaluate multiple polynomials at multiple points (transposed layout).
 *
 * Supports both scalar and SIMD paths.
 *
 * @tparam M_total Compile-time number of points (0 for runtime).
 * @tparam N_total Compile-time number of coefficients (0 for runtime).
 * @tparam simd_width SIMD width (0 for scalar).
 * @tparam Out Output value type.
 * @tparam In Input value type.
 * @param x Pointer to input points array (size M).
 * @param c Pointer to coefficients array (row-major: N rows × M columns).
 * @param out Pointer to output array (size M).
 * @param M Number of points (used if M_total == 0).
 * @param N Number of coefficients (used if N_total == 0).
 */
template<std::size_t M_total = 0, std::size_t N_total = 0, std::size_t simd_width = 0, bool aligned = false,
         typename Out, typename In>
PF_ALWAYS_INLINE constexpr void horner_transposed(const In *xin, const Out *coeffs, Out *out, std::size_t M = 0,
                                                  std::size_t N = 0) noexcept;

/**
 * @brief Evaluate a multivariate polynomial using N-dimensional Horner's method.
 *
 * @tparam DegCT Compile-time degree (0 for runtime).
 * @tparam OutT Output type (e.g., scalar or std::array).
 * @tparam InVec Input vector type (e.g., std::array).
 * @tparam Mdspan Coefficient tensor type (e.g., mdspan).
 * @param x Input vector of variables.
 * @param coeffs Coefficient tensor (last dimension is output dimension).
 * @param deg_rt Runtime degree (used if DegCT == 0).
 * @return Evaluated polynomial value(s).
 */
template<std::size_t DegCT = 0, bool SIMD = true, typename OutT, typename InVec, typename Mdspan>
PF_ALWAYS_INLINE constexpr OutT horner(const InVec &x, const Mdspan &coeffs, int deg_rt);

//------------------------------------------------------------------------------
// detail namespace
//------------------------------------------------------------------------------
namespace detail {

// helper to invoke coeffs(idx[0],…,idx[D‑1], d) in C++17
template<std::size_t D, typename MdspanType, std::size_t... Is>
PF_ALWAYS_INLINE constexpr auto call_coeffs_impl(const MdspanType &coeffs, const std::array<std::size_t, D> &idx,
                                                 std::index_sequence<Is...>, std::size_t d) noexcept {
    return coeffs(idx[Is]..., d);
}

template<std::size_t D, typename MdspanType>
PF_ALWAYS_INLINE constexpr auto call_coeffs(const MdspanType &coeffs, const std::array<std::size_t, D> &idx,
                                            std::size_t d) noexcept {
    return call_coeffs_impl<D, MdspanType>(coeffs, idx, std::make_index_sequence<D>{}, d);
}

//------------------------------------------------------------------------------
// detail::horner_nd_impl — tag-dispatched scalar / SIMD overloads
//------------------------------------------------------------------------------

struct scalar_t {};
struct simd_t {};

// Scalar overload — fully constexpr (no xsimd types)
template<std::size_t Level, std::size_t Dim, std::size_t DegCT, typename OutT, typename InVec, typename Mdspan>
PF_ALWAYS_INLINE constexpr OutT horner_nd_impl(scalar_t, const InVec &x, const Mdspan &coeffs,
                                               std::array<std::size_t, Dim> &idx, int deg_rt) {
    constexpr std::size_t axis = Dim - Level;
    constexpr std::size_t OUT = std::tuple_size_v<OutT>;

    OutT res{};

    auto step = [&](std::size_t k) {
        idx[axis] = k;

        OutT inner{};
        if constexpr (Level > 1) {
            inner = horner_nd_impl<Level - 1, Dim, DegCT, OutT>(scalar_t{}, x, coeffs, idx, deg_rt);
        } else {
            poet::static_for<OUT>([&](auto i) { inner[i] = detail::call_coeffs<Dim>(coeffs, idx, i); });
        }

        poet::static_for<OUT>([&](auto i) { res[i] = std::fma(res[i], x[axis], inner[i]); });
    };

    if constexpr (DegCT != 0) {
        poet::static_for<DegCT>([&](auto k) { step(k); });
    } else {
        for (int k = 0; k < deg_rt; ++k) step(static_cast<std::size_t>(k));
    }
    return res;
}

// SIMD overload — uses xsimd batch operations
template<std::size_t Level, std::size_t Dim, std::size_t DegCT, typename OutT, typename InVec, typename Mdspan>
PF_ALWAYS_INLINE OutT horner_nd_impl(simd_t, const InVec &x, const Mdspan &coeffs, std::array<std::size_t, Dim> &idx,
                                     int deg_rt) {
    constexpr std::size_t axis = Dim - Level;
    constexpr std::size_t OUT = std::tuple_size_v<OutT>;

    using Scalar = typename OutT::value_type;
    using batch = xsimd::make_sized_batch_t<Scalar, optimal_simd_width<Scalar, OUT>()>;
    alignas(batch::arch_type::alignment()) OutT res{0};
    const batch x_vec(x[axis]);

    auto step = [&](std::size_t k) {
        idx[axis] = k;
        alignas(batch::arch_type::alignment()) OutT inner{};
        if constexpr (Level > 1) {
            inner = horner_nd_impl<Level - 1, Dim, DegCT, OutT>(simd_t{}, x, coeffs, idx, deg_rt);
        } else {
            poet::static_for<OUT>([&](auto i) { inner[i] = call_coeffs<Dim>(coeffs, idx, i); });
        }

        poet::static_for<0, OUT, batch::size>([&](auto i) {
            if constexpr (i + batch::size <= OUT) {
                detail::fma(batch::load_aligned(res.data() + i), x_vec, batch::load_aligned(inner.data() + i))
                    .store_aligned(res.data() + i);
            } else {
                poet::static_for<i, OUT>([&](auto j) { res[j] = detail::fma(res[j], Scalar(x[axis]), inner[j]); });
            }
        });
    };

    if constexpr (DegCT != 0) {
        poet::static_for<DegCT>([&](auto k) { step(k); });
    } else {
        for (int k = 0; k < deg_rt; ++k) step(static_cast<std::size_t>(k));
    }
    return res;
}

} // namespace detail
// Compensated Horner (scalar only, for fitting accuracy)
// Tracks rounding errors at each FMA step using two_prod and two_sum,
// effectively doubling working precision: relative error O(u + n*u^2*cond)
// instead of O(n*u*cond). Only used during fitting (refinement residuals).

// Real overload (general OutputType)
template<std::size_t N_total = 0, typename OutputType, typename InputType>
PF_ALWAYS_INLINE constexpr OutputType compensated_horner(InputType x, const OutputType *c_ptr,
                                                         std::size_t c_size = 0) noexcept;

// Complex overload (more specialized, preferred by partial ordering)
template<std::size_t N_total = 0, typename T, typename InputType>
PF_ALWAYS_INLINE constexpr std::complex<T> compensated_horner(InputType x, const std::complex<T> *c_ptr,
                                                              std::size_t c_size = 0) noexcept;

} // namespace poly_eval

//------------------------------------------------------------------------------
// Implementation of Public API
//------------------------------------------------------------------------------

namespace poly_eval {

//------------------------------------------------------------------------------
// Horner (scalar, one-point)
//------------------------------------------------------------------------------

template<std::size_t N_total, typename OutputType, typename InputType>
PF_ALWAYS_INLINE constexpr OutputType horner(const InputType x, const OutputType *c_ptr,
                                             const std::size_t c_size) noexcept {
    if constexpr (N_total != 0) {
        // CT path: full unroll (N is known, serial chain)
        // GCC false positive: deep static_for inlining confuses -Wmaybe-uninitialized
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
        OutputType acc = c_ptr[0];
        poet::static_for<1, N_total>([&](auto i) { acc = detail::fma(acc, x, c_ptr[i]); });
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
        return acc;
    } else {
        OutputType acc = c_ptr[0];
        for (std::size_t k = 1; k < c_size; ++k) {
            acc = detail::fma(acc, x, c_ptr[k]);
        }
        return acc;
    }
}

//------------------------------------------------------------------------------
// Compensated Horner (scalar, one-point)
//------------------------------------------------------------------------------

// Real overload
template<std::size_t N_total, typename OutputType, typename InputType>
PF_ALWAYS_INLINE constexpr OutputType compensated_horner(const InputType x, const OutputType *c_ptr,
                                                         const std::size_t c_size) noexcept {
    namespace eft = polyfit::internal::helpers::eft;
    const std::size_t n = N_total != 0 ? N_total : c_size;
    OutputType acc = c_ptr[0];
    OutputType comp = OutputType(0);
    for (std::size_t k = 1; k < n; ++k) {
        auto [p, pi] = eft::two_prod(acc, OutputType(x));
        auto [s, sigma] = eft::two_sum(p, c_ptr[k]);
        comp = detail::fma(comp, OutputType(x), pi + sigma);
        acc = s;
    }
    return acc + comp;
}

// Complex overload: component-wise compensation
template<std::size_t N_total, typename T, typename InputType>
PF_ALWAYS_INLINE constexpr std::complex<T> compensated_horner(const InputType x, const std::complex<T> *c_ptr,
                                                              const std::size_t c_size) noexcept {
    namespace eft = polyfit::internal::helpers::eft;
    const std::size_t n = N_total != 0 ? N_total : c_size;
    T acc_re = c_ptr[0].real(), acc_im = c_ptr[0].imag();
    T comp_re = T(0), comp_im = T(0);
    const T xv = T(x);
    for (std::size_t k = 1; k < n; ++k) {
        auto [pr, pi_r] = eft::two_prod(acc_re, xv);
        auto [sr, sig_r] = eft::two_sum(pr, c_ptr[k].real());
        comp_re = detail::fma(comp_re, xv, pi_r + sig_r);
        acc_re = sr;

        auto [pi, pi_i] = eft::two_prod(acc_im, xv);
        auto [si, sig_i] = eft::two_sum(pi, c_ptr[k].imag());
        comp_im = detail::fma(comp_im, xv, pi_i + sig_i);
        acc_im = si;
    }
    return {acc_re + comp_re, acc_im + comp_im};
}

//------------------------------------------------------------------------------
// SIMD Horner (coeffs reversed)
//------------------------------------------------------------------------------

template<std::size_t N_monomials, bool pts_aligned, bool out_aligned, int UNROLL, typename InputType,
         typename OutputType, typename MapFunc>
PF_ALWAYS_INLINE constexpr void horner(const InputType *pts, OutputType *out, std::size_t num_points,
                                       const OutputType *monomials, std::size_t monomials_size,
                                       const MapFunc map_func) noexcept {
    PF_C23STATIC constexpr auto simd_size = xsimd::batch<InputType>::size;
    PF_C23STATIC constexpr auto UF = UNROLL > 0 ? UNROLL : detail::optimal_horner_uf<OutputType>();
    PF_C23STATIC constexpr auto block = simd_size * UF;

    using pts_mode = std::conditional_t<pts_aligned, xsimd::aligned_mode, xsimd::unaligned_mode>;
    using out_mode = std::conditional_t<out_aligned, xsimd::aligned_mode, xsimd::unaligned_mode>;

    std::size_t i = 0;
    if (num_points >= block)
        for (const std::size_t limit = num_points - block; i <= limit; i += block) {
            xsimd::batch<InputType> pt_batches[UF];
            xsimd::batch<OutputType> acc_batches[UF];

            // Load and init with highest-degree term
            poet::static_for<UF>([&](auto j) {
                pt_batches[j] = map_func(xsimd::load(pts + i + j * simd_size, pts_mode{}));
                acc_batches[j] = xsimd::batch<OutputType>(monomials[0]);
            });

            // Horner steps
            if constexpr (N_monomials != 0) {
                // CT path: full unroll (N is known)
                poet::static_for<1, N_monomials>([&](auto k) {
                    const auto ck = xsimd::batch<OutputType>(monomials[k]);
                    poet::static_for<UF>(
                        [&](auto j) { acc_batches[j] = detail::fma(acc_batches[j], pt_batches[j], ck); });
                });
            } else {
                for (std::size_t k = 1; k < monomials_size; ++k) {
                    const auto ck = xsimd::batch<OutputType>(monomials[k]);
                    poet::static_for<UF>(
                        [&](auto j) { acc_batches[j] = detail::fma(acc_batches[j], pt_batches[j], ck); });
                }
            }

            // Store results
            poet::static_for<UF>([&](auto j) { acc_batches[j].store(out + i + j * simd_size, out_mode{}); });
        }

    // Middle: single-batch SIMD for remaining full batches
    for (; i + simd_size <= num_points; i += simd_size) {
        auto pt = map_func(xsimd::load(pts + i, pts_mode{}));
        auto acc = xsimd::batch<OutputType>(monomials[0]);
        if constexpr (N_monomials != 0) {
            poet::static_for<1, N_monomials>(
                [&](auto k) { acc = detail::fma(acc, pt, xsimd::batch<OutputType>(monomials[k])); });
        } else {
            for (std::size_t k = 1; k < monomials_size; ++k) {
                acc = detail::fma(acc, pt, xsimd::batch<OutputType>(monomials[k]));
            }
        }
        acc.store(out + i, out_mode{});
    }

    // Final scalar remainder (< simd_size points)
    for (; i < num_points; ++i) {
        out[i] = horner<N_monomials>(map_func(pts[i]), monomials, monomials_size);
    }
}

//------------------------------------------------------------------------------
// horner_many
//------------------------------------------------------------------------------

template<std::size_t M_total, std::size_t N_total, bool scaling, typename OutputType, typename InputType>
PF_ALWAYS_INLINE constexpr void horner_many(const InputType xin, const OutputType *coeffs, OutputType *out,
                                            const std::size_t M, const std::size_t N, const InputType *low,
                                            const InputType *high) noexcept {
    const std::size_t m_lim = M_total ? M_total : M;
    const std::size_t n_lim = N_total ? N_total : N;

    if constexpr (M_total != 0) {
        poet::static_for<M_total>([&](auto m) {
            const auto xm = scaling ? (((InputType{2} * xin) - high[m]) * low[m]) : xin;
            out[m] = horner<N_total>(xm, coeffs + (m * n_lim), n_lim);
        });
    } else {
        for (std::size_t m = 0; m < m_lim; ++m) {
            const auto xm = scaling ? (((InputType{2} * xin) - high[m]) * low[m]) : xin;
            out[m] = horner<N_total>(xm, coeffs + (m * n_lim), n_lim);
        }
    }
}

//------------------------------------------------------------------------------
// horner_transposed — shared SIMD core + public dispatch
//------------------------------------------------------------------------------

namespace detail {

template<std::size_t N_total, std::size_t simd_width, bool is_aligned, typename Out, typename AccArr, typename XvecArr,
         typename ForChunks>
PF_ALWAYS_INLINE void horner_transposed_simd_core(const Out *coeffs, Out *out, std::size_t stride, std::size_t n_lim,
                                                  AccArr &acc, XvecArr &xvec, ForChunks for_chunks) noexcept {
    using batch_out = xsimd::make_sized_batch_t<Out, simd_width>;
    using mode = std::conditional_t<is_aligned, xsimd::aligned_mode, xsimd::unaligned_mode>;

    // Init accumulators from row 0
    for_chunks([&](auto ci) { acc[ci] = batch_out::load(coeffs + std::size_t(ci) * simd_width, mode{}); });

    // Horner steps
    auto step = [&](std::size_t k) {
        const auto *col = coeffs + k * stride;
        for_chunks([&](auto ci) {
            acc[ci] = detail::fma(acc[ci], xvec[ci], batch_out::load(col + std::size_t(ci) * simd_width, mode{}));
        });
    };

    if constexpr (N_total != 0)
        poet::static_for<1, N_total>([&](auto k) { step(k); });
    else
        for (std::size_t k = 1; k < n_lim; ++k) step(k);

    // Store results
    for_chunks([&](auto ci) { acc[ci].store(out + std::size_t(ci) * simd_width, mode{}); });
}

} // namespace detail

template<std::size_t M_total, std::size_t N_total, std::size_t simd_width, bool aligned, typename Out, typename In>
PF_ALWAYS_INLINE constexpr void horner_transposed(const In *xin, const Out *coeffs, Out *out, const std::size_t M,
                                                  const std::size_t N) noexcept {
    constexpr bool has_Mt = (M_total != 0);
    constexpr bool has_Nt = (N_total != 0);
    const std::size_t m_lim = has_Mt ? M_total : M;
    const std::size_t n_lim = has_Nt ? N_total : N;
    const std::size_t stride = m_lim;

    if constexpr (simd_width > 0) {
        using batch_in = xsimd::make_sized_batch_t<In, simd_width>;
        using mode = std::conditional_t<aligned, xsimd::aligned_mode, xsimd::unaligned_mode>;

        if constexpr (has_Mt) {
            constexpr std::size_t C = M_total / simd_width;
            std::array<xsimd::make_sized_batch_t<Out, simd_width>, C> acc;
            std::array<batch_in, C> xvec;
            poet::static_for<C>(
                [&](auto ci) { xvec[ci] = batch_in::load(xin + std::size_t(ci) * simd_width, mode{}); });
            detail::horner_transposed_simd_core<N_total, simd_width, aligned>(
                coeffs, out, stride, n_lim, acc, xvec,
                [](auto body) { poet::static_for<C>([&](auto ci) { body(ci); }); });
        } else {
            using batch_out = xsimd::make_sized_batch_t<Out, simd_width>;
            constexpr std::size_t max_chunks = 8; // supports up to 8*simd_width points
            const std::size_t C = m_lim / simd_width;
            assert(C <= max_chunks && "horner_transposed: C exceeds stack buffer capacity");
            batch_out acc[max_chunks];
            batch_in xvec[max_chunks];
            for (std::size_t ci = 0; ci < C; ++ci) xvec[ci] = batch_in::load(xin + ci * simd_width, mode{});
            detail::horner_transposed_simd_core<N_total, simd_width, aligned>(coeffs, out, stride, n_lim, acc, xvec,
                                                                              [C](auto body) {
                                                                                  for (std::size_t ci = 0; ci < C; ++ci)
                                                                                      body(ci);
                                                                              });
        }
    } else {
        // Scalar path
        if constexpr (has_Mt) {
            poet::static_for<M_total>([&](auto i) { out[i] = coeffs[i]; });
        } else {
            for (std::size_t i = 0; i < m_lim; ++i) out[i] = coeffs[i];
        }

        const auto step = [&](const std::size_t k) noexcept {
            const auto col = coeffs + (k * stride);
            if constexpr (has_Mt) {
                poet::static_for<M_total>([&](auto i) { out[i] = detail::fma(out[i], xin[i], col[i]); });
            } else {
                for (std::size_t i = 0; i < m_lim; ++i) {
                    out[i] = detail::fma(out[i], xin[i], col[i]);
                }
            }
        };

        if constexpr (has_Nt) {
            poet::static_for<1, N_total>([&](auto k) { step(k); });
        } else {
            for (std::size_t k = 1; k < n_lim; ++k) step(k);
        }
    }
}

//------------------------------------------------------------------------------
// ND Horner front-end
//------------------------------------------------------------------------------

template<std::size_t DegCT, bool SIMD, typename OutT, typename InVec, typename Mdspan>
PF_ALWAYS_INLINE constexpr OutT horner(const InVec &x, const Mdspan &coeffs, int deg_rt) {
    constexpr std::size_t Dim = Mdspan::rank() - 1;
    std::array<std::size_t, Dim> idx{};
    return detail::horner_nd_impl<Dim, Dim, DegCT, OutT>(std::conditional_t<SIMD, detail::simd_t, detail::scalar_t>{},
                                                         x, coeffs, idx, deg_rt);
}

} // namespace poly_eval
