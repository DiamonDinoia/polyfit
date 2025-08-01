#pragma once

#include "macros.h"
#include "utils.h"

#include <cassert>
#include <cstddef>
#include <xsimd/xsimd.hpp>

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
template <std::size_t N_total = 0, typename OutputType, typename InputType>
ALWAYS_INLINE constexpr OutputType horner(InputType x, const OutputType *c_ptr, std::size_t c_size = 0) noexcept;

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
template <std::size_t N_monomials = 0, bool pts_aligned = false, bool out_aligned = false, int UNROLL = 0,
          typename InputType, typename OutputType, typename MapFunc>
ALWAYS_INLINE constexpr void horner(
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
template <std::size_t M_total = 0, std::size_t N_total = 0, bool scaling = false, typename OutputType,
          typename InputType>
ALWAYS_INLINE constexpr void horner_many(InputType x, const OutputType *coeffs, OutputType *out, std::size_t M = 0,
                                         std::size_t N = 0, const InputType *low = nullptr,
                                         const InputType *hi = nullptr) noexcept;

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
template <std::size_t M_total = 0, std::size_t N_total = 0, std::size_t simd_width = 0, typename Out, typename In>
ALWAYS_INLINE constexpr void horner_transposed(const In *x, const Out *c, Out *out, std::size_t M = 0,
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
template <std::size_t DegCT = 0, bool SIMD = true, typename OutT, typename InVec, typename Mdspan>
ALWAYS_INLINE constexpr OutT horner(const InVec &x, const Mdspan &coeffs, int deg_rt);

//------------------------------------------------------------------------------
// detail namespace
//------------------------------------------------------------------------------
namespace detail {

// helper to invoke coeffs(idx[0],…,idx[D‑1], d) in C++17
template <std::size_t D, typename MdspanType, std::size_t... Is>
ALWAYS_INLINE constexpr auto call_coeffs_impl(const MdspanType &coeffs, const std::array<std::size_t, D> &idx,
                                              std::index_sequence<Is...>, std::size_t d) noexcept {
    return coeffs(idx[Is]..., d);
}

template <std::size_t D, typename MdspanType>
ALWAYS_INLINE constexpr auto call_coeffs(const MdspanType &coeffs, const std::array<std::size_t, D> &idx,
                                         std::size_t d) noexcept {
    return call_coeffs_impl<D, MdspanType>(coeffs, idx, std::make_index_sequence<D>{}, d);
}

//------------------------------------------------------------------------------
// detail::horner_nd_impl
//------------------------------------------------------------------------------

template <std::size_t Level, std::size_t Dim, std::size_t DegCT, bool SIMD, typename OutT, typename InVec,
          typename Mdspan>
ALWAYS_INLINE constexpr OutT horner_nd_impl(const InVec &x, const Mdspan &coeffs, std::array<std::size_t, Dim> &idx,
                                            const int deg_rt) {
    if constexpr (!SIMD) {
        constexpr std::size_t axis = Dim - Level; // current axis
        constexpr std::size_t OUT = std::tuple_size_v<OutT>;
        const int deg = DegCT ? int(DegCT) : deg_rt;

        OutT res{}; // zero-init

        /* ---------- descend k = 0 … deg-1 (highest → lowest power) ---------- */
        for (int k = 0; k < deg; ++k) {
            idx[axis] = static_cast<std::size_t>(k);

            OutT inner{};
            if constexpr (Level > 1) {
                inner = horner_nd_impl<Level - 1, Dim, DegCT, SIMD, OutT>(x, coeffs, idx, deg_rt);
            } else {
                //  last spatial level → read the coefficient vector
                detail::unroll_loop<OUT>([&](const auto i) {
                    constexpr auto d = decltype(i)::value;
                    inner[d] = [&]<std::size_t... Is>(std::index_sequence<Is...>) {
                        return coeffs(idx[Is]..., d); // Dim indices + output
                    }(std::make_index_sequence<Dim>{});
                });
            }

            /* ---------- Horner accumulate along this axis ---------- */
            detail::unroll_loop<OUT>([&](const auto i) {
                constexpr auto d = decltype(i)::value;
                res[d] = std::fma(res[d], x[axis], inner[d]);
            });
        }
        return res;
    }
    constexpr std::size_t axis = Dim - Level;
    constexpr std::size_t OUT = std::tuple_size_v<OutT>;
    const int deg = DegCT ? static_cast<int>(DegCT) : deg_rt;

    using batch =
        xsimd::make_sized_batch_t<typename OutT::value_type, optimal_simd_width<typename OutT::value_type, OUT>()>;
    alignas(batch::arch_type::alignment()) OutT res{0};
    const batch x_vec(x[axis]);

    for (int k = 0; k < deg; ++k) {
        idx[axis] = static_cast<std::size_t>(k);
        alignas(batch::arch_type::alignment()) OutT inner{};
        if constexpr (Level > 1) {
            inner = horner_nd_impl<Level - 1, Dim, DegCT, SIMD, OutT>(x, coeffs, idx, deg_rt);
        } else {
            detail::unroll_loop<OUT>([&]([[maybe_unused]] const auto I) {
                constexpr auto d = decltype(I)::value;
                inner[d] = call_coeffs<Dim>(coeffs, idx, d);
            });
        }

        detail::unroll_loop<0, OUT, /*Inc=*/batch::size>([&]([[maybe_unused]] const auto I) {
            constexpr auto d = decltype(I)::value;
            if constexpr (d + batch::size <= OUT) {
                auto in = batch::load_aligned(inner.data() + d);
                auto r = batch::load_aligned(res.data() + d);
                detail::fma(r, x_vec, in).store_aligned(res.data() + d);
            } else {
                detail::unroll_loop<d, OUT>([&]([[maybe_unused]] const auto J) {
                    constexpr auto last = decltype(J)::value;
                    res[last] = detail::fma(res[last], typename OutT::value_type(x[axis]), inner[last]);
                });
            }
        });
    }
    return res;
}

} // namespace detail
} // namespace poly_eval

//------------------------------------------------------------------------------
// Implementation of Public API
//------------------------------------------------------------------------------

namespace poly_eval {

//------------------------------------------------------------------------------
// Horner (scalar, one-point)
//------------------------------------------------------------------------------

template <std::size_t N_total, typename OutputType, typename InputType>
ALWAYS_INLINE constexpr OutputType horner(const InputType x, const OutputType *c_ptr,
                                          const std::size_t c_size) noexcept {
    if constexpr (N_total != 0) {
        // Compile-time unrolled Horner on reversed array
        OutputType acc = c_ptr[0];
        detail::unroll_loop<1, N_total>([&]([[maybe_unused]] const auto I) {
            constexpr auto k = decltype(I)::value;
            acc = detail::fma(acc, x, c_ptr[k]);
        });
        return acc;
    } else {
        // Runtime iterative Horner
        OutputType acc = c_ptr[0];
        for (std::size_t k = 1; k < c_size; ++k) {
            acc = detail::fma(acc, x, c_ptr[k]);
        }
        return acc;
    }
}

//------------------------------------------------------------------------------
// SIMD Horner (coeffs reversed)
//------------------------------------------------------------------------------

template <std::size_t N_monomials, bool pts_aligned, bool out_aligned, int UNROLL, typename InputType,
          typename OutputType, typename MapFunc>
ALWAYS_INLINE constexpr void horner(const InputType *pts, OutputType *out, std::size_t num_points,
                                    const OutputType *monomials, std::size_t monomials_size,
                                    const MapFunc map_func) noexcept {
    C23STATIC constexpr auto simd_size = xsimd::batch<InputType>::size;
    C23STATIC constexpr auto OuterUnrollFactor = UNROLL > 0 ? UNROLL : 32 / sizeof(OutputType);
    C23STATIC constexpr auto block = simd_size * OuterUnrollFactor;

    using pts_mode = std::conditional_t<pts_aligned, xsimd::aligned_mode, xsimd::unaligned_mode>;
    using out_mode = std::conditional_t<out_aligned, xsimd::aligned_mode, xsimd::unaligned_mode>;
    const std::size_t trunc = num_points & (-block);

    for (std::size_t i = 0; i < trunc; i += block) {
        xsimd::batch<InputType> pt_batches[OuterUnrollFactor];
        xsimd::batch<OutputType> acc_batches[OuterUnrollFactor];

        // Load and init with highest-degree term
        detail::unroll_loop<OuterUnrollFactor>([&]([[maybe_unused]] const auto I) {
            constexpr auto j = decltype(I)::value;
            pt_batches[j] = map_func(xsimd::load(pts + i + j * simd_size, pts_mode{}));
            acc_batches[j] = xsimd::batch<OutputType>(monomials[0]);
        });

        // Horner steps
        if constexpr (N_monomials != 0) {
            detail::unroll_loop<1, N_monomials>([&]([[maybe_unused]] const auto I) {
                constexpr auto k = decltype(I)::value;
                detail::unroll_loop<OuterUnrollFactor>([&]([[maybe_unused]] const auto J) {
                    constexpr auto j = decltype(J)::value;
                    acc_batches[j] = detail::fma(acc_batches[j], pt_batches[j], xsimd::batch<OutputType>(monomials[k]));
                });
            });
        } else {
            for (std::size_t k = 1; k < monomials_size; ++k) {
                detail::unroll_loop<OuterUnrollFactor>([&]([[maybe_unused]] const auto I) {
                    constexpr auto j = decltype(I)::value;
                    acc_batches[j] = detail::fma(acc_batches[j], pt_batches[j], xsimd::batch<OutputType>(monomials[k]));
                });
            }
        }

        // Store results
        detail::unroll_loop<OuterUnrollFactor>([&]([[maybe_unused]] const auto I) {
            constexpr auto j = decltype(I)::value;
            acc_batches[j].store(out + i + j * simd_size, out_mode{});
        });
    }

    // Remainder
    ASSUME((num_points - trunc) < block);
    for (std::size_t idx = trunc; idx < num_points; ++idx) {
        out[idx] = horner<N_monomials>(map_func(pts[idx]), monomials, monomials_size);
    }
}

//------------------------------------------------------------------------------
// horner_many
//------------------------------------------------------------------------------

template <std::size_t M_total, std::size_t N_total, bool scaling, typename OutputType, typename InputType>
ALWAYS_INLINE constexpr void horner_many(const InputType x, const OutputType *coeffs, OutputType *out,
                                         const std::size_t M, const std::size_t N, const InputType *low,
                                         const InputType *hi) noexcept {
    const std::size_t m_lim = M_total ? M_total : M;
    const std::size_t n_lim = N_total ? N_total : N;

    if constexpr (M_total != 0) {
        detail::unroll_loop<M_total>([&]([[maybe_unused]] const auto I) {
            constexpr auto m = decltype(I)::value;
            const auto xm = scaling ? (InputType{2} * x - hi[m]) * low[m] : x;
            out[m] = horner<N_total>(xm, coeffs + m * n_lim, n_lim);
        });
    } else {
        for (std::size_t m = 0; m < m_lim; ++m) {
            const auto xm = scaling ? (InputType{2} * x - hi[m]) * low[m] : x;
            out[m] = horner<N_total>(xm, coeffs + m * n_lim, n_lim);
        }
    }
}

//------------------------------------------------------------------------------
// horner_transposed
//------------------------------------------------------------------------------

template <std::size_t M_total, std::size_t N_total, std::size_t simd_width, typename Out, typename In>
ALWAYS_INLINE constexpr void horner_transposed(const In *x, const Out *c, Out *out, const std::size_t M,
                                               const std::size_t N) noexcept {
    constexpr bool has_Mt = (M_total != 0);
    constexpr bool has_Nt = (N_total != 0);
    constexpr bool do_simd = (simd_width > 0);
    const std::size_t m_lim = has_Mt ? M_total : M;
    const std::size_t n_lim = has_Nt ? N_total : N;
    const std::size_t stride = m_lim;

    if constexpr (do_simd) {
        using batch_in = xsimd::make_sized_batch_t<In, simd_width>;
        using batch_out = xsimd::make_sized_batch_t<Out, simd_width>;
        if constexpr (has_Mt) {
            constexpr std::size_t C = M_total / simd_width;
            std::array<batch_out, C> acc;
            std::array<batch_in, C> xvec;

            detail::unroll_loop<C>([&]([[maybe_unused]] const auto I) {
                constexpr auto ci = decltype(I)::value;
                constexpr auto base = ci * simd_width;
                xvec[ci] = batch_in::load_unaligned(x + base);
                acc[ci] = batch_out::load_unaligned(c + base);
            });

            if constexpr (has_Nt) {
                detail::unroll_loop<0, N_total>([&]([[maybe_unused]] const auto I) {
                    constexpr auto k = decltype(I)::value;
                    if constexpr (k > 0) {
                        const auto col = c + k * stride;
                        detail::unroll_loop<C>([&]([[maybe_unused]] const auto J) {
                            constexpr auto ci = decltype(J)::value;
                            constexpr auto base = ci * simd_width;
                            batch_out ck = batch_out::load_unaligned(col + base);
                            acc[ci] = detail::fma(acc[ci], xvec[ci], ck);
                        });
                    }
                });
            } else {
                for (std::size_t k = 1; k < n_lim; ++k) {
                    const auto col = c + k * stride;
                    detail::unroll_loop<C>([&]([[maybe_unused]] const auto I) {
                        constexpr auto ci = decltype(I)::value;
                        constexpr auto base = ci * simd_width;
                        batch_out ck = batch_out::load_unaligned(col + base);
                        acc[ci] = detail::fma(acc[ci], xvec[ci], ck);
                    });
                }
            }

            detail::unroll_loop<C>([&]([[maybe_unused]] const auto I) {
                constexpr auto ci = decltype(I)::value;
                constexpr auto base = ci * simd_width;
                acc[ci].store_unaligned(out + base);
            });
        } else {
            const std::size_t chunks = m_lim / simd_width;
            std::vector<batch_out> acc(chunks);
            std::vector<batch_in> xvec(chunks);

            for (std::size_t ci = 0; ci < chunks; ++ci) {
                const auto base = ci * simd_width;
                xvec[ci] = batch_in::load_unaligned(x + base);
                acc[ci] = batch_out::load_unaligned(c + base);
            }
            for (std::size_t k = 1; k < n_lim; ++k) {
                const auto col = c + k * stride;
                for (std::size_t ci = 0; ci < chunks; ++ci) {
                    const auto base = ci * simd_width;
                    batch_out ck = batch_out::load_unaligned(col + base);
                    acc[ci] = detail::fma(acc[ci], xvec[ci], ck);
                }
            }
            for (std::size_t ci = 0; ci < chunks; ++ci) {
                const auto base = ci * simd_width;
                acc[ci].store_unaligned(out + base);
            }
        }
    } else {
        // Scalar path
        if constexpr (has_Mt) {
            detail::unroll_loop<M_total>([&]([[maybe_unused]] const auto I) {
                constexpr auto i = decltype(I)::value;
                out[i] = c[i];
            });
        } else {
            for (std::size_t i = 0; i < m_lim; ++i)
                out[i] = c[i];
        }

        const auto step = [&](const std::size_t k) noexcept {
            const auto col = c + k * stride;
            if constexpr (has_Mt) {
                detail::unroll_loop<M_total>([&]([[maybe_unused]] const auto I) {
                    constexpr auto i = decltype(I)::value;
                    out[i] = detail::fma(out[i], x[i], col[i]);
                });
            } else {
                for (std::size_t i = 0; i < m_lim; ++i)
                    out[i] = detail::fma(out[i], x[i], col[i]);
            }
        };

        if constexpr (has_Nt) {
            detail::unroll_loop<0, N_total>([&]([[maybe_unused]] const auto I) {
                constexpr auto k = decltype(I)::value;
                if constexpr (k > 0)
                    step(k);
            });
        } else {
            for (std::size_t k = 1; k < n_lim; ++k)
                step(k);
        }
    }
}

//------------------------------------------------------------------------------
// ND Horner front-end
//------------------------------------------------------------------------------

template <std::size_t DegCT, bool SIMD, typename OutT, typename InVec, typename Mdspan>
ALWAYS_INLINE constexpr OutT horner(const InVec &x, const Mdspan &coeffs, int deg_rt) {
    constexpr std::size_t Dim = Mdspan::rank() - 1;
    std::array<std::size_t, Dim> idx{};
    return detail::horner_nd_impl<Dim, Dim, DegCT, SIMD, OutT>(x, coeffs, idx, deg_rt);
}

} // namespace poly_eval