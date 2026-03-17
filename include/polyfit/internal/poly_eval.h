#pragma once

#include "macros.h"
#include "simd_utils.h"
#include "utils.h"

#include <poet/poet.hpp>

#include <cassert>
#include <cstddef>

namespace poly_eval {

// Evaluate one polynomial at one point. Coefficients are in Horner order.
template<std::size_t NCOEFFS = 0, typename CoeffType, typename InputType>
PF_ALWAYS_INLINE constexpr CoeffType horner(InputType xin, const CoeffType *c_ptr, std::size_t c_size = 0) noexcept;

// Evaluate one polynomial across many points, optionally with SIMD and a point mapping step.
template<std::size_t NMONOMIALS = 0, bool PTS_ALIGNED = false, bool OUT_ALIGNED = false, int UNROLL = 0,
         typename InputType, typename OutputType, typename MapFunc>
PF_ALWAYS_INLINE constexpr void horner(
    const InputType *pts, OutputType *out, std::size_t num_points, const OutputType *monomials,
    std::size_t monomials_size, MapFunc map_func = [](auto v) { return v; }) noexcept;

// Evaluate several 1D polynomials at one point.
template<std::size_t MCOUNT = 0, std::size_t NCOEFFS = 0, typename OutputType, typename InputType>
PF_ALWAYS_INLINE constexpr void horner_many(InputType xin, const OutputType *coeffs, OutputType *out, std::size_t M = 0,
                                            std::size_t N = 0) noexcept;

// Evaluate many polynomials across many points in transposed coefficient layout.
template<std::size_t MCOUNT = 0, std::size_t NCOEFFS = 0, std::size_t SIMD_WIDTH = 0, bool ALIGNED = false,
         typename Out, typename In>
PF_ALWAYS_INLINE constexpr void horner_transposed(const In *xin, const Out *coeffs, Out *out, std::size_t M = 0,
                                                  std::size_t N = 0) noexcept;

// Evaluate an ND polynomial stored in an mdspan-backed coefficient tensor.
template<std::size_t NCOEFFS = 0, bool SIMD = true, typename OutT, typename InVec, typename Mdspan>
PF_ALWAYS_INLINE constexpr OutT horner(const InVec &x, const Mdspan &coeffs, int nCoeffsRt);

namespace detail {

template<std::size_t D, typename MdspanType, std::size_t... Is>
PF_ALWAYS_INLINE constexpr auto coeffAtImpl(const MdspanType &coeffs, const std::array<std::size_t, D> &idx,
                                            std::index_sequence<Is...>, std::size_t outputIndex) noexcept {
    return coeffs[std::array<std::size_t, D + 1>{idx[Is]..., outputIndex}];
}

template<std::size_t D, typename MdspanType>
PF_ALWAYS_INLINE constexpr auto coeffAt(const MdspanType &coeffs, const std::array<std::size_t, D> &idx,
                                        std::size_t outputIndex) noexcept {
    return coeffAtImpl<D, MdspanType>(coeffs, idx, std::make_index_sequence<D>{}, outputIndex);
}

template<std::size_t NCOEFFS, class Step>
PF_ALWAYS_INLINE constexpr void forEachCoeff(int nCoeffsRt, Step &&step) {
    if constexpr (NCOEFFS != 0) {
        poet::static_for<NCOEFFS>([&](auto k) { step(std::size_t(k)); });
    } else {
        for (int k = 0; k < nCoeffsRt; ++k) step(static_cast<std::size_t>(k));
    }
}

struct ScalarEvalTag {};
struct SimdEvalTag {};

template<std::size_t Level, std::size_t Dim, std::size_t NCOEFFS, typename OutT, typename InVec, typename Mdspan>
PF_ALWAYS_INLINE constexpr OutT horner_nd_impl(ScalarEvalTag, const InVec &x, const Mdspan &coeffs,
                                               std::array<std::size_t, Dim> &idx, int nCoeffsRt) {
    constexpr std::size_t axis = Dim - Level;
    constexpr std::size_t OUT = std::tuple_size_v<OutT>;

    OutT res{};

    auto step = [&](std::size_t k) {
        idx[axis] = k;

        OutT inner{};
        if constexpr (Level > 1) {
            inner = horner_nd_impl<Level - 1, Dim, NCOEFFS, OutT>(ScalarEvalTag{}, x, coeffs, idx, nCoeffsRt);
        } else {
            poet::static_for<OUT>([&](auto i) { inner[i] = detail::coeffAt<Dim>(coeffs, idx, i); });
        }

        poet::static_for<OUT>([&](auto i) { res[i] = detail::fma(res[i], x[axis], inner[i]); });
    };

    forEachCoeff<NCOEFFS>(nCoeffsRt, step);
    return res;
}

template<std::size_t Level, std::size_t Dim, std::size_t NCOEFFS, typename OutT, typename InVec, typename Mdspan>
PF_ALWAYS_INLINE OutT horner_nd_impl(SimdEvalTag, const InVec &x, const Mdspan &coeffs,
                                     std::array<std::size_t, Dim> &idx, int nCoeffsRt) {
    constexpr std::size_t axis = Dim - Level;
    constexpr std::size_t OUT = std::tuple_size_v<OutT>;

    using Scalar = typename OutT::value_type;
    using batch = xsimd::make_sized_batch_t<Scalar, optimalSimdWidth<Scalar, OUT>()>;
    alignas(batch::arch_type::alignment()) OutT res{0};
    const batch x_vec(x[axis]);

    auto step = [&](std::size_t k) {
        idx[axis] = k;
        alignas(batch::arch_type::alignment()) OutT inner{};
        if constexpr (Level > 1) {
            inner = horner_nd_impl<Level - 1, Dim, NCOEFFS, OutT>(SimdEvalTag{}, x, coeffs, idx, nCoeffsRt);
        } else {
            poet::static_for<OUT>([&](auto i) { inner[i] = coeffAt<Dim>(coeffs, idx, i); });
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

    forEachCoeff<NCOEFFS>(nCoeffsRt, step);
    return res;
}

} // namespace detail

template<std::size_t NCOEFFS = 0, typename OutputType, typename InputType>
constexpr OutputType compensated_horner(InputType x, const OutputType *c_ptr, std::size_t c_size = 0) noexcept;

template<std::size_t NCOEFFS = 0, typename T, typename InputType>
constexpr std::complex<T> compensated_horner(InputType x, const std::complex<T> *c_ptr, std::size_t c_size = 0) noexcept;

} // namespace poly_eval

namespace poly_eval::detail {

template<std::size_t NCOEFFS, typename EvalType, typename CoeffType, typename InputType>
PF_ALWAYS_INLINE constexpr EvalType horner_impl(const InputType x, const CoeffType *c_ptr,
                                                const std::size_t c_size) noexcept {
    if constexpr (NCOEFFS != 0) {
        // CT path: full unroll (N is known, serial chain)
        // GCC false positive: deep static_for inlining confuses -Wmaybe-uninitialized
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
        EvalType acc = static_cast<EvalType>(c_ptr[0]);
        poet::static_for<1, NCOEFFS>([&](auto i) { acc = detail::fma(acc, x, static_cast<EvalType>(c_ptr[i])); });
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
        return acc;
    } else {
        EvalType acc = static_cast<EvalType>(c_ptr[0]);
        for (std::size_t k = 1; k < c_size; ++k) {
            acc = detail::fma(acc, x, static_cast<EvalType>(c_ptr[k]));
        }
        return acc;
    }
}

} // namespace poly_eval::detail

namespace poly_eval {

template<std::size_t NCOEFFS, typename CoeffType, typename InputType>
PF_ALWAYS_INLINE constexpr CoeffType horner(const InputType x, const CoeffType *c_ptr, const std::size_t c_size) noexcept {
    return detail::horner_impl<NCOEFFS, CoeffType>(x, c_ptr, c_size);
}

template<std::size_t NCOEFFS, typename OutputType, typename InputType>
constexpr OutputType compensated_horner(const InputType x, const OutputType *c_ptr, const std::size_t c_size) noexcept {
    namespace eft = polyfit::internal::helpers::eft;
    const std::size_t n = NCOEFFS != 0 ? NCOEFFS : c_size;
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
template<std::size_t NCOEFFS, typename T, typename InputType>
constexpr std::complex<T> compensated_horner(const InputType x, const std::complex<T> *c_ptr, const std::size_t c_size) noexcept {
    namespace eft = polyfit::internal::helpers::eft;
    const std::size_t n = NCOEFFS != 0 ? NCOEFFS : c_size;
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
// SIMD Horner
//------------------------------------------------------------------------------

template<std::size_t NMONOMIALS, bool PTS_ALIGNED, bool OUT_ALIGNED, int UNROLL, typename InputType,
         typename OutputType, typename MapFunc>
PF_ALWAYS_INLINE constexpr void horner(const InputType *pts, OutputType *out, std::size_t num_points,
                                       const OutputType *monomials, std::size_t monomials_size,
                                       const MapFunc map_func) noexcept {
    PF_STATIC_CONSTEXPR_LOCAL auto simd_size = xsimd::batch<InputType>::size;
    PF_STATIC_CONSTEXPR_LOCAL auto UF = UNROLL > 0 ? UNROLL : detail::optimalHornerUf<OutputType>();
    PF_STATIC_CONSTEXPR_LOCAL auto block = simd_size * UF;

    using pts_mode = std::conditional_t<PTS_ALIGNED, xsimd::aligned_mode, xsimd::unaligned_mode>;
    using out_mode = std::conditional_t<OUT_ALIGNED, xsimd::aligned_mode, xsimd::unaligned_mode>;

    const auto tile_end = detail::roundDown<block>(num_points);
    const auto simd_end = detail::roundDown<simd_size>(num_points);

    auto for_k = [&](auto step) {
        if constexpr (NMONOMIALS != 0)
            poet::static_for<1, NMONOMIALS>(step);
        else
            for (std::size_t k = 1; k < monomials_size; ++k) step(k);
    };

    for (std::size_t i = 0; i < tile_end; i += block) {
        xsimd::batch<InputType> pt_batches[UF];
        xsimd::batch<OutputType> acc_batches[UF];

        poet::static_for<UF>([&](auto j) {
            pt_batches[j] = map_func(xsimd::load(pts + i + j * simd_size, pts_mode{}));
            acc_batches[j] = xsimd::batch<OutputType>(monomials[0]);
        });

        for_k([&](auto k) {
            const auto ck = xsimd::batch<OutputType>(monomials[k]);
            poet::static_for<UF>(
                [&](auto j) { acc_batches[j] = detail::fma(acc_batches[j], pt_batches[j], ck); });
        });

        poet::static_for<UF>([&](auto j) { acc_batches[j].store(out + i + j * simd_size, out_mode{}); });
    }

    for (std::size_t i = tile_end; i < simd_end; i += simd_size) {
        auto pt = map_func(xsimd::load(pts + i, pts_mode{}));
        auto acc = xsimd::batch<OutputType>(monomials[0]);
        for_k([&](auto k) { acc = detail::fma(acc, pt, xsimd::batch<OutputType>(monomials[k])); });
        acc.store(out + i, out_mode{});
    }

    for (std::size_t i = simd_end; i < num_points; ++i)
        out[i] = horner<NMONOMIALS>(map_func(pts[i]), monomials, monomials_size);
}

//------------------------------------------------------------------------------
// horner_many
//------------------------------------------------------------------------------

template<std::size_t MCOUNT, std::size_t NCOEFFS, typename OutputType, typename InputType>
PF_ALWAYS_INLINE constexpr void horner_many(const InputType xin, const OutputType *coeffs, OutputType *out,
                                            const std::size_t M, const std::size_t N) noexcept {
    const std::size_t n_lim = NCOEFFS ? NCOEFFS : N;

    if constexpr (MCOUNT != 0) {
        poet::static_for<MCOUNT>([&](auto m) {
            out[m] = horner<NCOEFFS>(xin, coeffs + (m * n_lim), n_lim);
        });
    } else {
        const std::size_t m_lim = M;
        constexpr auto UF = detail::optimalHornerManyUf();
        poet::dynamic_for<UF>(m_lim, [&](auto /*lane*/, std::size_t m) {
            out[m] = horner<NCOEFFS>(xin, coeffs + (m * n_lim), n_lim);
        });
    }
}

//------------------------------------------------------------------------------
// horner_transposed — shared SIMD core + public dispatch
//------------------------------------------------------------------------------

namespace detail {

template<std::size_t NCOEFFS, std::size_t SIMD_WIDTH, bool IS_ALIGNED, typename Out, typename AccArr, typename XvecArr,
         typename ForChunks>
PF_ALWAYS_INLINE void horner_transposed_simd_core(const Out *coeffs, Out *out, std::size_t stride, std::size_t n_lim,
                                                  AccArr &acc, XvecArr &xvec, ForChunks for_chunks) noexcept {
    using batch_out = xsimd::make_sized_batch_t<Out, SIMD_WIDTH>;
    using mode = std::conditional_t<IS_ALIGNED, xsimd::aligned_mode, xsimd::unaligned_mode>;

    for_chunks([&](auto ci) { acc[ci] = batch_out::load(coeffs + std::size_t(ci) * SIMD_WIDTH, mode{}); });

    auto step = [&](std::size_t k) {
        const auto *col = coeffs + k * stride;
        for_chunks([&](auto ci) {
            acc[ci] = detail::fma(acc[ci], xvec[ci], batch_out::load(col + std::size_t(ci) * SIMD_WIDTH, mode{}));
        });
    };

    if constexpr (NCOEFFS != 0)
        poet::static_for<1, NCOEFFS>([&](auto k) { step(k); });
    else
        for (std::size_t k = 1; k < n_lim; ++k) step(k);

    for_chunks([&](auto ci) { acc[ci].store(out + std::size_t(ci) * SIMD_WIDTH, mode{}); });
}

} // namespace detail

template<std::size_t MCOUNT, std::size_t NCOEFFS, std::size_t SIMD_WIDTH, bool ALIGNED, typename Out, typename In>
PF_ALWAYS_INLINE constexpr void horner_transposed(const In *xin, const Out *coeffs, Out *out, const std::size_t M,
                                                  const std::size_t N) noexcept {
    constexpr bool has_Mt = (MCOUNT != 0);
    constexpr bool has_Nt = (NCOEFFS != 0);
    const std::size_t m_lim = has_Mt ? MCOUNT : M;
    const std::size_t n_lim = has_Nt ? NCOEFFS : N;
    const std::size_t stride = m_lim;

    if constexpr (SIMD_WIDTH > 0) {
        using batch_in = xsimd::make_sized_batch_t<In, SIMD_WIDTH>;
        using mode = std::conditional_t<ALIGNED, xsimd::aligned_mode, xsimd::unaligned_mode>;

        if constexpr (has_Mt) {
            constexpr std::size_t C = MCOUNT / SIMD_WIDTH;
            std::array<xsimd::make_sized_batch_t<Out, SIMD_WIDTH>, C> acc;
            std::array<batch_in, C> xvec;
            poet::static_for<C>(
                [&](auto ci) { xvec[ci] = batch_in::load(xin + std::size_t(ci) * SIMD_WIDTH, mode{}); });
            detail::horner_transposed_simd_core<NCOEFFS, SIMD_WIDTH, ALIGNED>(
                coeffs, out, stride, n_lim, acc, xvec,
                [](auto body) { poet::static_for<C>([&](auto ci) { body(ci); }); });
        } else {
            using batch_out = xsimd::make_sized_batch_t<Out, SIMD_WIDTH>;
            constexpr std::size_t max_chunks = 8;
            const std::size_t C = m_lim / SIMD_WIDTH;
            assert(C <= max_chunks && "horner_transposed: C exceeds stack buffer capacity");
            batch_out acc[max_chunks];
            batch_in xvec[max_chunks];
            for (std::size_t ci = 0; ci < C; ++ci) xvec[ci] = batch_in::load(xin + ci * SIMD_WIDTH, mode{});
            detail::horner_transposed_simd_core<NCOEFFS, SIMD_WIDTH, ALIGNED>(coeffs, out, stride, n_lim, acc, xvec,
                                                                              [C](auto body) {
                                                                                  for (std::size_t ci = 0; ci < C; ++ci)
                                                                                      body(ci);
                                                                              });
        }
    } else {
        if constexpr (has_Mt) {
            poet::static_for<MCOUNT>([&](auto i) { out[i] = coeffs[i]; });
        } else {
            for (std::size_t i = 0; i < m_lim; ++i) out[i] = coeffs[i];
        }

        const auto step = [&](const std::size_t k) noexcept {
            const auto col = coeffs + (k * stride);
            if constexpr (has_Mt) {
                poet::static_for<MCOUNT>([&](auto i) { out[i] = detail::fma(out[i], xin[i], col[i]); });
            } else {
                for (std::size_t i = 0; i < m_lim; ++i) {
                    out[i] = detail::fma(out[i], xin[i], col[i]);
                }
            }
        };

        if constexpr (has_Nt) {
            poet::static_for<1, NCOEFFS>([&](auto k) { step(k); });
        } else {
            for (std::size_t k = 1; k < n_lim; ++k) step(k);
        }
    }
}

//------------------------------------------------------------------------------
// ND Horner front-end
//------------------------------------------------------------------------------

template<std::size_t NCOEFFS, bool SIMD, typename OutT, typename InVec, typename Mdspan>
PF_ALWAYS_INLINE constexpr OutT horner(const InVec &x, const Mdspan &coeffs, int nCoeffsRt) {
    constexpr std::size_t Dim = Mdspan::rank() - 1;
    std::array<std::size_t, Dim> idx{};
    return detail::horner_nd_impl<Dim, Dim, NCOEFFS, OutT>(
        std::conditional_t<SIMD, detail::SimdEvalTag, detail::ScalarEvalTag>{}, x, coeffs, idx, nCoeffsRt);
}

} // namespace poly_eval
