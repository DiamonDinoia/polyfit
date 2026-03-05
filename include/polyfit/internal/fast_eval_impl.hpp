#pragma once

#include <cassert>
#include <stdexcept>
#include <xsimd/xsimd.hpp>

#include "helpers.h"
#include "macros.h"
#include "poly_eval.h"
#include "utils.h"

#include <poet/poet.hpp>

namespace poly_eval {
namespace detail {

PF_ALWAYS_INLINE constexpr int validate_positive_degree(const int n) {
    if (n <= 0) {
        throw std::invalid_argument("Runtime polynomial degree must be positive (n > 0)");
    }
    return n;
}

} // namespace detail

// -----------------------------------------------------------------------------
// FuncEval Implementation (Runtime)
// -----------------------------------------------------------------------------

template<class Func, std::size_t N_compile_time, std::size_t Iters_compile_time>
template<std::size_t CurrentN, typename>
PF_C20CONSTEXPR FuncEval<Func, N_compile_time, Iters_compile_time>::FuncEval(Func F, const int n, const InputType a,
                                                                             const InputType b, const InputType *pts)
    : low(InputType(1) / (b - a)), hi(b + a) {
    monomials.resize(static_cast<std::size_t>(detail::validate_positive_degree(n)));
    initialize_monomials(F, pts);
}

template<class Func, std::size_t N_compile_time, std::size_t Iters_compile_time>
template<std::size_t CurrentN, typename>
PF_C20CONSTEXPR FuncEval<Func, N_compile_time, Iters_compile_time>::FuncEval(Func F, const InputType a,
                                                                             const InputType b, const InputType *pts)
    : low(InputType(1) / (b - a)), hi(b + a) {
    static_assert(CurrentN > 0, "Polynomial degree must be positive (template N > 0)");
    initialize_monomials(F, pts);
}

template<class Func, std::size_t N_compile_time, std::size_t Iters_compile_time>
template<bool>
constexpr typename FuncEval<Func, N_compile_time, Iters_compile_time>::OutputType PF_ALWAYS_INLINE
FuncEval<Func, N_compile_time, Iters_compile_time>::operator()(const InputType pt) const noexcept {
    const auto xi = map_from_domain(pt);
    return horner<N_compile_time>(xi, monomials.data(), monomials.size()); // Pass data pointer and size
}

// Batch evaluation implementation using SIMD and unrolling
PF_FAST_EVAL_BEGIN
template<class Func, std::size_t N_compile_time, std::size_t Iters_compile_time>
template<int OuterUnrollFactor, bool pts_aligned, bool out_aligned>
PF_ALWAYS_INLINE constexpr void FuncEval<Func, N_compile_time, Iters_compile_time>::horner_polyeval(
    const InputType *PF_RESTRICT pts, OutputType *PF_RESTRICT out, std::size_t num_points) const noexcept {
    return horner<N_compile_time, pts_aligned, out_aligned, OuterUnrollFactor>(
        pts, out, num_points, monomials.data(), monomials.size(),
        [this](const auto v) { return this->map_from_domain(v); });
}
PF_FAST_EVAL_END

// MUST be defined in a c++ source file
// This is a workaround for the compiler to not the inline the function passed to it.
template<typename F, typename... Args> PF_NO_INLINE static auto noinline(F &&f, Args &&...args) {
    return std::forward<F>(f)(std::forward<Args>(args)...);
}

PF_FAST_EVAL_BEGIN
template<class Func, std::size_t N_compile_time, std::size_t Iters_compile_time>
template<bool pts_aligned, bool out_aligned>
PF_ALWAYS_INLINE constexpr void FuncEval<Func, N_compile_time, Iters_compile_time>::operator()(
    const InputType *PF_RESTRICT pts, OutputType *PF_RESTRICT out, std::size_t num_points) const noexcept {
    constexpr auto alignment = xsimd::batch<OutputType>::arch_type::alignment();
#ifdef PF_OUTER_UNROLL
    PF_C23STATIC constexpr auto unroll_factor = PF_OUTER_UNROLL;
#else
    PF_C23STATIC constexpr auto unroll_factor = 0;
#endif

    if constexpr (pts_aligned) {
        if constexpr (out_aligned) {
            return horner_polyeval<unroll_factor, true, true>(pts, out, num_points);
        } else {
            return horner_polyeval<unroll_factor, true, false>(pts, out, num_points);
        }
    } else {
        // if the input and output types are not the same the alignment is tricky
        if constexpr (!std::is_same_v<InputType, OutputType>) {
            return horner_polyeval<unroll_factor, false, false>(pts, out, num_points);
        }
        const auto pts_alignment = detail::get_alignment(pts);
        const auto out_alignment = detail::get_alignment(out);
        if (pts_alignment != out_alignment) PF_UNLIKELY {
                if (pts_alignment >= alignment && out_alignment >= alignment) PF_UNLIKELY {
                        return noinline([this, pts, out, num_points] {
                            return horner_polyeval<unroll_factor, true, false>(pts, out, num_points);
                        });
                    }
                if (out_alignment >= alignment) PF_UNLIKELY {
                        return noinline([this, pts, out, num_points] {
                            return horner_polyeval<unroll_factor, false, true>(pts, out, num_points);
                        });
                    }
                return noinline([this, pts, out, num_points] {
                    return horner_polyeval<unroll_factor, false, false>(pts, out, num_points);
                });
            }

        const auto unaligned_points = std::min(
            ((alignment - pts_alignment) & (alignment - 1)) >> detail::countr_zero(sizeof(InputType)), num_points);

        // Maximum possible unaligned_points = alignment/sizeof(T) - 1 (e.g. 3 for double+AVX2).
        // PF_ASSUME gives the compiler a tight loop-count bound for the scalar prologue.
        constexpr std::size_t scalar_unroll = alignment / sizeof(InputType);
        static_assert(scalar_unroll > 0);
        PF_ASSUME(unaligned_points < scalar_unroll);
        // process scalar until we reach the first aligned point
        for (std::size_t i = 0; i < unaligned_points; ++i) {
            out[i] = operator()(pts[i]);
        }
        return horner_polyeval<unroll_factor, true, true>(pts + unaligned_points, out + unaligned_points,
                                                          num_points - unaligned_points);
    }
}
PF_FAST_EVAL_END

template<class Func, std::size_t N_compile_time, std::size_t Iters_compile_time>
PF_C20CONSTEXPR const Buffer<typename FuncEval<Func, N_compile_time, Iters_compile_time>::OutputType, N_compile_time> &
FuncEval<Func, N_compile_time, Iters_compile_time>::coeffs() const noexcept {
    return monomials;
}

template<class Func, std::size_t N_compile_time, std::size_t Iters_compile_time>
template<class T>
PF_ALWAYS_INLINE constexpr T
FuncEval<Func, N_compile_time, Iters_compile_time>::map_to_domain(const T T_arg) const noexcept {
    return polyfit::internal::helpers::map_to_domain_scalar(T_arg, low, hi);
}

template<class Func, std::size_t N_compile_time, std::size_t Iters_compile_time>
template<class T>
PF_ALWAYS_INLINE constexpr T
FuncEval<Func, N_compile_time, Iters_compile_time>::map_from_domain(const T T_arg) const noexcept {
    return polyfit::internal::helpers::map_from_domain_scalar(T_arg, low, hi);
}

template<class Func, std::size_t N_compile_time, std::size_t Iters_compile_time>
PF_C20CONSTEXPR void FuncEval<Func, N_compile_time, Iters_compile_time>::initialize_monomials(Func F,
                                                                                              const InputType *pts) {
    // 1) allocate
    auto grid = make_buffer<InputType, N_compile_time>(monomials.size());
    auto samples = make_buffer<OutputType, N_compile_time>(monomials.size());

    // 2) fill
    const auto n_terms = monomials.size();
    for (std::size_t k = 0; k < n_terms; ++k) {
        grid[k] =
            pts ? pts[k]
                : InputType(detail::cos((2.0 * double(k) + 1.0) * detail::constants::pi / (2.0 * double(n_terms))));
    }
    for (std::size_t i = 0; i < n_terms; ++i) {
        samples[i] = F(map_to_domain(grid[i]));
    }

    // 3) compute Newton → monomial
    auto newton = detail::bjorck_pereyra<N_compile_time, InputType, OutputType>(grid, samples);
    auto temp_monomial = detail::newton_to_monomial<N_compile_time, InputType, OutputType>(newton, grid);
    assert(temp_monomial.size() == monomials.size() && "size mismatch!");

    std::copy(temp_monomial.begin(), temp_monomial.end(), monomials.begin());

    // 4) optional refine
    refine(grid, samples);

    // 5) fuse domain mapping into coefficients: p(t) → q(x) = p(alpha*x + beta)
    //    where alpha = 2*low, beta = -hi*low. After this, evaluation skips per-point mapping.
    //    Uses compensated arithmetic (twoProd+twoSum) so coefficient precision is O(n^2*eps^2).
    //    Guard: the fused Horner evaluation's condition number is bounded by (|alpha|+|beta|)^deg.
    //    Skip fusion when this exceeds the type's useful precision.
    //
    // std::abs and std::log10 are not constexpr until C++26. Skip domain fusion in
    // constant-evaluated contexts; the polynomial remains correct but uses per-point mapping.
    if (!PF_IS_CONSTANT_EVALUATED()) {
        using Scalar = typename value_type_or_identity<InputType>::type;
        const auto alpha = Scalar(2) * static_cast<Scalar>(low);
        const auto beta = -static_cast<Scalar>(hi) * static_cast<Scalar>(low);
        const auto deg = static_cast<int>(monomials.size()) - 1;
        // Evaluation condition guard: (|alpha|+|beta|)^deg must leave sufficient digits
        const auto cond_base = std::abs(alpha) + std::abs(beta) + Scalar(1);
        constexpr auto max_log = Scalar(std::numeric_limits<Scalar>::digits10 - 3);
        if (deg > 0 && Scalar(deg) * std::log10(cond_base) < max_log) {
            polyfit::internal::helpers::fuse_linear_map(monomials.data(), monomials.size(), alpha, beta);
            low = InputType(0.5);
            hi = InputType(0);
        }
    }
}

template<class Func, std::size_t N_compile_time, std::size_t Iters_compile_time>
PF_C20CONSTEXPR void
FuncEval<Func, N_compile_time, Iters_compile_time>::refine(const Buffer<InputType, N_compile_time> &x_cheb_,
                                                           const Buffer<OutputType, N_compile_time> &y_cheb_) {
    const auto n = monomials.size();
    std::reverse(monomials.begin(), monomials.end()); // to Horner order once

    // Compensated Horner gives more accurate residuals, preventing refinement
    // divergence at high degrees where standard Horner rounding errors exceed
    // the correction magnitude from bjorck_pereyra + newton_to_monomial.
    // Extra passes for n > 32: the compensated residuals converge where standard ones diverge.
    const std::size_t total_iters = Iters_compile_time + (n > 32 ? 2 : 0);

    for (std::size_t pass = 0; pass < total_iters; ++pass) {
        auto r_cheb = make_buffer<OutputType, N_compile_time>(n);
        for (std::size_t i = 0; i < n; ++i) {
            auto p_val = poly_eval::compensated_horner<N_compile_time>(x_cheb_[i], monomials.data(), n);
            r_cheb[i] = y_cheb_[i] - p_val;
        }
        auto newton_r = detail::bjorck_pereyra<N_compile_time, InputType, OutputType>(x_cheb_, r_cheb);
        auto mono_r = detail::newton_to_monomial<N_compile_time, InputType, OutputType>(newton_r, x_cheb_);

        // mono_r is low-degree-first; monomials is high-degree-first — add reversed
        for (std::size_t j = 0; j < n; ++j) monomials[n - 1 - j] += mono_r[j];
    }
    // monomials are now in Horner order (high-degree-first) — done
}

// ------------------------------ truncate -------------------------------

template<class Func, std::size_t N_compile_time, std::size_t Iters_compile_time>
PF_C20CONSTEXPR void FuncEval<Func, N_compile_time, Iters_compile_time>::truncate(
    typename value_type_or_identity<OutputType>::type eps) noexcept {
    if constexpr (N_compile_time == 0) {
        // monomials are in Horner order: [0]=highest degree, [size-1]=constant
        std::size_t skip = 0;
        while (skip + 1 < monomials.size() && std::abs(monomials[skip]) < eps) ++skip;
        if (skip > 0) monomials.erase(monomials.begin(), monomials.begin() + static_cast<std::ptrdiff_t>(skip));
    }
}

// ------------------------------ ctor -----------------------------------

template<typename... EvalTypes> PF_C20CONSTEXPR FuncEvalMany<EvalTypes...>::FuncEvalMany(const EvalTypes &...evals) {
    /* Copy per‑poly scaling */
    auto tmp_low = std::array<InputType, kF>{evals.low...};
    auto tmp_hi = std::array<InputType, kF>{evals.hi...};
    for (std::size_t i = 0; i < kF; ++i) {
        low_[i] = tmp_low[i];
        hi_[i] = tmp_hi[i];
    }

    /* Degree‑dependent storage */
    if constexpr (deg_max_ctime_ == 0) {
        const std::size_t deg_max = std::max({std::size_t(evals.monomials.size())...});
        coeff_store_.assign(kF_pad * deg_max, OutputType{});
        coeffs_ = decltype(coeffs_){coeff_store_.data(), deg_max, kF_pad};
    } else {
        coeffs_ = decltype(coeffs_){coeff_store_.data(), deg_max_ctime_, kF_pad};
    }

    /* Copy coefficients and pad */
    copy_coeffs<0>(evals...);
    zero_pad_coeffs();
}

template<typename... EvalTypes>
PF_C20CONSTEXPR FuncEvalMany<EvalTypes...>::FuncEvalMany(const FuncEvalMany &other)
    : coeff_store_(other.coeff_store_), coeffs_{coeff_store_.data(), other.coeffs_.extents()}, low_(other.low_),
      hi_(other.hi_), truncated_(other.truncated_) {}

template<typename... EvalTypes>
PF_C20CONSTEXPR auto FuncEvalMany<EvalTypes...>::operator=(const FuncEvalMany &other) -> FuncEvalMany & {
    if (this != &other) {
        coeff_store_ = other.coeff_store_;
        low_ = other.low_;
        hi_ = other.hi_;
        truncated_ = other.truncated_;
        coeffs_ = decltype(coeffs_){coeff_store_.data(), other.coeffs_.extents()};
    }
    return *this;
}

template<typename... EvalTypes>
PF_C20CONSTEXPR FuncEvalMany<EvalTypes...>::FuncEvalMany(FuncEvalMany &&other) noexcept
    : coeff_store_(std::move(other.coeff_store_)), coeffs_{coeff_store_.data(), other.coeffs_.extents()},
      low_(std::move(other.low_)), hi_(std::move(other.hi_)), truncated_(other.truncated_) {}

template<typename... EvalTypes>
PF_C20CONSTEXPR auto FuncEvalMany<EvalTypes...>::operator=(FuncEvalMany &&other) noexcept -> FuncEvalMany & {
    if (this != &other) {
        coeff_store_ = std::move(other.coeff_store_);
        low_ = std::move(other.low_);
        hi_ = std::move(other.hi_);
        truncated_ = other.truncated_;
        coeffs_ = decltype(coeffs_){coeff_store_.data(), other.coeffs_.extents()};
    }
    return *this;
}

// ------------------------------ size / degree --------------------------

template<typename... EvalTypes> constexpr std::size_t FuncEvalMany<EvalTypes...>::size() const noexcept { return kF; }

template<typename... EvalTypes> constexpr std::size_t FuncEvalMany<EvalTypes...>::degree() const noexcept {
    return coeffs_.extent(0);
}

// ------------------------------ scalar eval ----------------------------

PF_FAST_EVAL_BEGIN
template<typename... EvalTypes>
auto FuncEvalMany<EvalTypes...>::operator()(InputType x) const noexcept -> std::array<OutputType, kF> {
    alignas(kAlignment) std::array<InputType, kF_pad> xu{};
    for (std::size_t i = 0; i < kF; ++i) xu[i] = xsimd::fms(InputType(2.0), x, hi_[i]) * low_[i];

    alignas(kAlignment) std::array<OutputType, kF_pad> res{};
    if (truncated_)
        horner_transposed<kF_pad, 0, vector_width, true>(xu.data(), coeffs_.data_handle(), res.data(), kF_pad,
                                                         static_cast<std::size_t>(coeffs_.extent(0)));
    else
        horner_transposed<kF_pad, deg_max_ctime_, vector_width, true>(
            xu.data(), coeffs_.data_handle(), res.data(), kF_pad, static_cast<std::size_t>(coeffs_.extent(0)));

    if constexpr (kF == kF_pad) {
        return res; // no padding
    }
    return extract_real(res);
}
PF_FAST_EVAL_END

// ------------------------------ array eval -----------------------------

PF_FAST_EVAL_BEGIN
template<typename... EvalTypes>
auto FuncEvalMany<EvalTypes...>::operator()(const std::array<InputType, kF> &xs) const noexcept
    -> std::array<OutputType, kF> {
    alignas(kAlignment) std::array<InputType, kF_pad> xu{};
    for (std::size_t i = 0; i < kF; ++i) xu[i] = xsimd::fms(InputType(2.0), xs[i], hi_[i]) * low_[i];

    alignas(kAlignment) std::array<OutputType, kF_pad> res{};
    if (truncated_)
        horner_transposed<kF_pad, 0, vector_width, true>(xu.data(), coeffs_.data_handle(), res.data(), kF_pad,
                                                         static_cast<std::size_t>(coeffs_.extent(0)));
    else
        horner_transposed<kF_pad, deg_max_ctime_, vector_width, true>(
            xu.data(), coeffs_.data_handle(), res.data(), kF_pad, static_cast<std::size_t>(coeffs_.extent(0)));
    return extract_real(res);
}
PF_FAST_EVAL_END

// ------------------------------ bulk eval ------------------------------

PF_FAST_EVAL_BEGIN
template<typename... EvalTypes>
void FuncEvalMany<EvalTypes...>::operator()(const InputType *x, OutputType *out,
                                            std::size_t num_points) const noexcept {
    constexpr std::size_t M = kF;
    constexpr std::size_t simd_size = xsimd::batch<InputType>::size;
    const std::size_t n_deg = static_cast<std::size_t>(coeffs_.extent(0));
    const std::size_t stride = kF_pad;

    // For each polynomial, evaluate all points using SIMD across points
    poet::static_for<M>([&](auto m) {
        const auto low_m = low_[m];
        const auto hi_m = hi_[m];

        // Gather this polynomial's coefficients (column m from transposed layout)
        const OutputType *col = coeffs_.data_handle();

        // SIMD loop: process simd_size points at a time
        std::size_t i = 0;
        for (; i + simd_size <= num_points; i += simd_size) {
            // Load simd_size x-values and apply domain mapping
            auto xv = xsimd::load_unaligned(x + i);
            auto xu = xsimd::fms(xsimd::batch<InputType>(InputType(2.0)), xv, xsimd::batch<InputType>(hi_m)) *
                      xsimd::batch<InputType>(low_m);

            // Horner evaluation across points
            auto acc = xsimd::batch<OutputType>(col[0 * stride + m]);
            if constexpr (deg_max_ctime_ != 0) {
                if (!truncated_) {
                    poet::static_for<1, deg_max_ctime_>([&](auto k) {
                        acc =
                            poly_eval::detail::fma(acc, xu, xsimd::batch<OutputType>(col[std::size_t(k) * stride + m]));
                    });
                } else {
                    for (std::size_t k = 1; k < n_deg; ++k)
                        acc = poly_eval::detail::fma(acc, xu, xsimd::batch<OutputType>(col[k * stride + m]));
                }
            } else {
                for (std::size_t k = 1; k < n_deg; ++k)
                    acc = poly_eval::detail::fma(acc, xu, xsimd::batch<OutputType>(col[k * stride + m]));
            }

            // Scatter to row-major output: out[i*M + m], out[(i+1)*M + m], ...
            alignas(xsimd::batch<OutputType>::arch_type::alignment()) OutputType tmp[simd_size];
            acc.store_aligned(tmp);
            for (std::size_t s = 0; s < simd_size; ++s) {
                out[(i + s) * M + m] = tmp[s];
            }
        }

        // Scalar remainder
        for (; i < num_points; ++i) {
            auto xu = (InputType(2.0) * x[i] - hi_m) * low_m;
            OutputType acc = col[0 * stride + m];
            if constexpr (deg_max_ctime_ != 0) {
                if (!truncated_) {
                    poet::static_for<1, deg_max_ctime_>(
                        [&](auto k) { acc = poly_eval::detail::fma(acc, xu, col[std::size_t(k) * stride + m]); });
                } else {
                    for (std::size_t k = 1; k < n_deg; ++k) acc = poly_eval::detail::fma(acc, xu, col[k * stride + m]);
                }
            } else {
                for (std::size_t k = 1; k < n_deg; ++k) acc = poly_eval::detail::fma(acc, xu, col[k * stride + m]);
            }
            out[i * M + m] = acc;
        }
    });
}
PF_FAST_EVAL_END

// ------------------------------ variadic convenience -------------------

template<typename... EvalTypes>
template<typename... Ts>
auto FuncEvalMany<EvalTypes...>::operator()(InputType first, Ts... rest) const noexcept -> std::array<OutputType, kF> {
    static_assert(sizeof...(Ts) + 1 == kF, "Incorrect number of arguments");
    return operator()(std::array<InputType, kF>{first, static_cast<InputType>(rest)...});
}

// ------------------------------ tuple convenience ----------------------

template<typename... EvalTypes>
template<typename... Ts>
auto FuncEvalMany<EvalTypes...>::operator()(const std::tuple<Ts...> &tup) const noexcept -> std::array<OutputType, kF> {
    static_assert(sizeof...(Ts) == kF, "Tuple size must equal number of polynomials");
    std::array<InputType, kF> xs{};
    std::apply([&](auto &&...e) { xs = {static_cast<InputType>(e)...}; }, tup);
    return operator()(xs);
}

// ------------------------------ copy_coeffs ----------------------------

template<typename... EvalTypes>
template<std::size_t I, typename FE, typename... Rest>
PF_C20CONSTEXPR void FuncEvalMany<EvalTypes...>::copy_coeffs(const FE &eval, const Rest &...rest) {
    for (std::size_t k = 0; k < eval.monomials.size(); ++k) coeffs_(k, I) = eval.monomials[k];
    for (std::size_t k = eval.monomials.size(); k < coeffs_.extent(0); ++k) coeffs_(k, I) = OutputType{};
    if constexpr (I + 1 < kF) copy_coeffs<I + 1>(rest...);
}

// ------------------------------ zero_pad_coeffs ------------------------

template<typename... EvalTypes> PF_C20CONSTEXPR void FuncEvalMany<EvalTypes...>::zero_pad_coeffs() {
    for (std::size_t j = kF; j < kF_pad; ++j)
        for (std::size_t k = 0; k < coeffs_.extent(0); ++k) coeffs_(k, j) = OutputType{};
}

// ------------------------------ truncate -------------------------------

template<typename... EvalTypes>
PF_C20CONSTEXPR void FuncEvalMany<EvalTypes...>::truncate(typename value_type_or_identity<OutputType>::type eps) {
    // Scan from highest-degree row downward
    std::size_t new_deg = coeffs_.extent(0);
    while (new_deg > 1) {
        OutputType row_max{};
        for (std::size_t j = 0; j < kF; ++j) row_max = std::max(row_max, std::abs(coeffs_(new_deg - 1, j)));
        if (row_max >= eps) break;
        --new_deg;
    }
    if (new_deg < coeffs_.extent(0)) {
        coeffs_ = decltype(coeffs_){coeff_store_.data(), new_deg, kF_pad};
        truncated_ = true;
    }
}

// ------------------------------ extract_real ---------------------------

template<typename... EvalTypes>
constexpr auto FuncEvalMany<EvalTypes...>::extract_real(const std::array<OutputType, kF_pad> &full) const noexcept
    -> std::array<OutputType, kF> {
    if constexpr (kF == kF_pad) {
        return full; // no padding
    }
    std::array<OutputType, kF> out{};
    for (std::size_t i = 0; i < kF; ++i) out[i] = full[i];
    return out;
}

// Constructor for static degree
template<class Func, std::size_t N_compile>
template<std::size_t C, typename>
constexpr FuncEvalND<Func, N_compile>::FuncEvalND(Func f, const InputType &a, const InputType &b)
    : coeffs_flat_(), coeffs_md_{coeffs_flat_.data(), extents_t{}} {
    compute_scaling(a, b);
    initialize(static_cast<int>(N_compile), f);
}

// Constructor for dynamic degree
template<class Func, std::size_t N_compile>
template<std::size_t C, typename>
constexpr FuncEvalND<Func, N_compile>::FuncEvalND(Func f, int n, const InputType &a, const InputType &b)
    : coeffs_flat_(storage_required(detail::validate_positive_degree(n))),
      coeffs_md_{coeffs_flat_.data(), make_ext(detail::validate_positive_degree(n))} {
    const auto checked_n = detail::validate_positive_degree(n);
    compute_scaling(a, b);
    initialize(checked_n, f);
}

template<class Func, std::size_t N_compile>
FuncEvalND<Func, N_compile>::FuncEvalND(const FuncEvalND &other)
    : low_(other.low_), hi_(other.hi_), coeffs_flat_(other.coeffs_flat_),
      coeffs_md_{coeffs_flat_.data(), other.coeffs_md_.extents()} {}

template<class Func, std::size_t N_compile>
auto FuncEvalND<Func, N_compile>::operator=(const FuncEvalND &other) -> FuncEvalND & {
    if (this != &other) {
        low_ = other.low_;
        hi_ = other.hi_;
        coeffs_flat_ = other.coeffs_flat_;
        coeffs_md_ = mdspan_t{coeffs_flat_.data(), other.coeffs_md_.extents()};
    }
    return *this;
}

template<class Func, std::size_t N_compile>
FuncEvalND<Func, N_compile>::FuncEvalND(FuncEvalND &&other) noexcept
    : low_(std::move(other.low_)), hi_(std::move(other.hi_)), coeffs_flat_(std::move(other.coeffs_flat_)),
      coeffs_md_{coeffs_flat_.data(), other.coeffs_md_.extents()} {}

template<class Func, std::size_t N_compile>
auto FuncEvalND<Func, N_compile>::operator=(FuncEvalND &&other) noexcept -> FuncEvalND & {
    if (this != &other) {
        low_ = std::move(other.low_);
        hi_ = std::move(other.hi_);
        coeffs_flat_ = std::move(other.coeffs_flat_);
        coeffs_md_ = mdspan_t{coeffs_flat_.data(), other.coeffs_md_.extents()};
    }
    return *this;
}

// Evaluate via Horner's method
PF_FAST_EVAL_BEGIN
template<class Func, std::size_t N_compile>
template<bool SIMD>
constexpr typename FuncEvalND<Func, N_compile>::OutputType
FuncEvalND<Func, N_compile>::operator()(const InputType &x) const {
    const int deg_rt = (N_compile ? static_cast<int>(N_compile) : static_cast<int>(coeffs_md_.extent(0)));
    return poly_eval::horner<N_compile, SIMD, OutputType>(map_from_domain(x), coeffs_md_, deg_rt);
}
PF_FAST_EVAL_END

// coeff_impl
template<class Func, std::size_t N_compile>
template<typename IdxArray, std::size_t... I>
constexpr typename FuncEvalND<Func, N_compile>::Scalar &
FuncEvalND<Func, N_compile>::coeff_impl(const IdxArray &idx, std::size_t k, std::index_sequence<I...>) noexcept {
    return coeffs_md_(static_cast<std::size_t>(idx[I])..., k);
}

template<class Func, std::size_t N_compile>
template<class IdxArray>
[[nodiscard]] constexpr typename FuncEvalND<Func, N_compile>::Scalar &
FuncEvalND<Func, N_compile>::coeff(const IdxArray &idx, std::size_t k) noexcept {
    return coeff_impl<IdxArray>(idx, k, std::make_index_sequence<dim_>{});
}

template<class Func, std::size_t N_compile> auto FuncEvalND<Func, N_compile>::make_ext(int n) noexcept -> extents_t {
    if constexpr (is_static) {
        return detail::make_static_extents<N_compile, dim_, outDim_>(std::make_index_sequence<dim_>{});
    } else {
        return make_ext(n, std::make_index_sequence<dim_ + 1>{});
    }
}

template<class Func, std::size_t N_compile>
template<std::size_t... Is>
auto FuncEvalND<Func, N_compile>::make_ext(int n, std::index_sequence<Is...>) noexcept -> extents_t {
    return extents_t{(Is < dim_ ? static_cast<std::size_t>(n) : static_cast<std::size_t>(outDim_))...};
}

template<class Func, std::size_t N_compile>
constexpr std::size_t FuncEvalND<Func, N_compile>::storage_required(const int n) noexcept {
    auto ext = make_ext(n);
    auto mapping = typename mdspan_t::mapping_type{ext};
    return mapping.required_span_size();
}
template<class Func, std::size_t N_compile> constexpr void FuncEvalND<Func, N_compile>::initialize(int n, Func f) {
    const auto n_nodes = static_cast<std::size_t>(n);
    auto nodes = make_buffer<Scalar, N_compile>(n_nodes);
    for (std::size_t k = 0; k < n_nodes; ++k)
        nodes[k] = detail::cos((2.0 * double(k) + 1.0) * detail::constants::pi / (2.0 * double(n)));

    std::array<int, dim_> ext_idx{};
    ext_idx.fill(n);

    // sample f on Chebyshev grid
    for_each_index<dim_>(ext_idx, [&](const std::array<int, dim_> &idx) {
        InputType x_dom{};
        for (std::size_t d = 0; d < dim_; ++d) x_dom[d] = nodes[static_cast<std::size_t>(idx[d])];
        OutputType y = f(map_to_domain(x_dom));
        for (std::size_t k = 0; k < outDim_; ++k) coeff(idx, k) = y[k];
    });

    convert_newton_to_monomial_axes(n, nodes);
    reverse_coefficients_axes(n);
}

template<class Func, std::size_t N_compile>
constexpr void FuncEvalND<Func, N_compile>::convert_newton_to_monomial_axes(int n,
                                                                            const Buffer<Scalar, N_compile> &nodes) {
    const auto n_nodes = static_cast<std::size_t>(n);
    auto rhs = make_buffer<Scalar, N_compile>(n_nodes);
    auto alpha = make_buffer<Scalar, N_compile>(n_nodes);
    auto mono = make_buffer<Scalar, N_compile>(n_nodes);

    std::array<int, dim_> ext_idx{};
    ext_idx.fill(n);
    std::array<int, dim_> base_idx{};
    for (std::size_t axis = 0; axis < dim_; ++axis) {
        auto inner_ext = ext_idx;
        inner_ext[axis] = 1;
        for_each_index<dim_>(inner_ext, [&](const std::array<int, dim_> &base) {
            for (std::size_t k = 0; k < outDim_; ++k) {
                for (std::size_t i = 0; i < n_nodes; ++i) {
                    base_idx = base;
                    base_idx[axis] = static_cast<int>(i);
                    rhs[i] = coeff(base_idx, k);
                }
                alpha = detail::bjorck_pereyra<N_compile, Scalar, Scalar>(nodes, rhs);
                mono = detail::newton_to_monomial<N_compile, Scalar, Scalar>(alpha, nodes);
                for (std::size_t i = 0; i < n_nodes; ++i) {
                    base_idx = base;
                    base_idx[axis] = static_cast<int>(i);
                    coeff(base_idx, k) = mono[i];
                }
            }
        });
    }
}

template<class Func, std::size_t N_compile>
constexpr void FuncEvalND<Func, N_compile>::reverse_coefficients_axes(int n) {
    std::array<int, dim_> ext_idx{};
    ext_idx.fill(n);
    std::array<int, dim_> base_idx{};
    for (std::size_t axis = 0; axis < dim_; ++axis) {
        auto inner_ext = ext_idx;
        inner_ext[axis] = 1;
        for_each_index<dim_>(inner_ext, [&](const std::array<int, dim_> &base) {
            for (std::size_t k = 0; k < outDim_; ++k) {
                int i = 0, j = n - 1;
                while (i < j) {
                    base_idx = base;
                    base_idx[axis] = i;
                    auto &a = coeff(base_idx, k);
                    base_idx[axis] = j;
                    auto &b = coeff(base_idx, k);
                    std::swap(a, b);
                    ++i;
                    --j;
                }
            }
        });
    }
}

template<class Func, std::size_t N_compile>
[[nodiscard]] constexpr typename FuncEvalND<Func, N_compile>::InputType
FuncEvalND<Func, N_compile>::map_to_domain(const InputType &t) const noexcept {
    return polyfit::internal::helpers::map_to_domain_array<Scalar, dim_>(t, low_, hi_);
}

template<class Func, std::size_t N_compile>
[[nodiscard]] constexpr typename FuncEvalND<Func, N_compile>::InputType
FuncEvalND<Func, N_compile>::map_from_domain(const InputType &x) const noexcept {
    return polyfit::internal::helpers::map_from_domain_array<Scalar, dim_>(x, low_, hi_);
}

template<class Func, std::size_t N_compile>
constexpr void FuncEvalND<Func, N_compile>::compute_scaling(const InputType &a, const InputType &b) noexcept {
    polyfit::internal::helpers::compute_scaling_array<Scalar, dim_>(a, b, low_, hi_);
}

template<class Func, std::size_t N_compile>
// Odometer-style multi-index iteration: visits every index tuple in the
// Cartesian product [0, ext[0]) × [0, ext[1]) × … × [0, ext[Rank-1]).
// The innermost (d=0) dimension increments fastest, like a little-endian
// counter. Each visited index tuple is passed to `body`.
template<std::size_t Rank, class F>
constexpr void FuncEvalND<Func, N_compile>::for_each_index(const std::array<int, Rank> &ext, F &&body) {
    std::array<int, Rank> idx{};
    while (true) {
        body(idx);
        for (std::size_t d = 0; d < Rank; ++d) {
            if (++idx[d] < ext[d]) break;
            if (d == Rank - 1) return;
            idx[d] = 0;
        }
    }
}

// -----------------------------------------------------------------------------
// make_func_eval API implementations (Runtime, C++20 compatible)
// -----------------------------------------------------------------------------
// Compile-time degree (1D or ND)
template<std::size_t N_compile_time, std::size_t Iters_compile_time, class Func>
[[nodiscard]] PF_C20CONSTEXPR auto make_func_eval(Func F, typename function_traits<Func>::arg0_type a,
                                                  typename function_traits<Func>::arg0_type b,
                                                  const typename function_traits<Func>::arg0_type *pts) {
    using InputType = typename function_traits<Func>::arg0_type;
    if constexpr (has_tuple_size_v<std::remove_cvref_t<InputType>>) {
        // ND
        return FuncEvalND<Func, N_compile_time>(F, a, b);
    } else {
        // 1D
        return FuncEval<Func, N_compile_time, Iters_compile_time>(F, a, b, pts);
    }
}

// Runtime degree (1D or ND)
template<std::size_t Iters_compile_time, class Func, typename IntType,
         std::enable_if_t<std::is_integral_v<std::remove_cvref_t<IntType>>, int>>
[[nodiscard]] PF_C20CONSTEXPR auto
make_func_eval(Func F, IntType n, typename function_traits<Func>::arg0_type a,
               typename function_traits<Func>::arg0_type b,
               const std::remove_reference_t<typename function_traits<Func>::arg0_type> *pts) {
    using InputType = typename function_traits<Func>::arg0_type;
    using RawInputType = std::remove_cvref_t<InputType>;

    if constexpr (has_tuple_size_v<RawInputType>) {
        // ND - pass by value for tuple types to avoid reference issues
        return FuncEvalND<Func, 0>(F, detail::validate_positive_degree(static_cast<int>(n)), a, b);
    } else {
        // 1D
        return FuncEval<Func, 0, Iters_compile_time>(F, detail::validate_positive_degree(static_cast<int>(n)), a, b,
                                                     pts);
    }
}

template<std::size_t N_compile_time, class Func,
         std::enable_if_t<has_tuple_size_v<std::remove_cvref_t<typename function_traits<Func>::arg0_type>>, int>>
[[nodiscard]] PF_C20CONSTEXPR auto make_func_eval(Func F, typename function_traits<Func>::arg0_type a,
                                                  typename function_traits<Func>::arg0_type b) {
    return FuncEvalND<Func, N_compile_time>(F, a, b);
}

// Runtime error tolerance (1D or ND) — fit once at MaxN, then truncate
template<std::size_t MaxN_val, std::size_t NumEvalPoints_val, std::size_t Iters_compile_time, class Func,
         typename FloatType, std::enable_if_t<std::is_floating_point_v<std::remove_cvref_t<FloatType>>, int>>
[[nodiscard]] PF_C20CONSTEXPR auto make_func_eval(Func F, FloatType eps, typename function_traits<Func>::arg0_type a,
                                                  typename function_traits<Func>::arg0_type b) {
    using RawInputType = std::remove_cvref_t<typename function_traits<Func>::arg0_type>;
    using evaluator_t =
        std::conditional_t<has_tuple_size_v<RawInputType>, FuncEvalND<Func, 0>, FuncEval<Func, 0, Iters_compile_time>>;
    const auto eval_pts = detail::linspace(a, b, int(NumEvalPoints_val));

    // 1. Fit once at MaxN
    evaluator_t evaluator(F, int(MaxN_val), a, b);

    // Helper to compute max error across eval points
    auto compute_max_err = [&]() {
        double max_err = 0.0;
        for (const auto &pt : eval_pts) {
            max_err = std::max(detail::relative_l2_norm(F(pt), evaluator(pt)), max_err);
        }
        return max_err;
    };

    // 2. Try truncating negligible high-degree coefficients
    if constexpr (!has_tuple_size_v<RawInputType>) {
        auto candidate = evaluator;
        candidate.truncate(static_cast<typename value_type_or_identity<typename evaluator_t::OutputType>::type>(eps));
        double err = 0.0;
        for (const auto &pt : eval_pts) err = std::max(detail::relative_l2_norm(F(pt), candidate(pt)), err);
        if (err <= eps) return candidate;
    }

    // 3. Full-degree fit — verify error
    if (compute_max_err() <= eps) return evaluator;

    throw std::runtime_error("No polynomial degree found for requested error tolerance. eps=" + std::to_string(eps) +
                             ", MaxN=" + std::to_string(MaxN_val) + ", max_err=" + std::to_string(compute_max_err()));
}

template<std::size_t N_compile_time, std::size_t Iters_compile_time, typename Func,
         std::enable_if_t<std::is_function_v<std::remove_pointer_t<std::decay_t<Func>>>, int>>
[[nodiscard]] PF_C20CONSTEXPR auto make_func_eval(Func *f, typename function_traits<Func *>::arg0_type a,
                                                  typename function_traits<Func *>::arg0_type b) {
    using InputType = typename function_traits<Func *>::arg0_type;
    if constexpr (has_tuple_size_v<std::remove_cvref_t<InputType>>) {
        // ND
        auto func_wrapper = [f](const InputType &in) {
            return f(in);
        };
        return FuncEvalND<decltype(func_wrapper), N_compile_time>(func_wrapper, a, b);
    } else {
        // 1D
        return FuncEval<Func *, N_compile_time, Iters_compile_time>(f, a, b, nullptr);
    }
}

#if PF_HAS_CONSTEXPR_EPS_OVERLOAD
template<double eps_val, auto a, auto b, std::size_t MaxN_val, std::size_t NumEvalPoints_val,
         std::size_t Iters_compile_time, class Func>
[[nodiscard]] constexpr auto make_func_eval(Func F) {
    using RawInputType = std::remove_cvref_t<typename function_traits<Func>::arg0_type>;
    static_assert(MaxN_val > 0, "Max polynomial degree must be positive.");
    static_assert(NumEvalPoints_val > 1, "Number of evaluation points must be greater than 1.");

    // Recursive constexpr search for minimal degree
    constexpr auto degree = [F] {
        constexpr auto compute_error = [F](const auto &evaluator) {
            constexpr auto eval_pts = detail::linspace<static_cast<int>(NumEvalPoints_val)>(a, b);
            double max_err = 0.0;
            for (const auto &pt : eval_pts) {
                const auto actual = F(pt);
                const auto approx = evaluator.template operator()<false>(pt);
                max_err = std::max(detail::relative_l2_norm(actual, approx), max_err);
            }
            return max_err;
        };
        int n = 0;
        poet::static_for<2, MaxN_val>([&](auto i) {
            if (n != 0) return; // already found minimal degree
            using evaluator_t = std::conditional_t<has_tuple_size_v<RawInputType>, FuncEvalND<Func, i>,
                                                   FuncEval<Func, i, Iters_compile_time>>;
            if constexpr (compute_error(evaluator_t(F, a, b)) <= eps_val) {
                n = i;
            }
        });
        return n;
    }();
    static_assert(degree != 0, "No polynomial degree found for requested error tolerance.");
    using evaluator_t = std::conditional_t<has_tuple_size_v<RawInputType>, FuncEvalND<Func, degree>,
                                           FuncEval<Func, degree, Iters_compile_time>>;
    return evaluator_t(F, a, b);
}
#endif // PF_HAS_CONSTEXPR_EPS_OVERLOAD

template<typename... EvalTypes>
[[nodiscard]] PF_C20CONSTEXPR FuncEvalMany<EvalTypes...> make_func_eval_many(EvalTypes... evals) noexcept {
    return FuncEvalMany<std::decay_t<EvalTypes>...>(std::forward<EvalTypes>(evals)...);
}

} // namespace poly_eval
