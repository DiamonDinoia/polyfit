#pragma once

#include <cassert>
#include <iomanip>
#include <iostream>
#include <xsimd/xsimd.hpp>

#include "macros.h"
#include "poly_eval.h"
#include "utils.h"

namespace poly_eval {

// -----------------------------------------------------------------------------
// FuncEval Implementation (Runtime)
// -----------------------------------------------------------------------------

template <class Func, std::size_t N_compile_time, std::size_t Iters_compile_time>
template <std::size_t CurrentN, typename>
PF_C20CONSTEXPR FuncEval<Func, N_compile_time, Iters_compile_time>::FuncEval(Func F, const int n, const InputType a,
                                                                            const InputType b, const InputType *pts)
    : n_terms(n), low(InputType(1) / (b - a)), hi(b + a) {
    // if constexpr (N_compile_time == 0) {
    //     assert(n_terms > 0 && "Polynomial degree must be positive");
    // }
    monomials.resize(n_terms);
    initialize_monomials(F, pts);
}

template <class Func, std::size_t N_compile_time, std::size_t Iters_compile_time>
template <std::size_t CurrentN, typename>
PF_C20CONSTEXPR FuncEval<Func, N_compile_time, Iters_compile_time>::FuncEval(Func F, const InputType a, const InputType b,
                                                                            const InputType *pts)
    : n_terms(static_cast<int>(CurrentN)), low(InputType(1) / (b - a)), hi(b + a) {
    assert(n_terms > 0 && "Polynomial degree must be positive (template N > 0)");
    initialize_monomials(F, pts);
}

template <class Func, std::size_t N_compile_time, std::size_t Iters_compile_time>
template <bool>
typename FuncEval<Func, N_compile_time, Iters_compile_time>::OutputType constexpr PF_ALWAYS_INLINE
FuncEval<Func, N_compile_time, Iters_compile_time>::operator()(const InputType pt) const noexcept {
    const auto xi = map_from_domain(pt);
    return horner<N_compile_time>(xi, monomials.data(), monomials.size()); // Pass data pointer and size
}

// Batch evaluation implementation using SIMD and unrolling
template <class Func, std::size_t N_compile_time, std::size_t Iters_compile_time>
template <int OuterUnrollFactor, bool pts_aligned, bool out_aligned>
PF_ALWAYS_INLINE constexpr void FuncEval<Func, N_compile_time, Iters_compile_time>::horner_polyeval(
    const InputType *PF_RESTRICT pts, OutputType *PF_RESTRICT out, std::size_t num_points) const noexcept {
    static_assert(OuterUnrollFactor >= 0 && (OuterUnrollFactor & (OuterUnrollFactor - 1)) == 0,
                  "OuterUnrollFactor must be a power of two greater than zero.");
    return horner<N_compile_time, pts_aligned, out_aligned, OuterUnrollFactor>(
        pts, out, num_points, monomials.data(), monomials.size(), [this](const auto v) { return map_from_domain(v); });
}

// MUST be defined in a c++ source file
// This is a workaround for the compiler to not the inline the function passed to it.
template <typename F, typename... Args> PF_NO_INLINE static auto noinline(F &&f, Args &&...args) {
    return std::forward<F>(f)(std::forward<Args>(args)...);
}

template <class Func, std::size_t N_compile_time, std::size_t Iters_compile_time>
template <bool pts_aligned, bool out_aligned>
PF_ALWAYS_INLINE void constexpr FuncEval<Func, N_compile_time, Iters_compile_time>::operator()(
    const InputType * PF_RESTRICT pts, OutputType * PF_RESTRICT out, std::size_t num_points) const noexcept {
    // find out the alignment of pts and out
    // constexpr auto simd_size = xsimd::batch<InputType>::size;
    constexpr auto alignment = xsimd::batch<OutputType>::arch_type::alignment();
    constexpr auto unroll_factor = 0;

    // const auto monomial_ptr = monomials.data();
    // const auto monomial_size = monomials.size();

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
        if (pts_alignment != out_alignment) [[unlikely]] {
            if (pts_alignment >= alignment && out_alignment >= alignment) [[unlikely]] {
                return noinline([this, pts, out, num_points] {
                    return horner_polyeval<unroll_factor, true, false>(pts, out, num_points);
                });
            }
            if (out_alignment >= alignment) [[unlikely]] {
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

        constexpr std::size_t min_align = alignof(std::max_align_t); // in bytes, typically 16
        constexpr std::size_t scalar_unroll = (alignment - min_align) / sizeof(InputType);

        // print alignment;
        PF_ASSUME(unaligned_points < scalar_unroll); // tells the compiler that this loop is at most alignment
        // process scalar until we reach the first aligned point
        for (std::size_t i = 0; i < unaligned_points; ++i) {
            out[i] = operator()(pts[i]);
        }
        return horner_polyeval<unroll_factor, true, true>(pts + unaligned_points, out + unaligned_points,
                                                          num_points - unaligned_points);
    }
}

template <class Func, std::size_t N_compile_time, std::size_t Iters_compile_time>
PF_C20CONSTEXPR const Buffer<typename FuncEval<Func, N_compile_time, Iters_compile_time>::OutputType, N_compile_time> &
FuncEval<Func, N_compile_time, Iters_compile_time>::coeffs() const noexcept {
    return monomials;
}

template <class Func, std::size_t N_compile_time, std::size_t Iters_compile_time>
template <class T>
PF_ALWAYS_INLINE constexpr T
FuncEval<Func, N_compile_time, Iters_compile_time>::map_to_domain(const T T_arg) const noexcept {
    return static_cast<T>(0.5 * (T_arg / low + hi));
}

template <class Func, std::size_t N_compile_time, std::size_t Iters_compile_time>
template <class T>
PF_ALWAYS_INLINE constexpr T
FuncEval<Func, N_compile_time, Iters_compile_time>::map_from_domain(const T T_arg) const noexcept {
    if constexpr (std::is_arithmetic_v<T>) {
        return static_cast<T>(std::fma(2.0, T_arg, -T(hi)) * low);
    }
    return static_cast<T>(xsimd::fms(T(2.0), T_arg, T(hi)) * low);
}

template <class Func, std::size_t N_compile_time, std::size_t Iters_compile_time>
PF_C20CONSTEXPR void FuncEval<Func, N_compile_time, Iters_compile_time>::initialize_monomials(Func F,
                                                                                              const InputType *pts) {
    // 1) allocate
    Buffer<InputType, N_compile_time> grid{};
    if constexpr (N_compile_time == 0)
        grid.resize(n_terms);

    Buffer<OutputType, N_compile_time> samples{};
    if constexpr (N_compile_time == 0)
        samples.resize(n_terms);

    // 2) fill
    for (auto k = 0; k < n_terms; ++k) {
        grid[k] = pts ? pts[k] : InputType(detail::cos((2.0 * InputType(k) + 1.0) * M_PI / (2.0 * n_terms)));
    }
    for (auto i = 0; i < n_terms; ++i) {
        samples[i] = F(map_to_domain(grid[i]));
    }

    // 3) compute Newton → monomial
    auto newton = detail::bjorck_pereyra<N_compile_time, InputType, OutputType>(grid, samples);
    auto temp_monomial = detail::newton_to_monomial<N_compile_time, InputType, OutputType>(newton, grid);
    assert(temp_monomial.size() == monomials.size() && "size mismatch!");

    std::copy(temp_monomial.begin(), temp_monomial.end(), monomials.begin());

    // 4) optional refine
    refine(grid, samples);
}

template <class Func, std::size_t N_compile_time, std::size_t Iters_compile_time>
PF_C20CONSTEXPR void
FuncEval<Func, N_compile_time, Iters_compile_time>::refine(const Buffer<InputType, N_compile_time> &x_cheb_,
                                                           const Buffer<OutputType, N_compile_time> &y_cheb_) {

    for (std::size_t pass = 0; pass < Iters_compile_time; ++pass) {
        // residuals
        Buffer<OutputType, N_compile_time> r_cheb;
        if constexpr (N_compile_time == 0) {
            r_cheb.resize(n_terms);
        }
        std::reverse(monomials.begin(), monomials.end());
        for (int i = 0; i < n_terms; ++i) {
            auto xi = x_cheb_[i];
            auto p_val = poly_eval::horner<N_compile_time>(xi, monomials.data(), monomials.size());
            r_cheb[i] = y_cheb_[i] - p_val;
        }
        std::reverse(monomials.begin(), monomials.end());
        // correction
        auto newton_r = detail::bjorck_pereyra<N_compile_time, InputType, OutputType>(x_cheb_, r_cheb);
        auto mono_r = detail::newton_to_monomial<N_compile_time, InputType, OutputType>(newton_r, x_cheb_);
        assert(mono_r.size() == monomials.size());

        for (std::size_t j = 0; j < monomials.size(); ++j) {
            monomials[j] += mono_r[j];
        }
    }
    std::reverse(monomials.begin(), monomials.end());
}

// ------------------------------ ctor -----------------------------------

template <typename... EvalTypes> FuncEvalMany<EvalTypes...>::FuncEvalMany(const EvalTypes &...evals) {
    /* Copy per‑poly scaling */
    auto tmp_low = std::array<InputType, kF>{evals.low...};
    auto tmp_hi = std::array<InputType, kF>{evals.hi...};
    for (std::size_t i = 0; i < kF; ++i) {
        low_[i] = tmp_low[i];
        hi_[i] = tmp_hi[i];
    }

    /* Degree‑dependent storage */
    if constexpr (deg_max_ctime_ == 0) {
        deg_max_ = std::max({evals.n_terms...});
        coeff_store_.assign(kF_pad * deg_max_, OutputType{});
        coeffs_ = decltype(coeffs_){coeff_store_.data(), deg_max_, kF_pad};
    } else {
        coeffs_ = decltype(coeffs_){coeff_store_.data(), deg_max_ctime_, kF_pad};
    }

    /* Copy coefficients and pad */
    copy_coeffs<0>(evals...);
    zero_pad_coeffs();
}

template <typename... EvalTypes>
FuncEvalMany<EvalTypes...>::FuncEvalMany(const FuncEvalMany &other)
    : coeff_store_(other.coeff_store_),
      coeffs_{coeff_store_.data(), other.coeffs_.extents()},
      low_(other.low_), hi_(other.hi_), deg_max_(other.deg_max_) {}

template <typename... EvalTypes>
auto FuncEvalMany<EvalTypes...>::operator=(const FuncEvalMany &other) -> FuncEvalMany & {
    if (this != &other) {
        coeff_store_ = other.coeff_store_;
        low_ = other.low_;
        hi_ = other.hi_;
        deg_max_ = other.deg_max_;
        coeffs_ = decltype(coeffs_){coeff_store_.data(), other.coeffs_.extents()};
    }
    return *this;
}

template <typename... EvalTypes>
FuncEvalMany<EvalTypes...>::FuncEvalMany(FuncEvalMany &&other) noexcept
    : coeff_store_(std::move(other.coeff_store_)),
      coeffs_{coeff_store_.data(), other.coeffs_.extents()},
      low_(std::move(other.low_)), hi_(std::move(other.hi_)), deg_max_(other.deg_max_) {}

template <typename... EvalTypes>
auto FuncEvalMany<EvalTypes...>::operator=(FuncEvalMany &&other) noexcept -> FuncEvalMany & {
    if (this != &other) {
        coeff_store_ = std::move(other.coeff_store_);
        low_ = std::move(other.low_);
        hi_ = std::move(other.hi_);
        deg_max_ = other.deg_max_;
        coeffs_ = decltype(coeffs_){coeff_store_.data(), other.coeffs_.extents()};
    }
    return *this;
}

// ------------------------------ size / degree --------------------------

template <typename... EvalTypes> std::size_t FuncEvalMany<EvalTypes...>::size() const noexcept { return kF; }

template <typename... EvalTypes> std::size_t FuncEvalMany<EvalTypes...>::degree() const noexcept { return deg_max_; }

// ------------------------------ scalar eval ----------------------------

template <typename... EvalTypes>
std::array<typename FuncEvalMany<EvalTypes...>::OutputType, FuncEvalMany<EvalTypes...>::kF>
FuncEvalMany<EvalTypes...>::operator()(InputType x) const noexcept {
    std::array<InputType, kF_pad> xu{};
    for (std::size_t i = 0; i < kF; ++i)
        xu[i] = xsimd::fms(InputType(2.0), x, hi_[i]) * low_[i];

    std::array<OutputType, kF_pad> res{};
    horner_transposed<kF_pad, deg_max_ctime_, vector_width>(xu.data(), coeffs_.data_handle(), res.data(), kF_pad,
                                                            deg_max_);

    if constexpr (kF == kF_pad) {
        return res; // no padding
    }
    return extract_real(res);
}

// ------------------------------ array eval -----------------------------

template <typename... EvalTypes>
std::array<typename FuncEvalMany<EvalTypes...>::OutputType, FuncEvalMany<EvalTypes...>::kF>
FuncEvalMany<EvalTypes...>::operator()(const std::array<InputType, kF> &xs) const noexcept {
    std::array<InputType, kF_pad> xu{};
    for (std::size_t i = 0; i < kF; ++i)
        xu[i] = xsimd::fms(InputType(2.0), xs[i], hi_[i]) * low_[i];

    std::array<OutputType, kF_pad> res{};
    horner_transposed<kF_pad, deg_max_ctime_, vector_width>(xu.data(), coeffs_.data_handle(), res.data(), kF_pad,
                                                            deg_max_);
    return extract_real(res);
}

// ------------------------------ bulk eval ------------------------------

template <typename... EvalTypes>
void FuncEvalMany<EvalTypes...>::operator()(const InputType *x, OutputType *out,
                                            std::size_t num_points) const noexcept {
    constexpr std::size_t M = kF;
    using extents_t =
        stdex::mdspan<OutputType, stdex::extents<std::size_t, stdex::dynamic_extent, M>, stdex::layout_right>;
    extents_t out_m{out, num_points};

    for (std::size_t i = 0; i < num_points; ++i) {
        auto vals = operator()(x[i]);
        detail::unroll_loop<M>([&](const auto I) { out_m(i, decltype(I)::value) = vals[decltype(I)::value]; });
    }
}

// ------------------------------ variadic convenience -------------------

template <typename... EvalTypes>
template <typename... Ts>
std::array<typename FuncEvalMany<EvalTypes...>::OutputType, FuncEvalMany<EvalTypes...>::kF>
FuncEvalMany<EvalTypes...>::operator()(InputType first, Ts... rest) const noexcept {
    static_assert(sizeof...(Ts) + 1 == kF, "Incorrect number of arguments");
    return operator()(std::array<InputType, kF>{first, static_cast<InputType>(rest)...});
}

// ------------------------------ tuple convenience ----------------------

template <typename... EvalTypes>
template <typename... Ts>
std::array<typename FuncEvalMany<EvalTypes...>::OutputType, FuncEvalMany<EvalTypes...>::kF>
FuncEvalMany<EvalTypes...>::operator()(const std::tuple<Ts...> &tup) const noexcept {
    static_assert(sizeof...(Ts) == kF, "Tuple size must equal number of polynomials");
    std::array<InputType, kF> xs{};
    std::apply([&](auto &&...e) { xs = {static_cast<InputType>(e)...}; }, tup);
    return operator()(xs);
}

// ------------------------------ copy_coeffs ----------------------------

template <typename... EvalTypes>
template <std::size_t I, typename FE, typename... Rest>
void FuncEvalMany<EvalTypes...>::copy_coeffs(const FE &fe, const Rest &...rest) {
    for (auto k = 0; k < fe.n_terms; ++k)
        coeffs_(k, I) = fe.monomials[k];
    for (auto k = std::size_t(fe.n_terms); k < deg_max_; ++k)
        coeffs_(k, I) = OutputType{};
    if constexpr (I + 1 < kF)
        copy_coeffs<I + 1>(rest...);
}

// ------------------------------ zero_pad_coeffs ------------------------

template <typename... EvalTypes> void FuncEvalMany<EvalTypes...>::zero_pad_coeffs() {
    for (std::size_t j = kF; j < kF_pad; ++j)
        for (std::size_t k = 0; k < deg_max_; ++k)
            coeffs_(k, j) = OutputType{};
}

// ------------------------------ extract_real ---------------------------

template <typename... EvalTypes>
std::array<typename FuncEvalMany<EvalTypes...>::OutputType, FuncEvalMany<EvalTypes...>::kF>
FuncEvalMany<EvalTypes...>::extract_real(const std::array<OutputType, kF_pad> &full) const noexcept {
    if constexpr (kF == kF_pad) {
        return full; // no padding
    }
    std::array<OutputType, kF> out{};
    for (std::size_t i = 0; i < kF; ++i)
        out[i] = full[i];
    return out;
}

// Constructor for static degree
template <class Func, std::size_t N_compile>
template <std::size_t C, typename>
constexpr FuncEvalND<Func, N_compile>::FuncEvalND(Func f, const InputType &a, const InputType &b)
    : func_{f}, degree_{static_cast<int>(N_compile)}, coeffs_flat_(), coeffs_md_{coeffs_flat_.data(), extents_t{}} {
    compute_scaling(a, b);
    initialize(static_cast<int>(N_compile));
}

// Constructor for dynamic degree
template <class Func, std::size_t N_compile>
template <std::size_t C, typename>
constexpr FuncEvalND<Func, N_compile>::FuncEvalND(Func f, int n, const InputType &a, const InputType &b)
    : func_{f}, degree_{n}, coeffs_flat_(storage_required(n)), coeffs_md_{coeffs_flat_.data(), make_ext(n)} {
    compute_scaling(a, b);
    initialize(n);
}

template <class Func, std::size_t N_compile>
FuncEvalND<Func, N_compile>::FuncEvalND(const FuncEvalND &other)
    : func_(other.func_), degree_(other.degree_), low_(other.low_), hi_(other.hi_),
      coeffs_flat_(other.coeffs_flat_),
      coeffs_md_{coeffs_flat_.data(), other.coeffs_md_.extents()} {}

template <class Func, std::size_t N_compile>
auto FuncEvalND<Func, N_compile>::operator=(const FuncEvalND &other) -> FuncEvalND & {
    if (this != &other) {
        func_ = other.func_;
        degree_ = other.degree_;
        low_ = other.low_;
        hi_ = other.hi_;
        coeffs_flat_ = other.coeffs_flat_;
        coeffs_md_ = mdspan_t{coeffs_flat_.data(), other.coeffs_md_.extents()};
    }
    return *this;
}

template <class Func, std::size_t N_compile>
FuncEvalND<Func, N_compile>::FuncEvalND(FuncEvalND &&other) noexcept
    : func_(std::move(other.func_)), degree_(other.degree_), low_(std::move(other.low_)), hi_(std::move(other.hi_)),
      coeffs_flat_(std::move(other.coeffs_flat_)),
      coeffs_md_{coeffs_flat_.data(), other.coeffs_md_.extents()} {}

template <class Func, std::size_t N_compile>
auto FuncEvalND<Func, N_compile>::operator=(FuncEvalND &&other) noexcept -> FuncEvalND & {
    if (this != &other) {
        func_ = std::move(other.func_);
        degree_ = other.degree_;
        low_ = std::move(other.low_);
        hi_ = std::move(other.hi_);
        coeffs_flat_ = std::move(other.coeffs_flat_);
        coeffs_md_ = mdspan_t{coeffs_flat_.data(), other.coeffs_md_.extents()};
    }
    return *this;
}

// Evaluate via Horner's method
template <class Func, std::size_t N_compile>
template <bool SIMD>
typename FuncEvalND<Func, N_compile>::OutputType constexpr FuncEvalND<Func, N_compile>::operator()(
    const InputType & x) const {
    return poly_eval::horner<N_compile, SIMD, OutputType>(map_from_domain(x), coeffs_md_, degree_);
}

// coeff_impl
template <class Func, std::size_t N_compile>
template <typename IdxArray, std::size_t... I>
constexpr typename FuncEvalND<Func, N_compile>::Scalar &
FuncEvalND<Func, N_compile>::coeff_impl(const IdxArray &idx, std::size_t k, std::index_sequence<I...>) noexcept {
    return coeffs_md_(static_cast<std::size_t>(idx[I])..., k);
}

template <class Func, std::size_t N_compile>
template <class IdxArray>
[[nodiscard]] constexpr typename FuncEvalND<Func, N_compile>::Scalar &
FuncEvalND<Func, N_compile>::coeff(const IdxArray &idx, std::size_t k) noexcept {
    return coeff_impl<IdxArray>(idx, k, std::make_index_sequence<dim_>{});
}

template <class Func, std::size_t N_compile> auto FuncEvalND<Func, N_compile>::make_ext(int n) noexcept -> extents_t {
    if constexpr (is_static) {
        return detail::make_static_extents<N_compile, dim_, outDim_>(std::make_index_sequence<dim_>{});
    } else {
        return make_ext(n, std::make_index_sequence<dim_ + 1>{});
    }
}

template <class Func, std::size_t N_compile>
template <std::size_t... Is>
auto FuncEvalND<Func, N_compile>::make_ext(int n, std::index_sequence<Is...>) noexcept -> extents_t {
    return extents_t{(Is < dim_ ? static_cast<std::size_t>(n) : static_cast<std::size_t>(outDim_))...};
}

template <class Func, std::size_t N_compile>
constexpr std::size_t FuncEvalND<Func, N_compile>::storage_required(const int n) noexcept {
    auto ext = make_ext(n);
    auto mapping = typename mdspan_t::mapping_type{ext};
    return mapping.required_span_size();
}

template <class Func, std::size_t N_compile> constexpr void FuncEvalND<Func, N_compile>::initialize(int n) {
    Buffer<Scalar, N_compile> nodes{};
    if constexpr (!N_compile)
        nodes.resize(n);
    for (int k = 0; k < n; ++k)
        nodes[k] = detail::cos((2.0 * double(k) + 1.0) * M_PI / (2.0 * n));

    std::array<int, dim_> ext_idx{};
    ext_idx.fill(n);

    // sample f on Chebyshev grid
    for_each_index<dim_>(ext_idx, [&](const std::array<int, dim_> &idx) {
        InputType x_dom{};
        for (std::size_t d = 0; d < dim_; ++d)
            x_dom[d] = nodes[idx[d]];
        OutputType y = func_(map_to_domain(x_dom));
        for (std::size_t k = 0; k < outDim_; ++k)
            coeff(idx, k) = y[k];
    });

    // convert Newton → monomial along each axis
    Buffer<Scalar, N_compile> rhs{}, alpha{}, mono{};
    if constexpr (!N_compile) {
        rhs.resize(n);
        alpha.resize(n);
        mono.resize(n);
    }

    std::array<int, dim_> base_idx{};
    for (std::size_t axis = 0; axis < dim_; ++axis) {
        auto inner_ext = ext_idx;
        inner_ext[axis] = 1;
        for_each_index<dim_>(inner_ext, [&](const std::array<int, dim_> &base) {
            for (std::size_t k = 0; k < outDim_; ++k) {
                for (int i = 0; i < n; ++i) {
                    base_idx = base;
                    base_idx[axis] = i;
                    rhs[i] = coeff(base_idx, k);
                }
                alpha = detail::bjorck_pereyra<N_compile, Scalar, Scalar>(nodes, rhs);
                mono = detail::newton_to_monomial<N_compile, Scalar, Scalar>(alpha, nodes);
                for (int i = 0; i < n; ++i) {
                    base_idx = base;
                    base_idx[axis] = i;
                    coeff(base_idx, k) = mono[i];
                }
            }
        });
    }

    // reverse coefficient order
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

template <class Func, std::size_t N_compile>
[[nodiscard]] constexpr typename FuncEvalND<Func, N_compile>::InputType
FuncEvalND<Func, N_compile>::map_to_domain(const InputType &t) const noexcept {
    InputType out{};
    for (std::size_t d = 0; d < dim_; ++d)
        out[d] = Scalar(0.5) * (t[d] / low_[d] + hi_[d]);
    return out;
}

template <class Func, std::size_t N_compile>
[[nodiscard]] constexpr typename FuncEvalND<Func, N_compile>::InputType
FuncEvalND<Func, N_compile>::map_from_domain(const InputType &x) const noexcept {
    InputType out{};
    for (std::size_t d = 0; d < dim_; ++d)
        out[d] = (Scalar(2) * x[d] - hi_[d]) * low_[d];
    return out;
}

template <class Func, std::size_t N_compile>
constexpr void FuncEvalND<Func, N_compile>::compute_scaling(const InputType &a, const InputType &b) noexcept {
    for (std::size_t d = 0; d < dim_; ++d) {
        low_[d] = Scalar(1) / (b[d] - a[d]);
        hi_[d] = b[d] + a[d];
    }
}

template <class Func, std::size_t N_compile>
template <std::size_t Rank, class F>
constexpr void FuncEvalND<Func, N_compile>::for_each_index(const std::array<int, Rank> &ext, F &&body) {
    std::array<int, Rank> idx{};
    while (true) {
        body(idx);
        for (std::size_t d = 0; d < Rank; ++d) {
            if (++idx[d] < ext[d])
                break;
            if (d == Rank - 1)
                return;
            idx[d] = 0;
        }
    }
}

// -----------------------------------------------------------------------------
// make_func_eval API implementations (Runtime, C++20 compatible)
// -----------------------------------------------------------------------------
// Compile-time degree (1D or ND)
template <std::size_t N_compile_time, std::size_t Iters_compile_time, class Func>
PF_C20CONSTEXPR auto make_func_eval(Func F, typename function_traits<Func>::arg0_type a,
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
template <std::size_t Iters_compile_time, class Func, typename IntType,
          std::enable_if_t<std::is_integral_v<std::remove_cvref_t<IntType>>, int>>
PF_C20CONSTEXPR auto make_func_eval(Func F, IntType n, typename function_traits<Func>::arg0_type a,
                                   typename function_traits<Func>::arg0_type b,
                                   const std::remove_reference_t<typename function_traits<Func>::arg0_type> *pts) {
    using InputType = typename function_traits<Func>::arg0_type;
    using RawInputType = std::remove_cvref_t<InputType>;

    if constexpr (has_tuple_size_v<RawInputType>) {
        // ND - pass by value for tuple types to avoid reference issues
        return FuncEvalND<Func, 0>(F, static_cast<int>(n), a, b);
    } else {
        // 1D
        return FuncEval<Func, 0, Iters_compile_time>(F, static_cast<int>(n), a, b, pts);
    }
}

template <std::size_t N_compile_time, class Func,
          std::enable_if_t<has_tuple_size_v<std::remove_cvref_t<typename function_traits<Func>::arg0_type>>, int>>
PF_C20CONSTEXPR auto make_func_eval(Func F, typename function_traits<Func>::arg0_type a,
                                   typename function_traits<Func>::arg0_type b) {
    return FuncEvalND<Func, N_compile_time>(F, a, b);
}

// Runtime error tolerance (1D or ND)
template <std::size_t MaxN_val, std::size_t NumEvalPoints_val, std::size_t Iters_compile_time, class Func,
          typename FloatType, std::enable_if_t<std::is_floating_point_v<std::remove_cvref_t<FloatType>>, int>>
PF_C20CONSTEXPR auto make_func_eval(Func F, FloatType eps, typename function_traits<Func>::arg0_type a,
                                   typename function_traits<Func>::arg0_type b) {
    using RawInputType = std::remove_cvref_t<typename function_traits<Func>::arg0_type>;
    using evaluator_t =
        std::conditional_t<has_tuple_size_v<RawInputType>, FuncEvalND<Func, 0>, FuncEval<Func, 0, Iters_compile_time>>;
    const auto eval_pts = detail::linspace(a, b, int(NumEvalPoints_val));
    double max_err;
    for (int n = 2; n <= int(MaxN_val); ++n) {
        const evaluator_t evaluator(F, n, a, b);
        max_err = 0.0;
        for (auto const &pt : eval_pts) {
            const auto actual = F(pt);
            const auto approx = evaluator(pt);
            max_err = std::max(detail::relative_l2_norm(actual, approx), max_err);
        }
        if (max_err <= eps) {
#ifdef DEBUG
            std::cout << "Converged with N=" << n << " (max err=" << std::scientific << max_err << ")\n";
#endif
            return evaluator;
        }
    }
    throw std::runtime_error("No polynomial degree found for requested error tolerance. eps=" + std::to_string(eps) +
                             ", MaxN=" + std::to_string(MaxN_val) + ", max_err=" + std::to_string(max_err));
}

template <std::size_t N_compile_time, std::size_t Iters_compile_time, typename Func,
          std::enable_if_t<std::is_function_v<std::remove_pointer_t<std::decay_t<Func>>>, int>>
PF_C20CONSTEXPR auto make_func_eval(Func *f, typename function_traits<Func *>::arg0_type a,
                                   typename function_traits<Func *>::arg0_type b) {
    using InputType = typename function_traits<Func *>::arg0_type;
    if constexpr (has_tuple_size_v<std::remove_cvref_t<InputType>>) {
        // ND
        auto func_wrapper = [f](const InputType &in) { return f(in); };
        return FuncEvalND<decltype(func_wrapper), N_compile_time>(func_wrapper, a, b);
    } else {
        // 1D
        return FuncEval<Func *, N_compile_time, Iters_compile_time>(f, a, b, nullptr);
    }
}

#if __cplusplus >= 202002L
template <double eps_val, auto a, auto b, std::size_t MaxN_val, std::size_t NumEvalPoints_val,
          std::size_t Iters_compile_time, class Func>
constexpr auto make_func_eval(Func F) {
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
        detail::unroll_loop<2, MaxN_val>([&](const auto I) {
            constexpr auto i = decltype(I)::value;
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
#endif

template <typename... EvalTypes>
PF_C20CONSTEXPR FuncEvalMany<EvalTypes...> make_func_eval_many(EvalTypes... evals) noexcept {
    return FuncEvalMany<std::decay_t<EvalTypes>...>(std::forward<EvalTypes>(evals)...);
}

} // namespace poly_eval