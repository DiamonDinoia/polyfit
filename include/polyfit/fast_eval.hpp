#pragma once

#include <array>
#include <cmath>
#include <type_traits>

#if __cpp_lib_mdspan >= 202310L
#include <mdspan>
namespace stdex = std;
#else
#include <experimental/mdspan>
namespace stdex = std::experimental;
#endif

#include "internal/macros.h"
#include "internal/poly_eval.h"
#include "internal/tags.h"

//
// Public API overview:
// - make_func_eval: factory producing FuncEval / FuncEvalND evaluators.
// - Prefer compile-time overloads when performance matters (use template N).
// - Adaptive overloads accept an error tolerance and search for a suitable degree.
// - ND (multi-dimensional) functions expect tuple-like / std::array inputs/outputs.
// - See docs/API.md and include/polyfit/internal/macros.h for implementation-only macro docs.
//
// Notes on template parameters:
// - N_compile_time: if >0, degree is fixed at compile time (faster evaluator).
// - Iters_compile_time: number of refinement passes performed after initial fit.
// - The public headers keep implementation details in include/polyfit/internal.

namespace poly_eval {

// Forward declarations
template<typename T> struct function_traits;
template<typename... EvalTypes> class FuncEvalMany;

// -----------------------------------------------------------------------------
// FuncEval: monomial least-squares fit using Chebyshev sampling
// (Runtime or Fixed-Size Compile-Time Storage, but fitting is runtime)
// -----------------------------------------------------------------------------
/**
 * @brief Evaluator for a single polynomial fit of a callable.
 *
 * Performs a (least-squares) polynomial approximation of the provided callable
 * over the domain [a,b] using Chebyshev sampling. The implementation stores
 * coefficients in monomial order and provides fast evaluation (Horner) paths
 * including a bulk SIMD-friendly evaluation overload.
 *
 * Template parameters:
 * @tparam Func Callable type (function, lambda, function pointer, functor).
 * @tparam N_compile_time If non-zero, the degree is fixed at compile-time
 *                        which may enable additional optimizations.
 * @tparam Iters_compile_time Number of refinement iterations performed after
 *                            the initial fit (defaults to 1).
 */
template<class Func, std::size_t N_compile_time, std::size_t Iters_compile_time,
         FusionMode Fusion>
class FuncEval {
  public:
    using InputType = typename function_traits<Func>::arg0_type;
    using OutputType = typename function_traits<Func>::result_type;

    static constexpr std::size_t kDegreeCompileTime = N_compile_time;
    static constexpr std::size_t kItersCompileTime = Iters_compile_time;

    /**
     * @brief Construct a compile-time degree evaluator.
     *
     * @param F Callable to approximate.
     * @param a Domain low endpoint.
     * @param b Domain high endpoint.
     * @param pts Optional pointer to evaluation points used for fitting.
     */
    template<std::size_t CurrentN = N_compile_time, typename = std::enable_if_t<CurrentN != 0>>
    PF_C20CONSTEXPR FuncEval(Func F, InputType a, InputType b, const InputType *pts = nullptr);

    /**
     * @brief Construct a runtime-degree evaluator.
     *
     * @param F Callable to approximate.
     * @param n Requested polynomial degree (runtime).
     * @param a Domain low endpoint.
     * @param b Domain high endpoint.
     * @param pts Optional pointer to evaluation points used for fitting.
     */
    template<std::size_t CurrentN = N_compile_time, typename = std::enable_if_t<CurrentN == 0>>
    PF_C20CONSTEXPR FuncEval(Func F, int n, InputType a, InputType b, const InputType *pts = nullptr);

    FuncEval(const FuncEval &) = default;
    FuncEval &operator=(const FuncEval &) = default;
    FuncEval(FuncEval &&) noexcept = default;
    FuncEval &operator=(FuncEval &&) noexcept = default;

    /**
     * @brief Evaluate the fitted polynomial at a single input point.
     *
     * @param pt Input value in the original domain [a,b].
     * @return OutputType Approximated value of the callable at `pt`.
     */
    template<bool = false> constexpr OutputType operator()(InputType pt) const noexcept;

    /**
     * @brief Bulk evaluation of multiple input points.
     *
     * This overload writes `num_points` outputs to `out` for the `pts`
     * array. `pts_aligned` and `out_aligned` hint alignment for SIMD paths.
     *
     * @tparam pts_aligned True if `pts` is SIMD-aligned.
     * @tparam out_aligned True if `out` is SIMD-aligned.
     * @param pts Pointer to input values.
     * @param out Pointer to output buffer (must have room for `num_points`).
     * @param num_points Number of points to evaluate.
     */
    template<bool pts_aligned = false, bool out_aligned = false>
    constexpr void operator()(const InputType *pts, OutputType *out, std::size_t num_points) const noexcept;

    /**
     * @brief Truncate negligible high-degree coefficients.
     *
     * Scans leading (high-degree) coefficients and removes those whose
     * absolute value is below `eps`. Only available for runtime-degree
     * evaluators (N_compile_time == 0).
     *
     * @param eps Threshold below which coefficients are considered negligible.
     */
    PF_C20CONSTEXPR void truncate(typename value_type_or_identity<OutputType>::type eps) noexcept;

    /**
     * @brief Access the fitted coefficients in monomial order.
     *
     * The returned `Buffer` contains coefficients ordered from the constant
     * term up to the highest-degree term used by this evaluator.
     */
    PF_C20CONSTEXPR const Buffer<OutputType, N_compile_time> &coeffs() const noexcept;

  private:
    InputType low, hi;
    PF_NO_UNIQUE_ADDRESS std::conditional_t<Fusion == FusionMode::auto_, bool, std::true_type> fused_{};
    Buffer<OutputType, N_compile_time> monomials;

    PF_C20CONSTEXPR void initialize_monomials(Func F, const InputType *pts);

    template<class T> PF_ALWAYS_INLINE constexpr T map_to_domain(T T_arg) const noexcept;
    template<class T> PF_ALWAYS_INLINE constexpr T map_from_domain(T T_arg) const noexcept;

    // Evaluate multiple points using SIMD with unrolling
    template<int OuterUnrollFactor, bool pts_aligned, bool out_aligned>
    constexpr void horner_polyeval(const InputType *pts, OutputType *out, std::size_t num_points) const noexcept;

    PF_C20CONSTEXPR void refine(const Buffer<InputType, N_compile_time> &x_cheb_,
                                const Buffer<OutputType, N_compile_time> &y_cheb_);

    // Friend declaration for FuncEvalMany to access private members
    template<typename... EvalTypes> friend class FuncEvalMany;
};

//======================================================================
//  FuncEvalMany – evaluates several FuncEval’s with SIMD-friendly layout
//======================================================================
/**
 * @brief Container evaluating multiple polynomial evaluators in a SIMD-friendly layout.
 *
 * `FuncEvalMany` packs coefficients for several `FuncEval` instances into a
 * layout amenable to vectorized evaluation. Use this to evaluate multiple
 * polynomials for the same input(s) efficiently.
 *
 * Template parameters:
 * @tparam EvalTypes List of `FuncEval<...>` types to pack together.
 *
 * Notes:
 * - `kF` is the number of polynomials and `kF_pad` is the padded storage width.
 * - The class exposes scalar and bulk evaluation overloads.
 */
template<typename... EvalTypes> class FuncEvalMany {
    static_assert(sizeof...(EvalTypes) > 0, "At least one FuncEval type is required");
    static constexpr bool has_runtime_degree_ = ((EvalTypes::kDegreeCompileTime == 0) || ...);
    static constexpr bool has_compile_time_degree_ = ((EvalTypes::kDegreeCompileTime != 0) || ...);
    static_assert(!(has_runtime_degree_ && has_compile_time_degree_),
                  "Mixing runtime-degree and compile-time-degree evaluators in FuncEvalMany is not supported");

  public:
    /* Public type aliases */
    using FirstEval = std::tuple_element_t<0, std::tuple<EvalTypes...>>;
    using InputType = typename FirstEval::InputType;
    using OutputType = typename FirstEval::OutputType;

    /* Compile‑time constants */
    static constexpr std::size_t kF = sizeof...(EvalTypes); // #polynomials
    static constexpr std::size_t kSimd = xsimd::batch<InputType>::size;
    static constexpr std::size_t kF_pad = ((kF + kSimd - 1) / kSimd) * kSimd;
    static constexpr std::size_t vector_width = kSimd;
    static constexpr std::size_t kAlignment = xsimd::batch<OutputType>::arch_type::alignment();

    static_assert(!std::is_void_v<xsimd::make_sized_batch_t<InputType, kSimd>>,
                  "SIMD width must be valid for the given type T");

    // Max compile‑time degree across EvalTypes (0 → runtime)
    static constexpr std::size_t deg_max_ctime_ = std::max({EvalTypes::kDegreeCompileTime...});

    /* Construction */
    explicit PF_C20CONSTEXPR FuncEvalMany(const EvalTypes &...evals);
    PF_C20CONSTEXPR FuncEvalMany(const FuncEvalMany &other);
    PF_C20CONSTEXPR FuncEvalMany &operator=(const FuncEvalMany &other);
    PF_C20CONSTEXPR FuncEvalMany(FuncEvalMany &&other) noexcept;
    PF_C20CONSTEXPR FuncEvalMany &operator=(FuncEvalMany &&other) noexcept;

    /* Introspection */
    [[nodiscard]] constexpr std::size_t size() const noexcept;
    [[nodiscard]] constexpr std::size_t degree() const noexcept;

    /* Evaluation – scalar or per‑poly inputs */
    std::array<OutputType, kF> operator()(InputType x) const noexcept;
    std::array<OutputType, kF> operator()(const std::array<InputType, kF> &xs) const noexcept;
    void operator()(const InputType *PF_RESTRICT x, OutputType *PF_RESTRICT out, std::size_t num_points) const noexcept;

    /* Convenience overloads */
    template<typename... Ts> std::array<OutputType, kF> operator()(InputType first, Ts... rest) const noexcept;

    template<typename... Ts> std::array<OutputType, kF> operator()(const std::tuple<Ts...> &tup) const noexcept;

    /**
     * @brief Truncate negligible high-degree coefficient rows.
     *
     * Scans from the highest-degree row downward and removes rows where
     * all coefficients are below `eps`. After truncation, evaluation uses
     * runtime-degree Horner.
     *
     * @param eps Threshold below which coefficient rows are considered negligible.
     */
    PF_C20CONSTEXPR void truncate(typename value_type_or_identity<OutputType>::type eps);

  private:
    /* Helper routines */
    template<std::size_t I, typename FE, typename... Rest>
    PF_C20CONSTEXPR void copy_coeffs(const FE &eval, const Rest &...rest);

    PF_C20CONSTEXPR void zero_pad_coeffs();
    constexpr std::array<OutputType, kF> extract_real(const std::array<OutputType, kF_pad> &full) const noexcept;

    /* Storage ---------------------------------------------------------- */
    static constexpr std::size_t dyn = stdex::dynamic_extent;
    using Ext = stdex::extents<std::size_t, (deg_max_ctime_ ? deg_max_ctime_ : dyn), kF_pad>;

    alignas(kAlignment) AlignedBuffer<OutputType, kF_pad * deg_max_ctime_, kAlignment> coeff_store_{};
    stdex::mdspan<OutputType, Ext> coeffs_{nullptr, 1, kF_pad};

    // Per‑polynomial scaling data
    std::array<InputType, kF_pad> low_{};
    std::array<InputType, kF_pad> hi_{};

};

/**
 * @brief Multi-dimensional polynomial evaluator.
 *
 * `FuncEvalND` extends `FuncEval` style fitting/evaluation to functions whose
 * input is a tuple/array (multi-dimensional). The evaluator stores a flat
 * coefficient buffer and exposes `operator()` that maps tuple inputs to
 * vector outputs.
 *
 * Template parameters:
 * @tparam Func Callable type accepting a tuple-like input.
 * @tparam N_compile If non-zero, the per-dimension degree is fixed at compile time.
 */
template<class Func, std::size_t N_compile = 0> class FuncEvalND {
  public:
    using Input0 = typename function_traits<Func>::arg0_type;
    using InputType = std::remove_cvref_t<Input0>;
    using OutputType = typename function_traits<Func>::result_type;
    using Scalar = typename OutputType::value_type;

    static constexpr std::size_t dim_ = std::tuple_size_v<InputType>;
    static constexpr std::size_t outDim_ = std::tuple_size_v<OutputType>;
    static constexpr bool is_static = (N_compile > 0);

    using extents_t = std::conditional_t<
        is_static, decltype(detail::make_static_extents<N_compile, dim_, outDim_>(std::make_index_sequence<dim_>{})),
        stdex::dextents<std::size_t, dim_ + 1>>;
    using mdspan_t = stdex::mdspan<Scalar, extents_t, stdex::layout_right>;

    template<std::size_t C = N_compile, typename = std::enable_if_t<(C != 0)>>
    constexpr FuncEvalND(Func f, const InputType &a, const InputType &b);

    template<std::size_t C = N_compile, typename = std::enable_if_t<(C == 0)>>
    constexpr FuncEvalND(Func f, int n, const InputType &a, const InputType &b);

    template<bool SIMD = true> constexpr OutputType operator()(const InputType &x) const;

    FuncEvalND(const FuncEvalND &other);
    FuncEvalND &operator=(const FuncEvalND &other);
    FuncEvalND(FuncEvalND &&other) noexcept;
    FuncEvalND &operator=(FuncEvalND &&other) noexcept;

  private:
    static constexpr std::size_t coeff_count = detail::storage_required<Scalar, N_compile, dim_, outDim_>();

    InputType low_{}, hi_{};
    alignas(xsimd::best_arch::alignment())
        AlignedBuffer<Scalar, coeff_count, xsimd::best_arch::alignment()> coeffs_flat_;
    mdspan_t coeffs_md_;

    template<typename IdxArray, std::size_t... I>
    constexpr Scalar &coeff_impl(const IdxArray &idx, std::size_t k, std::index_sequence<I...>) noexcept;

    template<class IdxArray> [[nodiscard]] constexpr Scalar &coeff(const IdxArray &idx, std::size_t k) noexcept;

    static extents_t make_ext(int n) noexcept;
    template<std::size_t... Is> static extents_t make_ext(int n, std::index_sequence<Is...>) noexcept;
    static constexpr std::size_t storage_required(int n) noexcept;

    constexpr void initialize(int n, Func f);
    constexpr void convert_newton_to_monomial_axes(int n, const Buffer<Scalar, N_compile> &nodes);
    constexpr void reverse_coefficients_axes(int n);
    [[nodiscard]] constexpr InputType map_to_domain(const InputType &t) const noexcept;
    [[nodiscard]] constexpr InputType map_from_domain(const InputType &x) const noexcept;
    constexpr void compute_scaling(const InputType &a, const InputType &b) noexcept;

    // utility for multidimensional loops
    template<std::size_t Rank, class F>
    static constexpr void for_each_index(const std::array<int, Rank> &ext, F &&body);
};

// Compile-time degree (1D or ND) — optional trailing tags: iters<N>{}, fusion::*{}
template<std::size_t N_compile_time, class Func, class... Tags,
         std::enable_if_t<(N_compile_time > 0) && detail::all_tags_v<Tags...>, int> = 0>
[[nodiscard]] PF_C20CONSTEXPR auto make_func_eval(Func F, typename function_traits<Func>::arg0_type a,
                                                  typename function_traits<Func>::arg0_type b, Tags...);

// Runtime degree (1D or ND) — optional trailing tags: iters<N>{}, fusion::*{}
template<class Func, typename IntType, class... Tags,
         std::enable_if_t<std::is_integral_v<std::remove_cvref_t<IntType>> && detail::all_tags_v<Tags...>, int> = 0>
[[nodiscard]] PF_C20CONSTEXPR auto
make_func_eval(Func F, IntType n, typename function_traits<Func>::arg0_type a,
               typename function_traits<Func>::arg0_type b, Tags...);

// Runtime error tolerance (1D or ND) — optional trailing tags: max_degree<N>{}, eval_pts<N>{}, iters<N>{}, fusion::*{}
template<class Func, typename FloatType, class... Tags,
         std::enable_if_t<std::is_floating_point_v<std::remove_cvref_t<FloatType>> && detail::all_tags_v<Tags...>,
                          int> = 0>
[[nodiscard]] PF_C20CONSTEXPR auto make_func_eval(Func F, FloatType eps, typename function_traits<Func>::arg0_type a,
                                                  typename function_traits<Func>::arg0_type b, Tags...);

// Function pointer overload — optional trailing tags
template<std::size_t N_compile_time, class... Tags, typename Func,
         std::enable_if_t<std::is_function_v<std::remove_pointer_t<std::decay_t<Func>>> && detail::all_tags_v<Tags...>,
                          int> = 0>
[[nodiscard]] PF_C20CONSTEXPR auto make_func_eval(Func *f, typename function_traits<Func *>::arg0_type a,
                                                  typename function_traits<Func *>::arg0_type b, Tags...);

#if __cplusplus >= 202002L

template<std::size_t N_compile_time, auto a, auto b, class Func> [[nodiscard]] constexpr auto make_func_eval(Func F) {
    // Forward to the appropriate constructor
    return FuncEvalND<Func, N_compile_time>(F, a, b);
}

#if PF_HAS_CONSTEXPR_EPS_OVERLOAD
template<double eps_val, auto a, auto b, std::size_t MaxN_val = 32, std::size_t NumEvalPoints_val = 100,
         std::size_t Iters_compile_time = 1, class Func>
[[nodiscard]] constexpr auto make_func_eval(Func F);
#endif // PF_HAS_CONSTEXPR_EPS_OVERLOAD
#endif

/**
 * @brief Pack multiple evaluators into a SIMD-friendly `FuncEvalMany`.
 *
 * @param evals Evaluator instances to pack.
 * @return FuncEvalMany<EvalTypes...> Packed evaluator object.
 */
template<typename... EvalTypes>
[[nodiscard]] PF_C20CONSTEXPR FuncEvalMany<EvalTypes...> make_func_eval_many(EvalTypes... evals) noexcept;
} // namespace poly_eval

// Include implementations
// ReSharper disable once CppUnusedIncludeDirective
#include "internal/fast_eval_impl.hpp"
#include "internal/macros_undef.h"
