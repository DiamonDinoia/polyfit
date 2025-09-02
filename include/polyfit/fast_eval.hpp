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

namespace poly_eval {

// Forward declarations
template <typename T> struct function_traits;
template <typename T, typename> struct is_tuple_like;
template <typename T> struct is_func_eval;
template <typename... EvalTypes> class FuncEvalMany;

// -----------------------------------------------------------------------------
// FuncEval: monomial least-squares fit using Chebyshev sampling
// (Runtime or Fixed-Size Compile-Time Storage, but fitting is runtime)
// -----------------------------------------------------------------------------
template <class Func, std::size_t N_compile_time = 0, std::size_t Iters_compile_time = 1> class FuncEval {
  public:
    using InputType = typename function_traits<Func>::arg0_type;
    using OutputType = typename function_traits<Func>::result_type;

    static constexpr std::size_t kDegreeCompileTime = N_compile_time;
    static constexpr std::size_t kItersCompileTime = Iters_compile_time;

    template <std::size_t CurrentN = N_compile_time, typename = std::enable_if_t<CurrentN != 0>>
    PF_C20CONSTEXPR FuncEval(Func F, InputType a, InputType b, const InputType *pts = nullptr);

    template <std::size_t CurrentN = N_compile_time, typename = std::enable_if_t<CurrentN == 0>>
    PF_C20CONSTEXPR FuncEval(Func F, int n, InputType a, InputType b, const InputType *pts = nullptr);

    FuncEval(const FuncEval &) = default;
    FuncEval &operator=(const FuncEval &) = default;
    FuncEval(FuncEval &&) noexcept = default;
    FuncEval &operator=(FuncEval &&) noexcept = default;

    template <bool = false> constexpr OutputType operator()(InputType pt) const noexcept;

    template <bool pts_aligned = false, bool out_aligned = false>
    constexpr void operator()(const InputType *pts, OutputType *out, std::size_t num_points) const noexcept;

    PF_C20CONSTEXPR const Buffer<OutputType, N_compile_time> &coeffs() const noexcept;

  private:
    int n_terms;
    InputType low, hi;
    Buffer<OutputType, N_compile_time> monomials;

    PF_C20CONSTEXPR void initialize_monomials(Func F, const InputType *pts);

    template <class T> PF_ALWAYS_INLINE constexpr T map_to_domain(T T_arg) const noexcept;
    template <class T> PF_ALWAYS_INLINE constexpr T map_from_domain(T T_arg) const noexcept;

    // Evaluate multiple points using SIMD with unrolling
    template <int OuterUnrollFactor, bool pts_aligned, bool out_aligned>
    constexpr void horner_polyeval(const InputType *pts, OutputType *out, std::size_t num_points) const noexcept;

    PF_C20CONSTEXPR void refine(const Buffer<InputType, N_compile_time> &x_cheb_,
                               const Buffer<OutputType, N_compile_time> &y_cheb_);

    // Friend declaration for FuncEvalMany to access private members
    template <typename... EvalTypes> friend class FuncEvalMany;
};

//======================================================================
//  FuncEvalMany – evaluates several FuncEval’s with SIMD-friendly layout
//======================================================================
template <typename... EvalTypes> class FuncEvalMany {
    static_assert(sizeof...(EvalTypes) > 0, "At least one FuncEval type is required");

  public:
    /* Public type aliases */
    using FirstEval = std::tuple_element_t<0, std::tuple<EvalTypes...>>;
    using InputType = typename FirstEval::InputType;
    using OutputType = typename FirstEval::OutputType;

    /* Compile‑time constants */
    static constexpr std::size_t kF = sizeof...(EvalTypes);            // #polynomials
    static constexpr std::size_t kSimd = 1;                            // SIMD width (1 → scalar)
    static constexpr std::size_t kF_pad = kF;                          // padded #polynomials
    static constexpr std::size_t vector_width = kSimd > 1 ? kSimd : 0; // 0 for scalar path

    static_assert(kSimd == 1 || !std::is_void_v<xsimd::make_sized_batch_t<InputType, kSimd>>,
                  "Best SIMD width must be valid for the given type T");

    // Max compile‑time degree across EvalTypes (0 → runtime)
    static constexpr std::size_t deg_max_ctime_ = std::max({EvalTypes::kDegreeCompileTime...});

    /* Construction */
    explicit FuncEvalMany(const EvalTypes &...evals);
    FuncEvalMany(const FuncEvalMany &other);
    FuncEvalMany &operator=(const FuncEvalMany &other);
    FuncEvalMany(FuncEvalMany &&other) noexcept;
    FuncEvalMany &operator=(FuncEvalMany &&other) noexcept;

    /* Introspection */
    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] std::size_t degree() const noexcept;

    /* Evaluation – scalar or per‑poly inputs */
    std::array<OutputType, kF> operator()(InputType x) const noexcept;
    std::array<OutputType, kF> operator()(const std::array<InputType, kF> &xs) const noexcept;
    void operator()(const InputType *x, OutputType *out, std::size_t num_points) const noexcept;

    /* Convenience overloads */
    template <typename... Ts> std::array<OutputType, kF> operator()(InputType first, Ts... rest) const noexcept;

    template <typename... Ts> std::array<OutputType, kF> operator()(const std::tuple<Ts...> &tup) const noexcept;

  private:
    /* Helper routines */
    template <std::size_t I, typename FE, typename... Rest> void copy_coeffs(const FE &fe, const Rest &...rest);

    void zero_pad_coeffs();
    std::array<OutputType, kF> extract_real(const std::array<OutputType, kF_pad> &full) const noexcept;

    /* Storage ---------------------------------------------------------- */
    static constexpr std::size_t dyn = stdex::dynamic_extent;
    using Ext = stdex::extents<std::size_t, (deg_max_ctime_ ? deg_max_ctime_ : dyn), kF_pad>;

    Buffer<OutputType, kF_pad * deg_max_ctime_> coeff_store_{};
    stdex::mdspan<OutputType, Ext> coeffs_{nullptr, 1, kF_pad};

    // Per‑polynomial scaling data
    std::array<InputType, kF_pad> low_{};
    std::array<InputType, kF_pad> hi_{};

    // Runtime max degree (≡ deg_max_ctime_ unless CT value is 0)
    std::size_t deg_max_ = deg_max_ctime_;
};

template <class Func, std::size_t N_compile = 0> class FuncEvalND {
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

    template <std::size_t C = N_compile, typename = std::enable_if_t<(C != 0)>>
    constexpr FuncEvalND(Func f, const InputType &a, const InputType &b);

    template <std::size_t C = N_compile, typename = std::enable_if_t<(C == 0)>>
    constexpr FuncEvalND(Func f, int n, const InputType &a, const InputType &b);

    template <bool SIMD = true> constexpr OutputType operator()(const InputType &x) const;

    FuncEvalND(const FuncEvalND &other);
    FuncEvalND &operator=(const FuncEvalND &other);
    FuncEvalND(FuncEvalND &&other) noexcept;
    FuncEvalND &operator=(FuncEvalND &&other) noexcept;

  private:
    static constexpr std::size_t coeff_count = detail::storage_required<Scalar, N_compile, dim_, outDim_>();

    Func func_;
    int degree_;
    InputType low_{}, hi_{};
    alignas(xsimd::best_arch::alignment())
        AlignedBuffer<Scalar, coeff_count, xsimd::best_arch::alignment()> coeffs_flat_;
    mdspan_t coeffs_md_;

    template <typename IdxArray, std::size_t... I>
    constexpr Scalar &coeff_impl(const IdxArray &idx, std::size_t k, std::index_sequence<I...>) noexcept;

    template <class IdxArray> [[nodiscard]] constexpr Scalar &coeff(const IdxArray &idx, std::size_t k) noexcept;

    static extents_t make_ext(int n) noexcept;
    template <std::size_t... Is> static extents_t make_ext(int n, std::index_sequence<Is...>) noexcept;
    static constexpr std::size_t storage_required(int n) noexcept;

    constexpr void initialize(int n);
    [[nodiscard]] constexpr InputType map_to_domain(const InputType &t) const noexcept;
    [[nodiscard]] constexpr InputType map_from_domain(const InputType &x) const noexcept;
    constexpr void compute_scaling(const InputType &a, const InputType &b) noexcept;

    // utility for multidimensional loops
    template <std::size_t Rank, class F>
    static constexpr void for_each_index(const std::array<int, Rank> &ext, F &&body);
};

// Compile-time degree (1D or ND)
template <std::size_t N_compile_time, std::size_t Iters_compile_time = 1, class Func>
PF_C20CONSTEXPR auto make_func_eval(Func F, typename function_traits<Func>::arg0_type a,
                                   typename function_traits<Func>::arg0_type b,
                                   const typename function_traits<Func>::arg0_type *pts = nullptr);

// Runtime degree (1D or ND, any integral type)
template <std::size_t Iters_compile_time = 1, class Func, typename IntType,
          std::enable_if_t<std::is_integral_v<std::remove_cvref_t<IntType>>, int> = 0>
PF_C20CONSTEXPR auto
make_func_eval(Func F, IntType n, typename function_traits<Func>::arg0_type a,
               typename function_traits<Func>::arg0_type b,
               const std::remove_reference_t<typename function_traits<Func>::arg0_type> *pts = nullptr);

// Runtime error tolerance (1D or ND, any floating-point type)
template <std::size_t MaxN_val = 32, std::size_t NumEvalPoints_val = 100, std::size_t Iters_compile_time = 1,
          class Func, typename FloatType,
          std::enable_if_t<std::is_floating_point_v<std::remove_cvref_t<FloatType>>, int> = 0>
PF_C20CONSTEXPR auto make_func_eval(Func F, FloatType eps, typename function_traits<Func>::arg0_type a,
                                   typename function_traits<Func>::arg0_type b);

template <std::size_t N_compile_time, class Func,
          std::enable_if_t<has_tuple_size_v<std::remove_cvref_t<typename function_traits<Func>::arg0_type>>, int> = 0>
PF_C20CONSTEXPR auto make_func_eval(Func F, typename function_traits<Func>::arg0_type a,
                                   typename function_traits<Func>::arg0_type b);

template <std::size_t N_compile_time, std::size_t Iters_compile_time = 0, typename Func,
          std::enable_if_t<std::is_function_v<std::remove_pointer_t<std::decay_t<Func>>>, int> = 0>
PF_C20CONSTEXPR auto make_func_eval(Func *f, typename function_traits<Func *>::arg0_type a,
                                   typename function_traits<Func *>::arg0_type b);
#if __cplusplus >= 202002L

template <std::size_t N_compile_time, auto a, auto b, class Func> constexpr auto make_func_eval(Func F) {
    // Forward to the appropriate constructor
    return FuncEvalND<Func, N_compile_time>(F, a, b);
}

template <double eps_val, auto a, auto b, std::size_t MaxN_val = 32, std::size_t NumEvalPoints_val = 100,
          std::size_t Iters_compile_time = 1, class Func>
constexpr auto make_func_eval(Func F);
#endif

template <typename... EvalTypes>
PF_C20CONSTEXPR FuncEvalMany<EvalTypes...> make_func_eval_many(EvalTypes... evals) noexcept;
} // namespace poly_eval

// Include implementations
// ReSharper disable once CppUnusedIncludeDirective
#include "internal/fast_eval_impl.hpp"
#include "internal/macros_undef.h"