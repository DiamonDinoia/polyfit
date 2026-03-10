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

#include "internal/api_types.hpp"
#include "internal/macros.h"
#include "internal/poly_eval.h"

namespace poly_eval {

template<typename... EvalTypes> class FuncEvalMany;

/**
 * @brief Evaluator for a single polynomial fit of a callable.
 *
 * The evaluator stores coefficients in Horner order, from highest-order term to
 * constant term.
 */
template<class Func, std::size_t NCoeffsCt, std::size_t RefineIters, FusionMode Fusion>
class FuncEval {
  public:
    using InputType = typename function_traits<Func>::arg0_type;
    using OutputType = typename function_traits<Func>::result_type;

    static constexpr std::size_t kNCoeffsCompileTime = NCoeffsCt;
    static constexpr std::size_t kItersCompileTime = RefineIters;

    template<std::size_t CurrentNCoeffs = NCoeffsCt, typename = std::enable_if_t<CurrentNCoeffs != 0>>
    PF_C20CONSTEXPR FuncEval(Func F, InputType a, InputType b, const InputType *pts = nullptr);

    template<std::size_t CurrentNCoeffs = NCoeffsCt, typename = std::enable_if_t<CurrentNCoeffs == 0>>
    PF_C20CONSTEXPR FuncEval(Func F, int nCoeffs, InputType a, InputType b, const InputType *pts = nullptr);

    FuncEval(const FuncEval &) = default;
    FuncEval &operator=(const FuncEval &) = default;
    FuncEval(FuncEval &&) noexcept = default;
    FuncEval &operator=(FuncEval &&) noexcept = default;

    template<bool = false> constexpr OutputType operator()(InputType pt) const noexcept;

    template<bool ptsAligned = false, bool outAligned = false>
    constexpr void operator()(const InputType *pts, OutputType *out, std::size_t numPoints) const noexcept;

    PF_C20CONSTEXPR void truncate(typename value_type_or_identity<OutputType>::type eps) noexcept;

    /**
     * @brief Access coefficients in Horner order.
     *
     * The returned buffer stores coefficients from highest-order term to constant
     * term so it can be passed directly to Horner-based evaluation helpers.
     */
    [[nodiscard]] PF_C20CONSTEXPR const Buffer<OutputType, NCoeffsCt> &coeffs() const noexcept;
    [[nodiscard]] constexpr std::size_t nCoeffs() const noexcept;

  private:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    InputType invSpan, sumEndpoints;
    PF_NO_UNIQUE_ADDRESS std::conditional_t<Fusion == FusionMode::auto_, bool, std::true_type> domainFused{};
    Buffer<OutputType, NCoeffsCt> coeffsBuf;

    PF_C20CONSTEXPR void initializeCoeffs(Func F, const InputType *pts);

    template<class T> PF_ALWAYS_INLINE constexpr T mapToDomain(T value) const noexcept;
    template<class T> PF_ALWAYS_INLINE constexpr T mapFromDomain(T value) const noexcept;

    template<int OuterUnrollFactor, bool ptsAligned, bool outAligned>
    constexpr void evalBatch(const InputType *pts, OutputType *out, std::size_t numPoints) const noexcept;

    PF_C20CONSTEXPR void refine(const Buffer<InputType, NCoeffsCt> &chebNodes,
                                const Buffer<OutputType, NCoeffsCt> &samples);

    template<typename... EvalTypes> friend class FuncEvalMany;
#endif
};

template<typename... EvalTypes> class FuncEvalMany {
    static_assert(sizeof...(EvalTypes) > 0, "At least one FuncEval type is required");
    static constexpr bool hasRuntimeNCoeffs = ((EvalTypes::kNCoeffsCompileTime == 0) || ...);
    static constexpr bool hasCompileTimeNCoeffs = ((EvalTypes::kNCoeffsCompileTime != 0) || ...);
    static_assert(!(hasRuntimeNCoeffs && hasCompileTimeNCoeffs),
                  "Mixing runtime-sized and compile-time-sized evaluators in FuncEvalMany is not supported");

  public:
    using FirstEval = std::tuple_element_t<0, std::tuple<EvalTypes...>>;
    using InputType = typename FirstEval::InputType;
    using OutputType = typename FirstEval::OutputType;

    static constexpr std::size_t kF = sizeof...(EvalTypes);
    static constexpr std::size_t kSimd = xsimd::batch<InputType>::size;
    static constexpr std::size_t kFPad = ((kF + kSimd - 1) / kSimd) * kSimd;
    static constexpr std::size_t vectorWidth = kSimd;
    static constexpr std::size_t kAlignment = xsimd::batch<OutputType>::arch_type::alignment();

    static_assert(!std::is_void_v<xsimd::make_sized_batch_t<InputType, kSimd>>,
                  "SIMD width must be valid for the given type T");

    static constexpr std::size_t maxNCoeffsCt = std::max({EvalTypes::kNCoeffsCompileTime...});

    explicit PF_C20CONSTEXPR FuncEvalMany(const EvalTypes &...evals);
    PF_C20CONSTEXPR FuncEvalMany(const FuncEvalMany &other);
    PF_C20CONSTEXPR FuncEvalMany &operator=(const FuncEvalMany &other);
    PF_C20CONSTEXPR FuncEvalMany(FuncEvalMany &&other) noexcept;
    PF_C20CONSTEXPR FuncEvalMany &operator=(FuncEvalMany &&other) noexcept;

    [[nodiscard]] constexpr std::size_t size() const noexcept;
    [[nodiscard]] constexpr std::size_t nCoeffs() const noexcept;

    std::array<OutputType, kF> operator()(InputType x) const noexcept;
    std::array<OutputType, kF> operator()(const std::array<InputType, kF> &xs) const noexcept;
    void operator()(const InputType *PF_RESTRICT x, OutputType *PF_RESTRICT out, std::size_t numPoints) const noexcept;

    template<typename... Ts> std::array<OutputType, kF> operator()(InputType first, Ts... rest) const noexcept;
    template<typename... Ts> std::array<OutputType, kF> operator()(const std::tuple<Ts...> &tup) const noexcept;

    PF_C20CONSTEXPR void truncate(typename value_type_or_identity<OutputType>::type eps);

  private:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    template<std::size_t I, typename FE, typename... Rest>
    PF_C20CONSTEXPR void copyCoeffs(const FE &eval, const Rest &...rest);

    PF_C20CONSTEXPR void zeroPadCoeffs();
    constexpr std::array<OutputType, kF> extractReal(const std::array<OutputType, kFPad> &full) const noexcept;

    static constexpr std::size_t dyn = stdex::dynamic_extent;
    using Ext = stdex::extents<std::size_t, (maxNCoeffsCt ? maxNCoeffsCt : dyn), kFPad>;

    alignas(kAlignment) AlignedBuffer<OutputType, kFPad * maxNCoeffsCt, kAlignment> coeffStore{};
    stdex::mdspan<OutputType, Ext> coeffs{nullptr, 1, kFPad};
    std::array<InputType, kFPad> invSpan{};
    std::array<InputType, kFPad> sumEndpoints{};
#endif
};

/**
 * @brief Multi-dimensional polynomial evaluator.
 */
template<class Func, std::size_t NCoeffsCt = 0> class FuncEvalND {
  public:
    using Input0 = typename function_traits<Func>::arg0_type;
    using InputType = std::remove_cvref_t<Input0>;
    using OutputType = typename function_traits<Func>::result_type;
    using Scalar = typename OutputType::value_type;

    static constexpr std::size_t dim = std::tuple_size_v<InputType>;
    static constexpr std::size_t outDim = std::tuple_size_v<OutputType>;
    static constexpr bool isStatic = (NCoeffsCt > 0);

    using Extents = std::conditional_t<
        isStatic, decltype(detail::make_static_extents<NCoeffsCt, dim, outDim>(std::make_index_sequence<dim>{})),
        stdex::dextents<std::size_t, dim + 1>>;
    using Mdspan = stdex::mdspan<Scalar, Extents, stdex::layout_right>;

    template<std::size_t C = NCoeffsCt, typename = std::enable_if_t<(C != 0)>>
    constexpr FuncEvalND(Func f, const InputType &a, const InputType &b);

    template<std::size_t C = NCoeffsCt, typename = std::enable_if_t<(C == 0)>>
    constexpr FuncEvalND(Func f, int nCoeffsPerAxis, const InputType &a, const InputType &b);

    template<bool SIMD = true> constexpr OutputType operator()(const InputType &x) const;
    [[nodiscard]] constexpr std::size_t nCoeffsPerAxis() const noexcept;

    FuncEvalND(const FuncEvalND &other);
    FuncEvalND &operator=(const FuncEvalND &other);
    FuncEvalND(FuncEvalND &&other) noexcept;
    FuncEvalND &operator=(FuncEvalND &&other) noexcept;

  private:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    static constexpr std::size_t coeffCount = detail::storage_required<Scalar, NCoeffsCt, dim, outDim>();

    InputType invSpan{}, sumEndpoints{};
    alignas(xsimd::best_arch::alignment())
        AlignedBuffer<Scalar, coeffCount, xsimd::best_arch::alignment()> coeffsFlat;
    Mdspan coeffsMd;

    template<typename IdxArray, std::size_t... I>
    constexpr Scalar &coeffImpl(const IdxArray &idx, std::size_t k, std::index_sequence<I...>) noexcept;

    template<class IdxArray> [[nodiscard]] constexpr Scalar &coeff(const IdxArray &idx, std::size_t k) noexcept;

    static Extents makeExtents(int nCoeffsPerAxis) noexcept;
    template<std::size_t... Is> static Extents makeExtents(int nCoeffsPerAxis, std::index_sequence<Is...>) noexcept;
    static constexpr std::size_t storage_required(int nCoeffsPerAxis) noexcept;

    constexpr void initialize(int nCoeffsPerAxis, Func f);
    constexpr void convertAxesToMonomial(int nCoeffsPerAxis, const Buffer<Scalar, NCoeffsCt> &nodes);
    constexpr void reverseAxes(int nCoeffsPerAxis);
    [[nodiscard]] constexpr InputType mapToDomain(const InputType &x) const noexcept;
    [[nodiscard]] constexpr InputType mapFromDomain(const InputType &x) const noexcept;
    constexpr void computeScaling(const InputType &a, const InputType &b) noexcept;

    template<std::size_t Rank, class F>
    static constexpr void forEachIndex(const std::array<int, Rank> &ext, F &&body);
#endif
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<std::size_t NCoeffsCt, class Func, class... Tags,
         std::enable_if_t<(NCoeffsCt > 0) && detail::all_tags_v<Tags...>, int> = 0>
[[nodiscard]] PF_C20CONSTEXPR auto make_func_eval(Func F, typename function_traits<Func>::arg0_type a,
                                                  typename function_traits<Func>::arg0_type b, Tags...);

template<class Func, typename IntType, class... Tags,
         std::enable_if_t<std::is_integral_v<std::remove_cvref_t<IntType>> && detail::all_tags_v<Tags...>, int> = 0>
[[nodiscard]] PF_C20CONSTEXPR auto
make_func_eval(Func F, IntType nCoeffs, typename function_traits<Func>::arg0_type a,
               typename function_traits<Func>::arg0_type b, Tags...);

template<class Func, typename FloatType, class... Tags,
         std::enable_if_t<std::is_floating_point_v<std::remove_cvref_t<FloatType>> && detail::all_tags_v<Tags...>,
                          int> = 0>
[[nodiscard]] PF_C20CONSTEXPR auto make_func_eval(Func F, FloatType eps, typename function_traits<Func>::arg0_type a,
                                                  typename function_traits<Func>::arg0_type b, Tags...);

template<std::size_t NCoeffsCt, class... Tags, typename Func,
         std::enable_if_t<std::is_function_v<std::remove_pointer_t<std::decay_t<Func>>> && detail::all_tags_v<Tags...>,
                          int> = 0>
[[nodiscard]] PF_C20CONSTEXPR auto make_func_eval(Func *f, typename function_traits<Func *>::arg0_type a,
                                                  typename function_traits<Func *>::arg0_type b, Tags...);
#endif

#if __cplusplus >= 202002L
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<std::size_t NCoeffsCt, auto a, auto b, class Func> [[nodiscard]] constexpr auto make_func_eval(Func F) {
    return FuncEvalND<Func, NCoeffsCt>(F, a, b);
}

#if PF_HAS_CONSTEXPR_EPS_OVERLOAD
template<double epsVal, auto a, auto b, std::size_t maxNCoeffsVal = 32, std::size_t numEvalPointsVal = 100,
         std::size_t refineIters = 1, class Func>
[[nodiscard]] constexpr auto make_func_eval(Func F);
#endif
#endif
#endif

template<typename... EvalTypes>
[[nodiscard]] PF_C20CONSTEXPR FuncEvalMany<EvalTypes...> make_func_eval_many(EvalTypes... evals) noexcept;

} // namespace poly_eval

#include "internal/fast_eval_impl.hpp"
#include "internal/macros_undef.h"
