#pragma once

#include <array>
#include <cmath>
#include <type_traits>
#if defined(__has_include)
#  if __has_include(<span>)
#    include <span>
#  endif
#endif

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
template<class Func, std::size_t NCOEFFS_CT, std::size_t ITERS_CT, FusionMode FUSION_MODE>
class FuncEval {
    using InputScalar = detail::value_type_or_t<fitInput_t<Func>>;
    using OutputScalar = detail::value_type_or_t<fitOutput_t<Func>>;
    static_assert(std::is_same_v<InputScalar, OutputScalar>,
                  "FuncEval requires matching input and output scalar types; for example double -> complex<double> is "
                  "supported, but float -> double is not");

  public:
    using InputType = fitInput_t<Func>;
    using OutputType = fitOutput_t<Func>;
    using InputBuffer = Buffer<InputType, NCOEFFS_CT>;
    using OutputBuffer = Buffer<OutputType, NCOEFFS_CT>;

    static constexpr std::size_t NCOEFFS = NCOEFFS_CT;
    static constexpr std::size_t ITERS = ITERS_CT;

    PF_CXX20_CONSTEXPR FuncEval(Func F, InputType a, InputType b, const InputType *pts = nullptr);

    PF_CXX20_CONSTEXPR FuncEval(Func F, int nCoeffs, InputType a, InputType b, const InputType *pts = nullptr);

    FuncEval(const FuncEval &) = default;
    FuncEval &operator=(const FuncEval &) = default;
    FuncEval(FuncEval &&) noexcept = default;
    FuncEval &operator=(FuncEval &&) noexcept = default;

    template<bool = false> constexpr OutputType operator()(InputType pt) const noexcept;

    template<bool = false, class V>
    constexpr auto operator()(V pt) const noexcept
        -> enable_if_t<!std::is_same_v<remove_cvref_t<V>, InputType> &&
                           std::is_constructible_v<remove_cvref_t<V>, OutputType>,
                       remove_cvref_t<V>>;

    template<bool ptsAligned = false, bool outAligned = false>
    constexpr void operator()(const InputType *pts, OutputType *out, std::size_t numPoints) const noexcept;

    PF_CXX20_CONSTEXPR void truncate(detail::value_type_or_t<OutputType> eps) noexcept;

    /**
     * @brief Access coefficients in Horner order.
     *
     * The returned buffer stores coefficients from highest-order term to constant
     * term so it can be passed directly to Horner-based evaluation helpers.
     */
    [[nodiscard]] PF_CXX20_CONSTEXPR const OutputBuffer &coeffs() const noexcept;
    [[nodiscard]] constexpr std::size_t nCoeffs() const noexcept;

  private:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    struct DomainParams {
        InputType invSpan{}, sumEndpoints{};
        bool identityMap = false;
    };
    static constexpr bool kStoresDomain = (FUSION_MODE != FusionMode::Always);
    struct EmptyDomain {};
    PF_NO_UNIQUE_ADDRESS std::conditional_t<kStoresDomain, DomainParams, EmptyDomain> domain_;
    OutputBuffer coeffsBuf;

    constexpr InputType domainInvSpan() const noexcept {
        if constexpr (kStoresDomain) return domain_.invSpan;
        else return InputType(0.5);
    }
    constexpr InputType domainSumEndpoints() const noexcept {
        if constexpr (kStoresDomain) return domain_.sumEndpoints;
        else return InputType(0);
    }
    constexpr bool domainIsIdentity() const noexcept {
        if constexpr (kStoresDomain) return domain_.identityMap;
        else return true;
    }

    PF_CXX20_CONSTEXPR void initialize(detail::CompileTimeCountTag, Func F, InputType a, InputType b,
                                    const InputType *pts, DomainParams &dp);
    PF_CXX20_CONSTEXPR void initialize(detail::RuntimeCountTag, Func F, int nCoeffs, InputType a, InputType b,
                                    const InputType *pts, DomainParams &dp);

    PF_CXX20_CONSTEXPR void buildNodeGrid(InputBuffer &grid, const InputType *pts) const;
    PF_CXX20_CONSTEXPR void sampleOnGrid(OutputBuffer &samples, const InputBuffer &grid, Func F,
                                       const DomainParams &dp) const;
    PF_CXX20_CONSTEXPR void computeMonomialCoeffs(const InputBuffer &grid, const OutputBuffer &samples);
    [[nodiscard]] PF_CXX20_CONSTEXPR bool shouldFuseDomain(const DomainParams &dp) const noexcept;
    PF_CXX20_CONSTEXPR void fuseDomain(DomainParams &dp);

    PF_CXX20_CONSTEXPR void initializeCoeffs(Func F, const InputType *pts, DomainParams &dp);

    template<class T> static PF_ALWAYS_INLINE constexpr T mapToDomain(const DomainParams &dp, T value) noexcept;
    template<class T> PF_ALWAYS_INLINE constexpr T mapFromDomain(T value) const noexcept;

    template<int OuterUnrollFactor, bool ptsAligned, bool outAligned>
    constexpr void evalBatch(const InputType *pts, OutputType *out, std::size_t numPoints) const noexcept;

    PF_CXX20_CONSTEXPR void refine(const InputBuffer &chebNodes, const OutputBuffer &samples);

    template<typename... EvalTypes> friend class FuncEvalMany;
#endif
};

template<typename... EvalTypes> class FuncEvalMany {
    static_assert(sizeof...(EvalTypes) > 0, "At least one FuncEval type is required");
    static constexpr bool hasRuntimeNCoeffs = ((EvalTypes::NCOEFFS == 0) || ...);
    static constexpr bool hasCompileTimeNCoeffs = ((EvalTypes::NCOEFFS != 0) || ...);
    static_assert(!(hasRuntimeNCoeffs && hasCompileTimeNCoeffs),
                  "Mixing runtime-sized and compile-time-sized evaluators in FuncEvalMany is not supported");

  public:
    using FirstEval = std::tuple_element_t<0, std::tuple<EvalTypes...>>;
    using InputType = typename FirstEval::InputType;
    using OutputType = typename FirstEval::OutputType;
    static_assert((std::is_same_v<InputType, typename EvalTypes::InputType> && ...),
                  "All evaluators packed into FuncEvalMany must have the same input type");
    static_assert((std::is_same_v<OutputType, typename EvalTypes::OutputType> && ...),
                  "All evaluators packed into FuncEvalMany must have the same output type");

    static constexpr std::size_t COUNT = sizeof...(EvalTypes);

  private:
    static constexpr std::size_t SIMD_WIDTH = xsimd::batch<InputType>::size;
    static constexpr std::size_t PADDED_COUNT = ((COUNT + SIMD_WIDTH - 1) / SIMD_WIDTH) * SIMD_WIDTH;
    static constexpr std::size_t VECTOR_WIDTH = SIMD_WIDTH;
    static constexpr std::size_t ALIGNMENT = xsimd::batch<OutputType>::arch_type::alignment();

    static_assert(!std::is_void_v<xsimd::make_sized_batch_t<InputType, SIMD_WIDTH>>,
                  "SIMD width must be valid for the given type T");

    static constexpr std::size_t MAX_NCOEFFS = std::max({EvalTypes::NCOEFFS...});

  public:

    explicit PF_CXX20_CONSTEXPR FuncEvalMany(const EvalTypes &...evals);
    PF_CXX20_CONSTEXPR FuncEvalMany(const FuncEvalMany &other);
    PF_CXX20_CONSTEXPR FuncEvalMany &operator=(const FuncEvalMany &other);
    PF_CXX20_CONSTEXPR FuncEvalMany(FuncEvalMany &&other) noexcept;
    PF_CXX20_CONSTEXPR FuncEvalMany &operator=(FuncEvalMany &&other) noexcept;

    [[nodiscard]] constexpr std::size_t size() const noexcept;
    [[nodiscard]] constexpr std::size_t nCoeffs() const noexcept;

    /**
     * @brief Access the coefficient at (coeffIndex, polyIndex) in Horner order.
     */
    [[nodiscard]] constexpr const OutputType &coeff(std::size_t coeffIndex, std::size_t polyIndex) const noexcept;

    std::array<OutputType, COUNT> operator()(InputType x) const noexcept;
    std::array<OutputType, COUNT> operator()(const std::array<InputType, COUNT> &xs) const noexcept;
    void operator()(const InputType *PF_RESTRICT x, OutputType *PF_RESTRICT out, std::size_t numPoints) const noexcept;

    template<typename... Ts> std::array<OutputType, COUNT> operator()(InputType first, Ts... rest) const noexcept;
    template<typename... Ts> std::array<OutputType, COUNT> operator()(const std::tuple<Ts...> &tup) const noexcept;

    PF_CXX20_CONSTEXPR void truncate(detail::value_type_or_t<OutputType> eps);

  private:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    template<class Step> constexpr void forEachCoeff(Step &&step) const noexcept;
    PF_CXX20_CONSTEXPR void rebindCoeffs(std::size_t coeffCount) noexcept;
    PF_CXX20_CONSTEXPR void bindCoeffView(std::size_t coeffCount);
    [[nodiscard]] constexpr OutputType &coeffRef(std::size_t coeffIndex, std::size_t polyIndex) noexcept;
    [[nodiscard]] constexpr const OutputType &coeffRef(std::size_t coeffIndex, std::size_t polyIndex) const noexcept;
    template<class Get> [[nodiscard]] constexpr std::array<InputType, PADDED_COUNT> gatherInputs(Get &&get) const noexcept {
        std::array<InputType, PADDED_COUNT> inputs{};
        poet::static_for<COUNT>([&](auto i) { inputs[i] = get(std::size_t(i)); });
        return inputs;
    }
    template<class Get> [[nodiscard]] constexpr std::array<InputType, PADDED_COUNT> mapInputs(Get &&get) const noexcept {
        auto inputs = gatherInputs(std::forward<Get>(get));
        poet::static_for<COUNT>([&](auto i) { inputs[i] = mapInput(std::size_t(i), inputs[i]); });
        return inputs;
    }
    [[nodiscard]] constexpr InputType mapInput(std::size_t polyIndex, InputType x) const noexcept {
        if (identityMap[polyIndex]) return x;
        return xsimd::fms(InputType(2.0), x, sumEndpoints[polyIndex]) * invSpan[polyIndex];
    }
    [[nodiscard]] constexpr std::array<InputType, PADDED_COUNT> mapInputs(InputType x) const noexcept {
        return mapInputs([=](std::size_t) { return x; });
    }
    [[nodiscard]] constexpr std::array<InputType, PADDED_COUNT> mapInputs(const std::array<InputType, COUNT> &xs) const noexcept {
        return mapInputs([&](std::size_t i) { return xs[i]; });
    }
    [[nodiscard]] std::array<OutputType, COUNT> evalInputs(const std::array<InputType, PADDED_COUNT> &xu) const noexcept;
    void scatterColumnBatch(xsimd::batch<OutputType> acc, OutputType *out, std::size_t base, std::size_t column) const noexcept;
    void evalColumn(std::size_t column, const InputType *x, OutputType *out, std::size_t numPoints) const noexcept;

    template<std::size_t I, typename FE, typename... Rest>
    PF_CXX20_CONSTEXPR void copyCoeffs(const FE &eval, const Rest &...rest);

    PF_CXX20_CONSTEXPR void zeroPadCoeffs();
    [[nodiscard]] constexpr std::array<OutputType, COUNT> extractReal(const std::array<OutputType, PADDED_COUNT> &full) const noexcept;

    static constexpr std::size_t dyn = stdex::dynamic_extent;
    using Ext = stdex::extents<std::size_t, (MAX_NCOEFFS ? MAX_NCOEFFS : dyn), PADDED_COUNT>;

    alignas(ALIGNMENT) AlignedBuffer<OutputType, PADDED_COUNT * MAX_NCOEFFS, ALIGNMENT> coeffStore{};
    stdex::mdspan<OutputType, Ext> coeffs{nullptr, 1, PADDED_COUNT};
    std::array<InputType, PADDED_COUNT> invSpan{};
    std::array<InputType, PADDED_COUNT> sumEndpoints{};
    std::array<bool, PADDED_COUNT> identityMap{};
    bool allIdentityMap = false;
#endif
};

/**
 * @brief Multi-dimensional polynomial evaluator.
 */
template<class Func, std::size_t NCOEFFS = 0> class FuncEvalND {
  public:
    using InputType = poly_eval::remove_cvref_t<fitInput_t<Func>>;
    using OutputType = poly_eval::remove_cvref_t<fitOutput_t<Func>>;
    using InputScalar = detail::index_value_t<InputType>;
    using Scalar = detail::index_value_t<OutputType>;

    static constexpr std::size_t DIM = detail::fixed_size_v<InputType>;
    static constexpr std::size_t OUT_DIM = detail::fixed_size_v<OutputType>;
    static constexpr bool IS_STATIC = (NCOEFFS > 0);
    using CanonicalInput = std::array<InputScalar, DIM>;
    using CanonicalOutput = std::array<Scalar, OUT_DIM>;

    using Extents = std::conditional_t<
        IS_STATIC, decltype(detail::makeStaticExtents<NCOEFFS, DIM, OUT_DIM>(std::make_index_sequence<DIM>{})),
        stdex::dextents<std::size_t, DIM + 1>>;
    using Mdspan = stdex::mdspan<Scalar, Extents, stdex::layout_right>;
    template<class... Coords> static constexpr bool acceptsCoords_v = detail::isCompatibleNdCoordPack_v<DIM, InputScalar, Coords...>;
    template<class Points, class Outputs>
    static constexpr bool acceptsBatchContainers_v =
        detail::isCompatibleNdBatchContainers_v<Points, Outputs, DIM, InputScalar, OutputType>;

    constexpr FuncEvalND(Func f, const InputType &a, const InputType &b);

    constexpr FuncEvalND(Func f, int nCoeffsPerAxis, const InputType &a, const InputType &b);

    template<bool SIMD = true> PF_FLATTEN constexpr OutputType operator()(const InputType &x) const;
    template<bool SIMD = true, class Point,
             class = enable_if_t<!std::is_same_v<remove_cvref_t<Point>, InputType> &&
                                  detail::isCompatibleNdPoint_v<Point, DIM, InputScalar>>>
    PF_FLATTEN constexpr OutputType operator()(const Point &x) const;
    template<class... Coords, class = enable_if_t<acceptsCoords_v<Coords...>>>
    PF_FLATTEN constexpr OutputType operator()(Coords... coords) const;
    constexpr void operator()(const CanonicalInput *pts, CanonicalOutput *out, std::size_t count) const noexcept;
#if defined(__cpp_lib_span) && (__cpp_lib_span >= 202002L)
    PF_FLATTEN constexpr void operator()(std::span<const CanonicalInput> pts, std::span<CanonicalOutput> out) const;
#endif
    template<class Points, class Outputs, class = enable_if_t<acceptsBatchContainers_v<Points, Outputs>>>
    PF_FLATTEN constexpr void operator()(const Points &pts, Outputs &out) const;
    [[nodiscard]] constexpr std::size_t nCoeffsPerAxis() const noexcept;

    /**
     * @brief Access the coefficient at spatial index @p idx for output dimension @p k.
     *
     * Coefficients are in Horner order (highest-degree first along each axis).
     */
    template<class IdxArray>
    [[nodiscard]] constexpr const Scalar &coeff(const IdxArray &idx, std::size_t k) const noexcept;

    FuncEvalND(const FuncEvalND &other);
    FuncEvalND &operator=(const FuncEvalND &other);
    FuncEvalND(FuncEvalND &&other) noexcept;
    FuncEvalND &operator=(FuncEvalND &&other) noexcept;

  private:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    static constexpr std::size_t COEFF_STORAGE = detail::storageRequired<Scalar, NCOEFFS, DIM, OUT_DIM>();

    CanonicalInput invSpan{}, sumEndpoints{};
    bool identityDomain = false;
    alignas(xsimd::best_arch::alignment())
        AlignedBuffer<Scalar, COEFF_STORAGE, xsimd::best_arch::alignment()> coeffsFlat;
    Mdspan coeffsMd;

    constexpr void initialize(detail::CompileTimeCountTag, Func f, const InputType &a, const InputType &b);
    constexpr void initialize(detail::RuntimeCountTag, Func f, int nCoeffsPerAxis, const InputType &a,
                              const InputType &b);

    template<typename IdxArray, std::size_t... I>
    constexpr Scalar &coeffImpl(const IdxArray &idx, std::size_t k, std::index_sequence<I...>) noexcept;
    template<typename IdxArray, std::size_t... I>
    constexpr const Scalar &coeffImpl(const IdxArray &idx, std::size_t k, std::index_sequence<I...>) const noexcept;

    template<class IdxArray> [[nodiscard]] constexpr Scalar &coeff(const IdxArray &idx, std::size_t k) noexcept;

    static Extents makeExtents(int nCoeffsPerAxis) noexcept;
    template<std::size_t... Is> static Extents makeExtents(int nCoeffsPerAxis, std::index_sequence<Is...>) noexcept;
    static constexpr std::size_t storageRequired(int nCoeffsPerAxis) noexcept;

    constexpr void buildCoeffs(int nCoeffsPerAxis, Func f);
    constexpr void convertAxesToMonomialBasis(int nCoeffsPerAxis, const Buffer<Scalar, NCOEFFS> &nodes);
    constexpr void reverseCoefficientOrder(int nCoeffsPerAxis);
    template<bool SIMD, class Point> [[nodiscard]] constexpr OutputType evalPoint(const Point &x) const;
    template<class Point> [[nodiscard]] static constexpr CanonicalInput toCanonicalInput(const Point &x) noexcept;
    [[nodiscard]] static constexpr InputType fromCanonicalInput(const CanonicalInput &x) noexcept;
    [[nodiscard]] static constexpr CanonicalOutput toCanonicalOutput(const OutputType &x) noexcept;
    [[nodiscard]] static constexpr OutputType fromCanonicalOutput(const CanonicalOutput &x) noexcept;
    template<bool SIMD = true> [[nodiscard]] constexpr CanonicalOutput evalCanonical(const CanonicalInput &x) const noexcept;
    [[nodiscard]] constexpr CanonicalInput mapToDomain(const CanonicalInput &x) const noexcept;
    [[nodiscard]] constexpr CanonicalInput mapFromDomain(const CanonicalInput &x) const noexcept;
    constexpr void computeScaling(const InputType &a, const InputType &b) noexcept;

    template<std::size_t Rank, class F>
    static constexpr void forEachIndex(const std::array<int, Rank> &ext, F &&body);
#endif
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<std::size_t NCOEFFS, class Func, class... Tags>
[[nodiscard]] PF_CXX20_CONSTEXPR auto fit(Func F, fitInput_t<Func> a, fitInput_t<Func> b, Tags...);

template<class Func, class Spec, class... Tags>
[[nodiscard]] PF_CXX20_CONSTEXPR auto fit(Func F, Spec spec, fitInput_t<Func> a, fitInput_t<Func> b, Tags...);
#endif

#if __cplusplus >= 202002L
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<std::size_t NCOEFFS, auto a, auto b, class Func> [[nodiscard]] constexpr auto fit(Func F) {
    return FuncEvalND<Func, NCOEFFS>(F, a, b);
}

#if PF_HAS_CONSTEXPR_EPS_OVERLOAD
template<double EPS, auto a, auto b, std::size_t MAX_NCOEFFS = 32, std::size_t EVAL_POINTS = 100,
         std::size_t ITERS = 1, class Func>
[[nodiscard]] constexpr auto fit(Func F);
#endif
#endif
#endif

template<typename... EvalTypes>
[[nodiscard]] PF_CXX20_CONSTEXPR FuncEvalMany<EvalTypes...> pack(EvalTypes... evals) noexcept;

} // namespace poly_eval

#include "internal/fast_eval_impl.hpp"
#include "internal/macros_undef.h"
