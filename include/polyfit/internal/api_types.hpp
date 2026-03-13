#pragma once

#include <array>
#include <complex>
#include <functional>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "tags.h"

namespace poly_eval {

template<class Func, std::size_t, std::size_t, FusionMode = FusionMode::Auto> class FuncEval;
template<class Func, std::size_t> class FuncEvalND;
template<typename T> using remove_cvref_t = std::remove_cv_t<std::remove_reference_t<T>>;
template<bool Condition, class T = int> using enable_if_t = std::enable_if_t<Condition, T>;
template<class...> inline constexpr bool alwaysFalse_v = false;
template<class T> inline constexpr bool isIntegralLike_v = std::is_integral_v<remove_cvref_t<T>>;
template<class T> inline constexpr bool isFloatingPointLike_v = std::is_floating_point_v<remove_cvref_t<T>>;

template<typename T, std::size_t NCOMPILE>
using Buffer = std::conditional_t<NCOMPILE == 0, std::vector<T>, std::array<T, NCOMPILE>>;

template<typename T, std::size_t NCOMPILE>
constexpr Buffer<T, NCOMPILE> makeBuffer(std::size_t runtimeSize) {
    Buffer<T, NCOMPILE> buf{};
    if constexpr (NCOMPILE == 0) {
        buf.resize(runtimeSize);
    }
    return buf;
}

template<typename R, typename Arg> struct FunctionSig {
    using result_type = R;
    using arg0_type = remove_cvref_t<Arg>;
};

template<typename T> struct FunctionTraits : FunctionTraits<decltype(&remove_cvref_t<T>::operator())> {};

template<typename R, typename Arg> struct FunctionTraits<R(Arg)> : FunctionSig<R, Arg> {};

template<typename R, typename Arg> struct FunctionTraits<R (*)(Arg)> : FunctionTraits<R(Arg)> {};

template<typename R, typename Arg> struct FunctionTraits<std::function<R(Arg)>> : FunctionTraits<R(Arg)> {};

template<typename F, typename R, typename Arg> struct FunctionTraits<R (F::*)(Arg)> : FunctionTraits<R(Arg)> {};

template<typename F, typename R, typename Arg> struct FunctionTraits<R (F::*)(Arg) const> : FunctionTraits<R(Arg)> {};

template<class Func> using fitInput_t = typename FunctionTraits<Func>::arg0_type;
template<class Func> using fitOutput_t = typename FunctionTraits<Func>::result_type;

namespace detail {

template<typename T> inline constexpr bool isComplexBase_v = false;
template<typename T> inline constexpr bool isComplexBase_v<std::complex<T>> = true;
template<typename T> inline constexpr bool isComplex_v = isComplexBase_v<remove_cvref_t<T>>;

template<typename T, typename = void> inline constexpr bool hasTupleSizeBase_v = false;
template<typename T> inline constexpr bool hasTupleSizeBase_v<T, std::void_t<decltype(std::tuple_size<T>::value)>> = true;
template<typename T> inline constexpr bool hasTupleSize_v = hasTupleSizeBase_v<remove_cvref_t<T>>;

template<typename T, typename = void> inline constexpr bool hasSubscriptBase_v = false;
template<typename T>
inline constexpr bool
    hasSubscriptBase_v<T, std::void_t<decltype(std::declval<T &>()[std::size_t{}]),
                                      decltype(std::declval<const T &>()[std::size_t{}])>> = true;
template<typename T> inline constexpr bool hasSubscript_v = hasSubscriptBase_v<remove_cvref_t<T>>;

template<typename T, typename = void> inline constexpr bool hasSizeMethodBase_v = false;
template<typename T> inline constexpr bool hasSizeMethodBase_v<T, std::void_t<decltype(std::declval<const T &>().size())>> = true;
template<typename T> inline constexpr bool hasSizeMethod_v = hasSizeMethodBase_v<remove_cvref_t<T>>;

template<typename T, typename = void> inline constexpr bool hasDataBase_v = false;
template<typename T> inline constexpr bool hasDataBase_v<T, std::void_t<decltype(std::declval<T &>().data())>> = true;
template<typename T> inline constexpr bool hasData_v = hasDataBase_v<remove_cvref_t<T>>;
template<typename T> inline constexpr bool isDataBackedContainer_v = hasData_v<T> && hasSizeMethod_v<T>;

// Detects constexpr .size() — requires default-constructibility. Types without
// a default constructor fall back to tuple_size or are treated as non-fixed.
template<typename T, typename = void> struct StaticSizeExpr : std::false_type {};
template<typename T>
struct StaticSizeExpr<T, std::void_t<decltype(std::integral_constant<std::size_t, remove_cvref_t<T>{}.size()>{})>>
    : std::true_type {};

template<typename T> constexpr std::size_t fixedSize() {
    using V = remove_cvref_t<T>;
    if constexpr (hasTupleSize_v<V>) {
        return std::tuple_size_v<V>;
    } else if constexpr (StaticSizeExpr<V>::value) {
        return V{}.size();
    } else {
        return 0;
    }
}
template<typename T> inline constexpr std::size_t fixed_size_v = fixedSize<T>();

template<typename T> using index_value_t = remove_cvref_t<decltype(std::declval<const remove_cvref_t<T> &>()[std::size_t{}])>;
template<typename T> using data_value_t = remove_cvref_t<decltype(std::declval<remove_cvref_t<T> &>().data()[std::size_t{}])>;
template<typename T>
inline constexpr bool isFixedIndexable_v =
    (fixed_size_v<T> != 0) && hasSubscript_v<T> && !isComplex_v<T> && !std::is_arithmetic_v<remove_cvref_t<T>>;

template<class T> struct FixedContainerAccess {
    using container_type = remove_cvref_t<T>;
    using value_type = index_value_t<container_type>;

    static constexpr std::size_t size = fixed_size_v<container_type>;
    static constexpr bool has_data = hasData_v<container_type>;

    [[nodiscard]] static constexpr decltype(auto) get(container_type &value, std::size_t i) noexcept { return value[i]; }
    [[nodiscard]] static constexpr decltype(auto) get(const container_type &value, std::size_t i) noexcept {
        return value[i];
    }

    [[nodiscard]] static constexpr auto data(container_type &value) noexcept {
        static_assert(has_data, "Container does not expose data()");
        return value.data();
    }
    [[nodiscard]] static constexpr auto data(const container_type &value) noexcept {
        static_assert(has_data, "Container does not expose data()");
        return value.data();
    }
};

template<class To, class From, enable_if_t<isFixedIndexable_v<To> && isFixedIndexable_v<From>, int> = 0>
[[nodiscard]] constexpr To fixedContainerCast(const From &from) noexcept {
    using ToAccess = FixedContainerAccess<To>;
    using FromAccess = FixedContainerAccess<From>;
    static_assert(ToAccess::size == FromAccess::size, "Fixed-size containers must have matching rank");

    To out{};
    for (std::size_t i = 0; i < ToAccess::size; ++i) {
        ToAccess::get(out, i) = static_cast<typename ToAccess::value_type>(FromAccess::get(from, i));
    }
    return out;
}

template<class Point, std::size_t Dim, class Scalar>
inline constexpr bool isCompatibleNdPoint_v =
    isFixedIndexable_v<Point> && fixed_size_v<Point> == Dim && std::is_convertible_v<index_value_t<Point>, Scalar>;

template<std::size_t Dim, class Scalar, class... Coords>
inline constexpr bool isCompatibleNdCoordPack_v =
    sizeof...(Coords) == Dim && (std::is_convertible_v<Coords, Scalar> && ...);

template<class Points, std::size_t Dim, class Scalar, class = void> struct IsDataBackedNdPointContainer : std::false_type {};
template<class Points, std::size_t Dim, class Scalar>
struct IsDataBackedNdPointContainer<Points, Dim, Scalar, std::void_t<decltype(std::declval<remove_cvref_t<Points> &>().data())>>
    : std::bool_constant<isDataBackedContainer_v<Points> && isCompatibleNdPoint_v<data_value_t<Points>, Dim, Scalar>> {};
template<class Points, std::size_t Dim, class Scalar>
inline constexpr bool isDataBackedNdPointContainer_v = IsDataBackedNdPointContainer<Points, Dim, Scalar>::value;

template<class Outputs, class Output, class = void> struct IsDataBackedNdOutputContainer : std::false_type {};
template<class Outputs, class Output>
struct IsDataBackedNdOutputContainer<Outputs, Output, std::void_t<decltype(std::declval<remove_cvref_t<Outputs> &>().data())>>
    : std::bool_constant<isDataBackedContainer_v<Outputs> && std::is_assignable_v<data_value_t<Outputs> &, Output>> {};
template<class Outputs, class Output>
inline constexpr bool isDataBackedNdOutputContainer_v = IsDataBackedNdOutputContainer<Outputs, Output>::value;

template<class Points, class Outputs, std::size_t Dim, class InputScalar, class Output>
inline constexpr bool isCompatibleNdBatchContainers_v = isDataBackedNdPointContainer_v<Points, Dim, InputScalar> &&
                                                        isDataBackedNdOutputContainer_v<Outputs, Output>;

template<typename T, typename = void> struct ValueType {
    using type = remove_cvref_t<T>;
};

template<typename T> struct ValueType<T, std::void_t<typename remove_cvref_t<T>::value_type>> {
    using type = typename remove_cvref_t<T>::value_type;
};

template<typename T> using value_type_or_t = typename ValueType<T>::type;

} // namespace detail

template<class Func>
inline constexpr bool takesNdInput_v =
    detail::isFixedIndexable_v<fitInput_t<Func>> && detail::isFixedIndexable_v<fitOutput_t<Func>>;

template<class Func, std::size_t NCOEFFS, std::size_t ITERS = 1, FusionMode FUSION = FusionMode::Auto>
using FitEvaluator =
    std::conditional_t<takesNdInput_v<Func>, FuncEvalND<Func, NCOEFFS>, FuncEval<Func, NCOEFFS, ITERS, FUSION>>;

} // namespace poly_eval
