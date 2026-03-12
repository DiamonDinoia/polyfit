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

template<typename T, typename = void> struct ValueType {
    using type = remove_cvref_t<T>;
};

template<typename T> struct ValueType<T, std::void_t<typename remove_cvref_t<T>::value_type>> {
    using type = typename remove_cvref_t<T>::value_type;
};

template<typename T> using value_type_or_t = typename ValueType<T>::type;

} // namespace detail

template<class Func> inline constexpr bool takesTupleInput_v = detail::hasTupleSize_v<fitInput_t<Func>>;

template<class Func, std::size_t NCOEFFS, std::size_t ITERS = 1, FusionMode FUSION = FusionMode::Auto>
using FitEvaluator =
    std::conditional_t<takesTupleInput_v<Func>, FuncEvalND<Func, NCOEFFS>, FuncEval<Func, NCOEFFS, ITERS, FUSION>>;

} // namespace poly_eval
