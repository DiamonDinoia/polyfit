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

template<class Func, std::size_t, std::size_t, FusionMode = FusionMode::auto_> class FuncEval;
template<typename T> using remove_cvref_t = std::remove_cv_t<std::remove_reference_t<T>>;

template<typename T, std::size_t N_compile_time>
using Buffer = std::conditional_t<N_compile_time == 0, std::vector<T>, std::array<T, N_compile_time>>;

template<typename T, std::size_t N>
constexpr Buffer<T, N> make_buffer(std::size_t runtime_size) {
    Buffer<T, N> buf{};
    if constexpr (N == 0) {
        buf.resize(runtime_size);
    }
    return buf;
}

template<typename T> struct function_traits : function_traits<decltype(&remove_cvref_t<T>::operator())> {};

template<typename R, typename Arg> struct function_traits<R (*)(Arg)> {
    using result_type = R;
    using arg0_type = Arg;
};

template<typename R, typename Arg> struct function_traits<std::function<R(Arg)>> {
    using result_type = R;
    using arg0_type = Arg;
};

template<typename F, typename R, typename Arg> struct function_traits<R (F::*)(Arg) const> {
    using result_type = R;
    using arg0_type = Arg;
};

template<typename F, typename R, typename Arg> struct function_traits<R (F::*)(Arg)> {
    using result_type = R;
    using arg0_type = Arg;
};

template<typename R, typename T> struct function_traits<R (*)(const T &)> {
    using result_type = R;
    using arg0_type = T;
};

template<typename T, typename = void> struct is_tuple_like : std::false_type {};

template<typename T>
struct is_tuple_like<T, std::void_t<decltype(std::tuple_size_v<remove_cvref_t<T>>)>> : std::true_type {};

#if __cpp_concepts >= 201907L
template<typename T>
concept tuple_like = is_tuple_like<T>::value;
#endif

template<typename T, typename = void> struct tuple_size_or_zero : std::integral_constant<std::size_t, 0> {};

template<typename T>
struct tuple_size_or_zero<T, std::void_t<decltype(std::tuple_size_v<remove_cvref_t<T>>)>>
    : std::integral_constant<std::size_t, std::tuple_size_v<remove_cvref_t<T>>> {};

template<typename T> struct is_func_eval : std::false_type {};
template<typename... Args> struct is_func_eval<FuncEval<Args...>> : std::true_type {};

} // namespace poly_eval

template<typename T> struct is_complex : std::false_type {};
template<typename T> struct is_complex<std::complex<T>> : std::true_type {};
template<typename T> inline constexpr bool is_complex_v = is_complex<poly_eval::remove_cvref_t<T>>::value;

template<typename, typename = void> struct has_tuple_size : std::false_type {};
template<typename T> struct has_tuple_size<T, std::void_t<decltype(std::tuple_size<T>::value)>> : std::true_type {};
template<typename T> inline constexpr bool has_tuple_size_v = has_tuple_size<poly_eval::remove_cvref_t<T>>::value;

template<typename T, typename = void> struct value_type_or_identity {
    using type = T;
};

template<typename T> struct value_type_or_identity<T, std::void_t<typename T::value_type>> {
    using type = typename T::value_type;
};
