#pragma once

#include <type_traits>

#if defined(__has_include)
#  if __has_include(<version>)
#    include <version>
#  endif
#endif

#if defined(_MSVC_LANG) && (_MSVC_LANG > __cplusplus)
#define PF_CPLUSPLUS _MSVC_LANG
#else
#define PF_CPLUSPLUS __cplusplus
#endif

#if defined(__has_builtin)
#define PF_HAS_BUILTIN(x) __has_builtin(x)
#else
#define PF_HAS_BUILTIN(x) 0
#endif

#if defined(__has_cpp_attribute)
#define PF_HAS_CPP_ATTRIBUTE(x) __has_cpp_attribute(x)
#else
#define PF_HAS_CPP_ATTRIBUTE(x) 0
#endif

#define PF_HAS_CXX20 (PF_CPLUSPLUS >= 202002L)
#define PF_HAS_CXX23 (PF_CPLUSPLUS >= 202302L)
#define PF_HAS_CXX26 (PF_CPLUSPLUS > 202302L)

#if defined(__cpp_if_consteval) && (__cpp_if_consteval >= 202106L)
#define PF_HAS_IF_CONSTEVAL 1
#else
#define PF_HAS_IF_CONSTEVAL 0
#endif

#if defined(__cpp_lib_constexpr_cmath)
#define PF_HAS_CONSTEXPR_CMATH 1
#else
#define PF_HAS_CONSTEXPR_CMATH 0
#endif

#if defined(__cpp_constexpr) && (__cpp_constexpr >= 202211L)
#define PF_HAS_CONSTEXPR_STATIC_LOCAL 1
#else
#define PF_HAS_CONSTEXPR_STATIC_LOCAL 0
#endif

#if (PF_HAS_CPP_ATTRIBUTE(likely) >= 201803L) && (!defined(_MSC_VER) || PF_HAS_CXX20)
#define PF_HAS_ATTRIBUTE_LIKELY 1
#else
#define PF_HAS_ATTRIBUTE_LIKELY 0
#endif

#if PF_HAS_CPP_ATTRIBUTE(no_unique_address) >= 201803L
#define PF_HAS_ATTRIBUTE_NO_UNIQUE_ADDRESS 1
#else
#define PF_HAS_ATTRIBUTE_NO_UNIQUE_ADDRESS 0
#endif

#if defined(_MSC_VER) || PF_HAS_BUILTIN(__builtin_unreachable) || defined(__GNUC__) || defined(__GNUG__)
#define PF_HAS_BUILTIN_UNREACHABLE 1
#else
#define PF_HAS_BUILTIN_UNREACHABLE 0
#endif

#if defined(__cpp_lib_unreachable)
#define PF_HAS_STD_UNREACHABLE 1
#else
#define PF_HAS_STD_UNREACHABLE 0
#endif

#if PF_HAS_CXX20
#define PF_CXX20_CONSTEXPR constexpr
#else
#define PF_CXX20_CONSTEXPR
#endif

#if PF_HAS_CXX20
#define PF_CXX20_CONSTEVAL consteval
#else
#define PF_CXX20_CONSTEVAL constexpr
#endif

#if PF_HAS_CONSTEXPR_STATIC_LOCAL
#define PF_STATIC_CONSTEXPR_LOCAL static constexpr
#else
#define PF_STATIC_CONSTEXPR_LOCAL constexpr
#endif

// The eps-based constexpr fit uses a function parameter in an if-constexpr
// condition, which is accepted on the supported GCC path but rejected by Clang.
#if PF_HAS_CXX20 && defined(__GNUC__) && !defined(__clang__)
#define PF_HAS_CONSTEXPR_EPS_OVERLOAD 1
#else
#define PF_HAS_CONSTEXPR_EPS_OVERLOAD 0
#endif

#if PF_HAS_CXX20
#define PF_IS_CONSTANT_EVALUATED() std::is_constant_evaluated()
#elif defined(_MSC_VER) && _MSC_VER >= 1925
#define PF_IS_CONSTANT_EVALUATED() __builtin_is_constant_evaluated()
#elif PF_HAS_BUILTIN(__builtin_is_constant_evaluated)
#define PF_IS_CONSTANT_EVALUATED() __builtin_is_constant_evaluated()
#else
#define PF_IS_CONSTANT_EVALUATED() false
#endif

#if PF_HAS_IF_CONSTEVAL && !defined(__cppcheck__)
#define PF_IF_CONSTEVAL if consteval
#define PF_IF_NOT_CONSTEVAL if !consteval
#else
#define PF_IF_CONSTEVAL if (PF_IS_CONSTANT_EVALUATED())
#define PF_IF_NOT_CONSTEVAL if (!PF_IS_CONSTANT_EVALUATED())
#endif
