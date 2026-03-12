#pragma once

#include <utility>

#if defined(__GNUC__) || defined(__clang__)
#define PF_ALWAYS_INLINE __attribute__((always_inline)) inline
#elif defined(_MSC_VER)
#define PF_ALWAYS_INLINE __forceinline
#else
#define PF_ALWAYS_INLINE inline
#endif

#if defined(__GNUC__) || defined(__clang__)
#define PF_NO_INLINE __attribute__((noinline))
#elif defined(_MSC_VER)
#define PF_NO_INLINE __declspec(noinline)
#else
#define PF_NO_INLINE
#endif

#if defined(__GNUC__) || defined(__clang__)
#define PF_RESTRICT __restrict__
#elif defined(_MSC_VER)
#define PF_RESTRICT __restrict
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)
#define PF_RESTRICT restrict
#else
#define PF_RESTRICT
#endif

#if PF_HAS_ATTRIBUTE_LIKELY
#define PF_UNLIKELY [[unlikely]]
#define PF_LIKELY [[likely]]
#else
#define PF_UNLIKELY
#define PF_LIKELY
#endif

#if defined(_MSC_VER)
#define PF_ASSUME(cond) __assume(cond)
#elif PF_HAS_BUILTIN(__builtin_assume)
#define PF_ASSUME(cond) __builtin_assume(cond)
#elif defined(__GNUC__) || defined(__GNUG__)
#define PF_ASSUME(cond)                                                                                                \
    do {                                                                                                               \
        if (!(cond)) __builtin_unreachable();                                                                          \
    } while (0)
#else
#define PF_ASSUME(cond) ((void)0)
#endif

#if defined(_MSC_VER)
#define PF_NO_UNIQUE_ADDRESS [[msvc::no_unique_address]]
#elif PF_HAS_ATTRIBUTE_NO_UNIQUE_ADDRESS
#define PF_NO_UNIQUE_ADDRESS [[no_unique_address]]
#else
#define PF_NO_UNIQUE_ADDRESS
#endif

#if defined(_MSC_VER)
#define PF_UNREACHABLE() __assume(0)
#elif PF_HAS_STD_UNREACHABLE
#define PF_UNREACHABLE() std::unreachable()
#elif PF_HAS_BUILTIN_UNREACHABLE
#define PF_UNREACHABLE() __builtin_unreachable()
#else
#define PF_UNREACHABLE() ((void)0)
#endif

// Per-function optimization pragmas for hot evaluation paths.
#if defined(__clang__)
#define PF_FAST_EVAL_BEGIN
#define PF_FAST_EVAL_END
#elif defined(__GNUC__)
#ifdef PF_DISABLE_FAST_EVAL
#define PF_FAST_EVAL_PUSH
#define PF_FAST_EVAL_OPTIMIZE
#define PF_FAST_EVAL_EXTRA
#define PF_FAST_EVAL_POP
#define PF_FAST_EVAL_BEGIN
#define PF_FAST_EVAL_END
#else
#define PF_FAST_EVAL_PUSH _Pragma("GCC push_options")
#define PF_FAST_EVAL_OPTIMIZE                                                                                            \
    _Pragma("GCC optimize(\"tree-vectorize,fp-contract=fast\")")                                                         \
    _Pragma("GCC optimize(\"-fira-hoist-pressure\")")                                                                     \
    _Pragma("GCC optimize(\"-fno-ira-share-spill-slots\")")                                                              \
    _Pragma("GCC optimize(\"-frename-registers\")")
#ifdef PF_EXTRA_FAST_EVAL
#define PF_FAST_EVAL_EXTRA _Pragma(PF_EXTRA_FAST_EVAL)
#else
#define PF_FAST_EVAL_EXTRA
#endif
#define PF_FAST_EVAL_POP _Pragma("GCC pop_options")
#define PF_FAST_EVAL_BEGIN PF_FAST_EVAL_PUSH PF_FAST_EVAL_OPTIMIZE PF_FAST_EVAL_EXTRA
#define PF_FAST_EVAL_END PF_FAST_EVAL_POP
#endif
#elif defined(_MSC_VER)
#define PF_FAST_EVAL_BEGIN __pragma(optimize("gt", on))
#define PF_FAST_EVAL_END __pragma(optimize("", on))
#else
#define PF_FAST_EVAL_BEGIN
#define PF_FAST_EVAL_END
#endif
