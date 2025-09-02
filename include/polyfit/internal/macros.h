#pragma once

// --- ALWAYS INLINE MACRO ---
#if defined(__GNUC__) || defined(__clang__)
// GCC and Clang support __attribute__((always_inline))
#define PF_ALWAYS_INLINE __attribute__((always_inline)) inline
#elif defined(_MSC_VER)
// MSVC supports __forceinline
#define PF_ALWAYS_INLINE __forceinline
#else
// Fallback for other compilers: just use inline.
// This is a weaker hint, but the best we can do generically.
#define PF_ALWAYS_INLINE inline
#endif

// --- NO INLINE MACRO ---
#if defined(__GNUC__) || defined(__clang__)
// GCC and Clang support __attribute__((noinline))
#define PF_NO_INLINE __attribute__((noinline))
#elif defined(_MSC_VER)
// MSVC supports __declspec(noinline)
#define PF_NO_INLINE __declspec(noinline)
#else
// Fallback for other compilers: no specific attribute.
// The compiler will decide whether to inline based on its heuristics.
// This essentially means "don't force inlining or disallow it explicitly".
#define PF_NO_INLINE
#endif

// Define RESTRICT based on the compiler
#if defined(__GNUC__) || defined(__clang__)
// GCC and Clang compilers support __restrict__
#define PF_RESTRICT __restrict__
#elif defined(_MSC_VER)
// Microsoft Visual C++ compiler supports __restrict
#define PF_RESTRICT __restrict
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)
// C99 standard introduced 'restrict' keyword
#define PF_RESTRICT restrict
#else
// Fallback for other compilers or older standards
// In this case, RESTRICT expands to nothing, effectively disabling the keyword.
// This ensures compilation but without the optimization benefits of restrict.
#define PF_RESTRICT
#endif

// Define a c++20 cconstexpr macro
#if __cplusplus >= 202002L
#define PF_C20CONSTEXPR constexpr
#else
#define PF_C20CONSTEXPR

#endif

#if __cplusplus >= 202309L
#define PF_C23STATIC static
#else
#define PF_C23STATIC
#endif

#if defined(_MSC_VER)
#define PF_ASSUME(cond) __assume(cond)
#elif defined(__clang__)
#define PF_ASSUME(cond) __builtin_assume(cond)
#elif defined(__GNUC__) || defined(__GNUG__)
#define PF_ASSUME(cond)                                                                                                  \
    do {                                                                                                               \
        if (!(cond))                                                                                                   \
            __builtin_unreachable();                                                                                   \
    } while (0)
#else
#define PF_ASSUME(cond) ((void)0)
#endif

