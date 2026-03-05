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

// Define a c++20 constexpr macro
#if __cplusplus >= 202002L
#define PF_C20CONSTEXPR constexpr
#else
#define PF_C20CONSTEXPR
#endif

// C++20 branch-prediction attributes — silently dropped in C++17
#if __cplusplus >= 202002L
#define PF_UNLIKELY [[unlikely]]
#define PF_LIKELY [[likely]]
#else
#define PF_UNLIKELY
#define PF_LIKELY
#endif

#if __cplusplus >= 202302L
#define PF_C23CONSTEXPR constexpr
#define PF_C23CONSTEVAL consteval
#define PF_C23STATIC static
#else
#define PF_C23CONSTEXPR
#define PF_C23CONSTEVAL constexpr
#define PF_C23STATIC
#endif

#if defined(_MSC_VER)
#define PF_ASSUME(cond) __assume(cond)
#elif defined(__clang__)
#define PF_ASSUME(cond) __builtin_assume(cond)
#elif defined(__GNUC__) || defined(__GNUG__)
#define PF_ASSUME(cond)                                                                                                \
    do {                                                                                                               \
        if (!(cond)) __builtin_unreachable();                                                                          \
    } while (0)
#else
#define PF_ASSUME(cond) ((void)0)
#endif

// The eps-based constexpr make_func_eval uses a function parameter in an
// if-constexpr condition — a GCC extension that Clang strictly rejects.
#if __cplusplus >= 202002L && defined(__GNUC__) && !defined(__clang__)
#define PF_HAS_CONSTEXPR_EPS_OVERLOAD 1
#else
#define PF_HAS_CONSTEXPR_EPS_OVERLOAD 0
#endif

// PF_IS_CONSTANT_EVALUATED(): portable wrapper for std::is_constant_evaluated().
// In C++20+ uses the standard API; in C++17 falls back to the compiler builtin
// (supported by GCC ≥9 and Clang ≥9 as an extension).
#if __cplusplus >= 202002L
#define PF_IS_CONSTANT_EVALUATED() std::is_constant_evaluated()
#elif __has_builtin(__builtin_is_constant_evaluated)
#define PF_IS_CONSTANT_EVALUATED() __builtin_is_constant_evaluated()
#else
#define PF_IS_CONSTANT_EVALUATED() false
#endif

// C++26: std::fma, std::abs, std::log10, std::sqrt become constexpr.
#if __cplusplus >= 202602L
#define PF_C26CONSTEXPR constexpr
#else
#define PF_C26CONSTEXPR
#endif

// C++23 native `if consteval` — cleaner than std::is_constant_evaluated() and
// works correctly inside consteval functions (std::is_constant_evaluated() always
// returns false inside a consteval function).
#if __cplusplus >= 202302L
#define PF_IF_CONSTEVAL     if consteval
#define PF_IF_NOT_CONSTEVAL if !consteval
#else
#define PF_IF_CONSTEVAL     if (PF_IS_CONSTANT_EVALUATED())
#define PF_IF_NOT_CONSTEVAL if (!PF_IS_CONSTANT_EVALUATED())
#endif

// Per-function "fast eval" optimization pragmas.
// These push per-function optimize pragmas that layer on top of any global
// compiler flags (so they effectively "push" extra optimizations for hot paths).
// The macros are intentionally conservative: they enable loop-unrolling,
// vectorization and fp-contract by default and only enable aggressive
// fast-math if PF_ALLOW_FAST_MATH is defined by the consumer.
//
// Usage:
//   PF_FAST_EVAL_BEGIN
//   // function definition
//   PF_FAST_EVAL_END
//
#if defined(__clang__)
// Clang does not support GCC push/pop_options pragmas.
// Per-function optimization pragmas are left empty; rely on global flags.
#define PF_FAST_EVAL_BEGIN
#define PF_FAST_EVAL_END
#elif defined(__GNUC__)
// Allow users to disable per-function pragmas entirely via PF_DISABLE_FAST_EVAL.
#ifdef PF_DISABLE_FAST_EVAL
#define PF_FAST_EVAL_PUSH
#define PF_FAST_EVAL_OPTIMIZE
#define PF_FAST_EVAL_EXTRA
#define PF_FAST_EVAL_POP
#define PF_FAST_EVAL_BEGIN
#define PF_FAST_EVAL_END
#else
// Push per-function optimize options that layer on top of user flags.
// Keep them conservative: enable unroll/vectorize and allow fp-contract.
#define PF_FAST_EVAL_PUSH _Pragma("GCC push_options")
#define PF_FAST_EVAL_OPTIMIZE _Pragma("GCC optimize (\"unroll-loops,tree-vectorize,fp-contract=fast\")")
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
// MSVC optimize pragma accepts a short option string; "gt" requests favoring
// speed and global optimizations. We keep this minimal and let global flags
// take precedence where needed.
#define PF_FAST_EVAL_BEGIN __pragma(optimize("gt", on))
#define PF_FAST_EVAL_END __pragma(optimize("", on))

#else
#define PF_FAST_EVAL_BEGIN
#define PF_FAST_EVAL_END
#endif
