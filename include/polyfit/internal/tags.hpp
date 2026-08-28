#pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace poly_eval {

enum class FusionMode : std::uint8_t { Auto, Always, Never };

using FuseAuto = std::integral_constant<FusionMode, FusionMode::Auto>;
using FuseAlways = std::integral_constant<FusionMode, FusionMode::Always>;
using FuseNever = std::integral_constant<FusionMode, FusionMode::Never>;

/// Single-point scalar evaluation kernel, selected at the `FuncEval` /
/// `FuncEvalND` type level.
///
/// - `Hybrid` (default): mixed Estrin/Horner. Best when the single-point
///   `operator()` of the evaluator is the ILP-limited hot path; the evaluator
///   provides its own parallelism.
/// - `Horner`          : the classical Horner chain. Lowest latency and
///   smallest code; best when the caller already has ILP, for example an outer
///   unrolled loop over independent evaluators.
enum class ScalarKernel : std::uint8_t {
    Hybrid,
    Horner,
};

/// Policy hint for `optimal_block_size<>` picking the Hybrid block size.
///
/// - `Latency`   : minimise the dependency chain through one evaluation. The
///   right call when a single point is the bottleneck and the caller cannot
///   expose more ILP at a higher level.
/// - `Throughput`: maximise per-call ILP when the caller already saturates wide
///   SIMD lanes with independent points. K degenerates to one long chain so the
///   FMA pipes stay full across lanes.
/// - `Balanced`  : compromise between the two, clamped by register pressure.
enum class EvalPolicy : std::uint8_t { Latency, Throughput, Balanced };

template<std::size_t N> struct Iters : std::integral_constant<std::size_t, N> {};
template<std::size_t N> struct MaxCoeffs : std::integral_constant<std::size_t, N> {};
template<std::size_t N> struct EvalPts : std::integral_constant<std::size_t, N> {};
/// Compile-time opt-in override for the Hybrid block size K. The default
/// `HybridK<0>` keeps the heuristic of `hybrid_block_size()`; any other value
/// pins K for every Hybrid kernel instantiation downstream.
template<std::size_t N> struct HybridK : std::integral_constant<std::size_t, N> {};

namespace detail {

struct CompileTimeCountTag {};
struct RuntimeCountTag {};
template<class T> using tag_type_t = std::remove_cv_t<std::remove_reference_t<T>>;

template<class T> inline constexpr bool isTag_v = false;
template<std::size_t N> inline constexpr bool isTag_v<Iters<N>> = true;
template<std::size_t N> inline constexpr bool isTag_v<MaxCoeffs<N>> = true;
template<std::size_t N> inline constexpr bool isTag_v<EvalPts<N>> = true;
template<std::size_t N> inline constexpr bool isTag_v<HybridK<N>> = true;
template<> inline constexpr bool isTag_v<FuseAuto> = true;
template<> inline constexpr bool isTag_v<FuseAlways> = true;
template<> inline constexpr bool isTag_v<FuseNever> = true;

template<class... Tags> inline constexpr bool allTags = (... && isTag_v<tag_type_t<Tags>>);

template<class T> inline constexpr bool isFusionTag_v = false;
template<> inline constexpr bool isFusionTag_v<FuseAuto> = true;
template<> inline constexpr bool isFusionTag_v<FuseAlways> = true;
template<> inline constexpr bool isFusionTag_v<FuseNever> = true;
template<class T> inline constexpr FusionMode fusionModeFor = FusionMode::Auto;
template<> inline constexpr FusionMode fusionModeFor<FuseAuto> = FusionMode::Auto;
template<> inline constexpr FusionMode fusionModeFor<FuseAlways> = FusionMode::Always;
template<> inline constexpr FusionMode fusionModeFor<FuseNever> = FusionMode::Never;

template<template<std::size_t> class Tag, class T> inline constexpr bool isCountTag_v = false;
template<template<std::size_t> class Tag, std::size_t N> inline constexpr bool isCountTag_v<Tag, Tag<N>> = true;

template<template<std::size_t> class Tag, class T> inline constexpr std::size_t tagArg_v = 0;
template<template<std::size_t> class Tag, std::size_t N> inline constexpr std::size_t tagArg_v<Tag, Tag<N>> = N;

template<template<std::size_t> class Tag, std::size_t Default> constexpr std::size_t tagValueOf() {
    return Default;
}

template<template<std::size_t> class Tag, std::size_t Default, class T, class... Tags> constexpr std::size_t tagValueOf() {
    using U = tag_type_t<T>;
    if constexpr (isCountTag_v<Tag, U>) {
        return tagArg_v<Tag, U>;
    } else {
        return tagValueOf<Tag, Default, Tags...>();
    }
}

template<template<std::size_t> class Tag, std::size_t Default, class... Tags>
inline constexpr std::size_t tagValue = tagValueOf<Tag, Default, Tags...>();

template<template<std::size_t> class Tag, class... Tags>
inline constexpr std::size_t tagCount = (std::size_t{0} + ... + std::size_t(isCountTag_v<Tag, tag_type_t<Tags>>));

template<class... Tags> struct FusionModeValue : std::integral_constant<FusionMode, FusionMode::Auto> {};

template<class T, class... Tags>
struct FusionModeValue<T, Tags...>
    : std::conditional_t<isFusionTag_v<tag_type_t<T>>, std::integral_constant<FusionMode, fusionModeFor<tag_type_t<T>>>,
                         FusionModeValue<Tags...>> {};

template<class... Tags> inline constexpr FusionMode fusionModeValue = FusionModeValue<Tags...>::value;

template<class... Tags>
inline constexpr std::size_t fusionTagCount =
    (std::size_t{0} + ... + std::size_t(isFusionTag_v<tag_type_t<Tags>>));

template<class... Tags> struct FitOptions {
    static constexpr bool UNIQUE = tagCount<Iters, Tags...> <= 1 && tagCount<MaxCoeffs, Tags...> <= 1
        && tagCount<EvalPts, Tags...> <= 1 && tagCount<HybridK, Tags...> <= 1 && fusionTagCount<Tags...> <= 1;
    static constexpr bool VALID = allTags<Tags...> && UNIQUE;
    static constexpr std::size_t ITERS = tagValue<Iters, 1, Tags...>;
    static constexpr std::size_t MAX_NCOEFFS = tagValue<MaxCoeffs, 32, Tags...>;
    static constexpr std::size_t EVAL_POINTS = tagValue<EvalPts, 100, Tags...>;
    /// 0 means "use the heuristic"; any other value pins K in the Hybrid kernel.
    static constexpr std::size_t HYBRID_K = tagValue<HybridK, 0, Tags...>;
    static constexpr FusionMode FUSION_MODE = fusionModeValue<Tags...>;
};

} // namespace detail
} // namespace poly_eval
