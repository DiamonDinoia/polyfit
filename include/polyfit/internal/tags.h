#pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace poly_eval {

enum class FusionMode : std::uint8_t { Auto, Always, Never };

struct FuseAuto : std::integral_constant<FusionMode, FusionMode::Auto> {};
struct FuseAlways : std::integral_constant<FusionMode, FusionMode::Always> {};
struct FuseNever : std::integral_constant<FusionMode, FusionMode::Never> {};

template<std::size_t N> struct Iters {};
template<std::size_t N> struct MaxCoeffs {};
template<std::size_t N> struct EvalPts {};

namespace detail {

struct CompileTimeCountTag {};
struct RuntimeCountTag {};
template<class T> using tag_type_t = std::remove_cv_t<std::remove_reference_t<T>>;

template<class T> struct IsTag : std::false_type {};
template<std::size_t N> struct IsTag<Iters<N>> : std::true_type {};
template<std::size_t N> struct IsTag<MaxCoeffs<N>> : std::true_type {};
template<std::size_t N> struct IsTag<EvalPts<N>> : std::true_type {};
template<> struct IsTag<FuseAuto> : std::true_type {};
template<> struct IsTag<FuseAlways> : std::true_type {};
template<> struct IsTag<FuseNever> : std::true_type {};

template<class... Tags>
inline constexpr bool allTags = (true && ... && IsTag<tag_type_t<Tags>>::value);

template<class T> struct IsFusionTag : std::false_type {};
template<> struct IsFusionTag<FuseAuto> : std::true_type {};
template<> struct IsFusionTag<FuseAlways> : std::true_type {};
template<> struct IsFusionTag<FuseNever> : std::true_type {};

template<template<std::size_t> class Tag, class T> struct TagValue {
    static constexpr bool matched = false;
    static constexpr std::size_t val = 0;
};
template<template<std::size_t> class Tag, std::size_t V> struct TagValue<Tag, Tag<V>> {
    static constexpr bool matched = true;
    static constexpr std::size_t val = V;
};

template<template<std::size_t> class Tag, std::size_t Default, class... Tags>
inline constexpr std::size_t tagValue = [] {
    std::size_t value = Default;
    [[maybe_unused]] bool matched = false;
    ((matched = matched || (TagValue<Tag, tag_type_t<Tags>>::matched
                                ? (value = TagValue<Tag, tag_type_t<Tags>>::val, true)
                                : false)),
     ...);
    return value;
}();

template<template<std::size_t> class Tag, class... Tags>
inline constexpr std::size_t tagCount =
    (std::size_t{0} + ... + std::size_t(TagValue<Tag, tag_type_t<Tags>>::matched));

template<class T>
inline constexpr FusionMode fusionModeFor =
    std::is_same_v<tag_type_t<T>, FuseAlways> ? FusionMode::Always
    : std::is_same_v<tag_type_t<T>, FuseNever> ? FusionMode::Never
                                               : FusionMode::Auto;

template<class... Tags>
inline constexpr FusionMode fusionModeValue = [] {
    FusionMode mode = FusionMode::Auto;
    [[maybe_unused]] bool matched = false;
    ((matched = matched || (fusionModeFor<Tags> != FusionMode::Auto ? (mode = fusionModeFor<Tags>, true) : false)), ...);
    return mode;
}();

template<class... Tags>
inline constexpr std::size_t fusionTagCount =
    (std::size_t{0} + ... + std::size_t(IsFusionTag<tag_type_t<Tags>>::value));

template<class... Tags> struct FitOptions {
    static constexpr bool UNIQUE = tagCount<Iters, Tags...> <= 1 && tagCount<MaxCoeffs, Tags...> <= 1
        && tagCount<EvalPts, Tags...> <= 1 && fusionTagCount<Tags...> <= 1;
    static constexpr bool VALID = allTags<Tags...> && UNIQUE;
    static constexpr std::size_t ITERS = tagValue<Iters, 1, Tags...>;
    static constexpr std::size_t MAX_NCOEFFS = tagValue<MaxCoeffs, 32, Tags...>;
    static constexpr std::size_t EVAL_POINTS = tagValue<EvalPts, 100, Tags...>;
    static constexpr FusionMode FUSION_MODE = fusionModeValue<Tags...>;
};

} // namespace detail
} // namespace poly_eval
