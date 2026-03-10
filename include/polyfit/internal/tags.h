#pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace poly_eval {

// ---------------------------------------------------------------------------
// Fusion control
// ---------------------------------------------------------------------------
enum class FusionMode : std::uint8_t { auto_, always, never };

struct fuse_auto   : std::integral_constant<FusionMode, FusionMode::auto_>  {};
struct fuse_always : std::integral_constant<FusionMode, FusionMode::always> {};
struct fuse_never  : std::integral_constant<FusionMode, FusionMode::never>  {};

// ---------------------------------------------------------------------------
// Value-carrying tags
// ---------------------------------------------------------------------------
template<std::size_t N> struct iters {};
template<std::size_t N> struct maxNCoeffs {};
template<std::size_t N> struct evalPts {};

namespace detail {

// ---------------------------------------------------------------------------
// Tag detection
// ---------------------------------------------------------------------------
template<class T> struct is_polyfit_tag : std::false_type {};
template<std::size_t N> struct is_polyfit_tag<iters<N>> : std::true_type {};
template<std::size_t N> struct is_polyfit_tag<maxNCoeffs<N>> : std::true_type {};
template<std::size_t N> struct is_polyfit_tag<evalPts<N>> : std::true_type {};
template<> struct is_polyfit_tag<fuse_auto> : std::true_type {};
template<> struct is_polyfit_tag<fuse_always> : std::true_type {};
template<> struct is_polyfit_tag<fuse_never> : std::true_type {};

template<class... Tags>
inline constexpr bool all_tags_v =
    (true && ... && is_polyfit_tag<std::remove_cv_t<std::remove_reference_t<Tags>>>::value);

// ---------------------------------------------------------------------------
// Value extraction from tag pack
// ---------------------------------------------------------------------------
template<template<std::size_t> class Tag, class T> struct tag_match {
    static constexpr bool matched = false;
    static constexpr std::size_t val = 0;
};
template<template<std::size_t> class Tag, std::size_t V> struct tag_match<Tag, Tag<V>> {
    static constexpr bool matched = true;
    static constexpr std::size_t val = V;
};

template<template<std::size_t> class Tag, std::size_t Default>
constexpr std::size_t extract_tag() {
    return Default;
}
template<template<std::size_t> class Tag, std::size_t Default, class First, class... Rest>
constexpr std::size_t extract_tag() {
    using F = std::remove_cv_t<std::remove_reference_t<First>>;
    if constexpr (tag_match<Tag, F>::matched)
        return tag_match<Tag, F>::val;
    else
        return extract_tag<Tag, Default, Rest...>();
}

// ---------------------------------------------------------------------------
// Fusion mode extraction
// ---------------------------------------------------------------------------
template<class... Ts> struct fusion_mode_impl {
    static constexpr FusionMode value = FusionMode::auto_;
};
template<class First, class... Rest> struct fusion_mode_impl<First, Rest...> {
    using F = std::remove_cv_t<std::remove_reference_t<First>>;
    static constexpr FusionMode value = std::is_same_v<F, fuse_always> ? FusionMode::always
                                        : std::is_same_v<F, fuse_never> ? FusionMode::never
                                                                        : fusion_mode_impl<Rest...>::value;
};

template<class... Tags>
inline constexpr FusionMode extract_fusion_v = fusion_mode_impl<Tags...>::value;

} // namespace detail
} // namespace poly_eval
