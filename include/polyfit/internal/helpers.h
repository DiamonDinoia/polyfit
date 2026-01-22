// Helper utilities for internal use (domain mapping, scaling)
#ifndef POLYFIT_INTERNAL_HELPERS_H
#define POLYFIT_INTERNAL_HELPERS_H

#include <array>
#include <type_traits>
#include <cmath>
#include <xsimd/xsimd.hpp>

namespace polyfit::internal::helpers {

// Scalar mapping to canonical domain [-1,1]
template <class ArgT, class ScalarT>
PF_ALWAYS_INLINE constexpr ArgT map_to_domain_scalar(const ArgT arg, const ScalarT low, const ScalarT hi) noexcept {
    if constexpr (std::is_arithmetic_v<ArgT>) {
        return static_cast<ArgT>(ArgT(0.5) * ((arg / static_cast<ArgT>(low)) + static_cast<ArgT>(hi)));
    }
    return static_cast<ArgT>(ArgT(0.5) * ((arg / ArgT(low)) + ArgT(hi)));
}

// Scalar/batch mapping from canonical domain back to original
template <class ArgT, class ScalarT>
PF_ALWAYS_INLINE constexpr ArgT map_from_domain_scalar(const ArgT arg, const ScalarT low, const ScalarT hi) noexcept {
    if constexpr (std::is_arithmetic_v<ArgT>) {
        return static_cast<ArgT>((std::fma(ArgT(2.0), arg, -ArgT(hi))) * ArgT(low));
    }
    return static_cast<ArgT>((xsimd::fms(ArgT(2.0), arg, ArgT(hi))) * ArgT(low));
}

// Array (std::array) mapping overloads
template <class T, std::size_t N>
PF_ALWAYS_INLINE constexpr std::array<T, N>
map_to_domain_array(const std::array<T, N> &t, const std::array<T, N> &low, const std::array<T, N> &hi) noexcept {
    std::array<T, N> out{};
    for (std::size_t d = 0; d < N; ++d) {
        out[d] = map_to_domain_scalar<T, T>(t[d], low[d], hi[d]);
    }
    return out;
}

template <class T, std::size_t N>
PF_ALWAYS_INLINE constexpr std::array<T, N>
map_from_domain_array(const std::array<T, N> &x, const std::array<T, N> &low, const std::array<T, N> &hi) noexcept {
    std::array<T, N> out{};
    for (std::size_t d = 0; d < N; ++d) {
        out[d] = map_from_domain_scalar<T, T>(x[d], low[d], hi[d]);
    }
    return out;
}

template <class T, std::size_t N>
PF_ALWAYS_INLINE constexpr void
compute_scaling_array(const std::array<T, N> &from, const std::array<T, N> &to,
                      std::array<T, N> &low, std::array<T, N> &high) noexcept {
    for (std::size_t dim = 0; dim < N; ++dim) {
        low[dim] = T(1) / (to[dim] - from[dim]);
        high[dim] = (to[dim] + from[dim]);
    }
}

} // namespace polyfit::internal::helpers

#endif // POLYFIT_INTERNAL_HELPERS_H
