#pragma once

#include <array>
#include <cmath>
#include <vector>

#if defined(__has_include)
#  if __has_include(<span>)
#    include <span>
#  endif
#endif

#include "polyfit/internal/feature_macros.h"
#include "polyfit/polyfit.hpp"

constexpr int kNumPoints = 1000;

#if PF_HAS_CONSTEXPR_CMATH
#define PF_TEST_ND_CONSTEXPR constexpr
#else
#define PF_TEST_ND_CONSTEXPR inline
#endif

template<typename Array, typename Output = Array> constexpr auto sumCos(const Array &x) {
    double s = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) s += std::cos(x[i]);
    Output y{};
    y.fill(s);
    for (std::size_t i = 1; i < y.size(); ++i) y[i] += y[i - 1];
    return y;
}

template<class T, std::size_t N> struct FixedVec {
    using value_type = T;
    std::array<T, N> values{};

    constexpr std::size_t size() const noexcept { return N; }
    constexpr void fill(const T &value) noexcept { values.fill(value); }
    constexpr T *data() noexcept { return values.data(); }
    constexpr const T *data() const noexcept { return values.data(); }
    constexpr T &operator[](std::size_t i) noexcept { return values[i]; }
    constexpr const T &operator[](std::size_t i) const noexcept { return values[i]; }
    constexpr auto begin() noexcept { return values.begin(); }
    constexpr auto begin() const noexcept { return values.begin(); }
    constexpr auto end() noexcept { return values.end(); }
    constexpr auto end() const noexcept { return values.end(); }
};

using Array2 = std::array<double, 2>;
using Fixed2 = FixedVec<double, 2>;

PF_TEST_ND_CONSTEXPR auto evalArray2(const Array2 &p) {
    return Array2{std::cos(p[0]) + std::sin(p[1]), p[0] * p[1]};
}

PF_TEST_ND_CONSTEXPR auto evalFixed2(const Fixed2 &p) {
    Fixed2 y{};
    y[0] = std::cos(p[0]) + std::sin(p[1]);
    y[1] = p[0] * p[1];
    return y;
}

inline auto makeArrayApprox2() { return poly_eval::fit(evalArray2, 10, Array2{-1.0, -1.0}, Array2{1.0, 1.0}); }

inline auto makeFixedApprox2() { return poly_eval::fit(evalFixed2, 10, Fixed2{{-1.0, -1.0}}, Fixed2{{1.0, 1.0}}); }

template<class Vec> inline auto samplePoints2() {
    return std::vector<Vec>{{{0.25, -0.5}}, {{-0.2, 0.75}}, {{0.8, -0.1}}};
}

#undef PF_TEST_ND_CONSTEXPR
