// test_1D_helpers.h — shared helpers for split test_1D_* translation units
#pragma once

#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <random>

#include "polyfit/polyfit.hpp"

static constexpr std::size_t kNumRandomTests = 100;

inline auto double_func = [](const double x) { return std::cos(x); };
inline auto float_func  = [](const float  x) { return std::cos(x); };
inline auto complex_func = [](const double x) {
    return std::complex<double>(x * x, std::sin(x));
};

template <typename T, typename X = double, typename F>
void batch_verify(const F &f, const std::vector<X> &xs, const std::vector<T> &ys,
                  const double tol) {
    for (std::size_t i = 0; i < xs.size(); ++i) {
        EXPECT_LE(poly_eval::detail::relativeL2Norm(ys[i], f(xs[i])), tol)
            << "Failed at x=" << xs[i] << ": expected " << f(xs[i]) << ", got " << ys[i];
    }
}
