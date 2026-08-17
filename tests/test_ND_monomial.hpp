#pragma once

// One instantiation per translation unit: each monomial fit is a heavyweight
// template instantiation (peak RSS of several GiB with sanitizers enabled, so
// the DegrCart product TESTs live one per .cpp).

#include <random>

#include <gtest/gtest.h>

#include "test_ND_shared.hpp"

template<std::size_t IN_DIM, std::size_t OUT_DIM, std::size_t NCOEFFS> void runMonomialTest() {
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    using Input = std::array<double, IN_DIM>;
    using Output = std::array<double, OUT_DIM>;
    const double tol = std::pow(10.0, -static_cast<double>(NCOEFFS - 3));
    Input a{};
    a.fill(-1.0);
    Input b{};
    b.fill(1.0);

    auto approx = poly_eval::fit<NCOEFFS>(sumCos<Input, Output>, a, b);

    for (int i = 0; i < kNumPoints; ++i) {
        Input x;
        for (auto &xi : x) xi = dist(gen);
        auto expected = sumCos<Input, Output>(x);
        auto actual = approx(x);
        for (std::size_t j = 0; j < OUT_DIM; ++j) ASSERT_NEAR(actual[j], expected[j], tol);
    }
}
