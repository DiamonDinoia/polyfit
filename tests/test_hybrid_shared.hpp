#pragma once

// Shared helpers for the hybrid Estrin/Horner evaluator tests. The suite
// lives one-TYPED_TEST-per-TU because the compile-time `hybrid<D>` and
// `hybrid_transposed<N>` paths each contain several layers of nested
// `poet::static_for`, and instantiating them for D = 1..32 in a single TU both
// trips MSVC's C1202 template-instantiation-context limit and inflates cc1plus
// memory. The fixture and type list are at global scope: gtest macros paste
// the names into identifiers, and the fixture type must be identical across
// the TUs.

#include <array>
#include <gtest/gtest.h>
#include <limits>
#include <random>
#include <type_traits>
#include <vector>

#include "polyfit/polyfit.hpp"

#include <poet/poet.hpp>

using std::size_t;

namespace test_hybrid_detail {

template<typename T> constexpr T eps = std::numeric_limits<T>::epsilon() * T(100);

inline std::mt19937 rng(42);
inline std::uniform_real_distribution<double> uni_dist(-1.0, 1.0);

template<typename T> T uni() { return static_cast<T>(uni_dist(rng)); }

template<typename T> std::vector<T> random_vector(size_t n) {
    std::vector<T> v(n);
    for (auto &x : v) x = uni<T>();
    return v;
}

template<typename Output, typename Input>
constexpr Output naive_horner_scalar(Input x, const Output *c, size_t n) {
    Output acc = c[0];
    for (size_t i = 1; i < n; ++i) acc = acc * x + c[i];
    return acc;
}

} // namespace test_hybrid_detail

template<typename T> class HybridTyped : public testing::Test {};
typedef testing::Types<float, double> FloatingTypes;
TYPED_TEST_SUITE(HybridTyped, FloatingTypes);
