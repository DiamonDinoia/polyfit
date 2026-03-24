// Correctness tests for the portable trig example.
#include <array>
#include <cmath>
#include <cstdint>
#include <gtest/gtest.h>
#include <limits>
#include <type_traits>
#include <vector>

#include "../examples/portable_trig.hpp"

namespace {

template<typename T>
using BitPattern = std::conditional_t<std::is_same<T, float>::value, std::uint32_t, std::uint64_t>;

template<typename T> BitPattern<T> to_bits(T value) noexcept {
    return polyfit_examples::portable_trig::detail::bitwise_cast_runtime<BitPattern<T>>(value);
}

template<typename T> bool same_bits(T lhs, T rhs) noexcept {
    return to_bits(lhs) == to_bits(rhs);
}

template<typename T> using Pair = std::pair<T, T>;

template<typename T,
         int TolDigits = polyfit_examples::portable_trig::default_approx_digits<T>::value>
Pair<T> eval_angle(T angle) noexcept {
    return polyfit_examples::portable_trig::sincos<T, TolDigits>(angle);
}

template<typename T> Pair<T> remap(unsigned quadrant, T s, T c) noexcept {
    return polyfit_examples::portable_trig::detail::remap_quadrant(quadrant, s, c);
}

template<typename T> constexpr T sincos_tolerance() noexcept {
    if constexpr (std::is_same<T, float>::value) {
        return static_cast<T>(2.0e-7);
    } else {
        return static_cast<T>(3.0e-13);
    }
}

template<typename T> std::vector<T> boundary_angles() {
    std::vector<T> angles = {
        T(0),
        T(-0.0),
        polyfit_examples::portable_trig::detail::pi<T>() / T(4),
        -polyfit_examples::portable_trig::detail::pi<T>() / T(4),
        polyfit_examples::portable_trig::detail::pi<T>(),
        -polyfit_examples::portable_trig::detail::pi<T>(),
        polyfit_examples::portable_trig::detail::pi<T>() * T(256),
        -polyfit_examples::portable_trig::detail::pi<T>() * T(256),
    };

    const T pio2 = polyfit_examples::portable_trig::detail::pi_over_2<T>();
    const T inf = std::numeric_limits<T>::infinity();
    for (int k = -64; k <= 64; ++k) {
        const T base = static_cast<T>(k) * pio2;
        angles.push_back(base);
        angles.push_back(std::nextafter(base, -inf));
        angles.push_back(std::nextafter(base, inf));
    }

    return angles;
}

template<typename T> class PortableTrigTyped : public testing::Test {};
using FloatingTypes = testing::Types<float, double>;
TYPED_TEST_SUITE(PortableTrigTyped, FloatingTypes);

TYPED_TEST(PortableTrigTyped, QuadrantMaskMatchesModuloFourForNegativeAndPositiveInputs) {
    for (long long qi = -19; qi <= 19; ++qi) {
        int expected = static_cast<int>(qi % 4LL);
        if (expected < 0) expected += 4;
        EXPECT_EQ(polyfit_examples::portable_trig::detail::quadrant_from_rounded_multiple(qi),
                  static_cast<unsigned>(expected));
    }
}

TYPED_TEST(PortableTrigTyped, RemapQuadrantMatchesReferenceOutputs) {
    using T = TypeParam;
    constexpr std::array<Pair<T>, 4> expected = {
        Pair<T>{T(2), T(3)},
        Pair<T>{T(3), T(-2)},
        Pair<T>{T(-2), T(-3)},
        Pair<T>{T(-3), T(2)},
    };

    for (unsigned quadrant = 0; quadrant < expected.size(); ++quadrant) {
        EXPECT_EQ(remap<T>(quadrant, T(2), T(3)), expected[quadrant]);
    }
}

TYPED_TEST(PortableTrigTyped, RemapQuadrantHandlesZeroSignBitsCorrectly) {
    using T = TypeParam;
    const T neg_zero = T(-0.0);
    const T pos_zero = T(0.0);

    const Pair<T> q0 = remap<T>(0U, neg_zero, T(1));
    const Pair<T> q1 = remap<T>(1U, neg_zero, T(1));
    const Pair<T> q2 = remap<T>(2U, neg_zero, T(1));
    const Pair<T> q3 = remap<T>(3U, neg_zero, T(1));

    EXPECT_TRUE(same_bits(q0.first, neg_zero));
    EXPECT_TRUE(same_bits(q0.second, T(1)));
    EXPECT_TRUE(same_bits(q1.first, T(1)));
    EXPECT_TRUE(same_bits(q1.second, pos_zero));
    EXPECT_TRUE(same_bits(q2.first, pos_zero));
    EXPECT_TRUE(same_bits(q2.second, T(-1)));
    EXPECT_TRUE(same_bits(q3.first, T(-1)));
    EXPECT_TRUE(same_bits(q3.second, neg_zero));
}

TYPED_TEST(PortableTrigTyped, SincosMatchesStdOnBoundaryAndLargeAngles) {
    using T = TypeParam;
    const T tol = sincos_tolerance<T>();

    for (T angle : boundary_angles<T>()) {
        const Pair<T> approx = eval_angle(angle);
        EXPECT_NEAR(approx.first, std::sin(angle), tol);
        EXPECT_NEAR(approx.second, std::cos(angle), tol);
    }
}

TYPED_TEST(PortableTrigTyped, SincosMatchesStdOnDenseSweep) {
    using T = TypeParam;
    const T pi = polyfit_examples::portable_trig::detail::pi<T>();
    const T tol = sincos_tolerance<T>();

    for (int i = 0; i <= 10000; ++i) {
        const T alpha = static_cast<T>(i) / static_cast<T>(10000);
        const T angle = (alpha * T(2) - T(1)) * T(256) * pi;
        const Pair<T> approx = eval_angle(angle);
        EXPECT_NEAR(approx.first, std::sin(angle), tol);
        EXPECT_NEAR(approx.second, std::cos(angle), tol);
    }
}

TEST(PortableTrigSpecialValues, NonFiniteInputsReturnNaNPairs) {
    using polyfit_examples::portable_trig::sincos;

    const auto inf_pair = sincos<double>(std::numeric_limits<double>::infinity());
    const auto nan_pair = sincos<double>(std::numeric_limits<double>::quiet_NaN());

    EXPECT_TRUE(std::isnan(inf_pair.first));
    EXPECT_TRUE(std::isnan(inf_pair.second));
    EXPECT_TRUE(std::isnan(nan_pair.first));
    EXPECT_TRUE(std::isnan(nan_pair.second));
}

TEST(PortableTrigApi, ExplicitDigitsTemplateArgumentCompilesInSecondSlot) {
    const auto default_digits = polyfit_examples::portable_trig::sincos<double>(0.5);
    const auto explicit_digits = polyfit_examples::portable_trig::sincos<double, 14>(0.5);

    EXPECT_TRUE(std::isfinite(default_digits.first));
    EXPECT_TRUE(std::isfinite(default_digits.second));
    EXPECT_TRUE(std::isfinite(explicit_digits.first));
    EXPECT_TRUE(std::isfinite(explicit_digits.second));
}

TEST(PortableTrigAccuracy, DefaultHighDigitPathMatchesStdOnPrimaryRangeForDouble) {
    const double pi = polyfit_examples::portable_trig::detail::pi<double>();

    for (int i = 0; i <= 4000; ++i) {
        const double alpha = static_cast<double>(i) / 4000.0;
        const double angle = (alpha * 2.0 - 1.0) * pi;
        const Pair<double> approx = eval_angle<double, 15>(angle);
        EXPECT_NEAR(approx.first, std::sin(angle), 5.0e-15);
        EXPECT_NEAR(approx.second, std::cos(angle), 5.0e-15);
    }
}

TEST(PortableTrigAccuracy, HighDigitPathHandlesVeryWideDoubleRanges) {
    const double pi = polyfit_examples::portable_trig::detail::pi<double>();
    constexpr double wide_scale = 1000000.123;

    for (int i = 0; i <= 4000; ++i) {
        const double alpha = static_cast<double>(i) / 4000.0;
        const double angle = (alpha * 2.0 - 1.0) * wide_scale * pi;
        const Pair<double> approx = eval_angle<double, 15>(angle);
        EXPECT_NEAR(approx.first, std::sin(angle), 5.0e-13);
        EXPECT_NEAR(approx.second, std::cos(angle), 5.0e-13);
    }
}

TEST(PortableTrigAccuracy, HighDigitPathHandlesMultiMillionPiRanges) {
    const double pi = polyfit_examples::portable_trig::detail::pi<double>();
    constexpr double wide_scale = 4000000.0;

    for (int i = 0; i <= 4000; ++i) {
        const double alpha = static_cast<double>(i) / 4000.0;
        const double angle = (alpha * 2.0 - 1.0) * wide_scale * pi;
        const Pair<double> approx = eval_angle<double, 15>(angle);
        EXPECT_NEAR(approx.first, std::sin(angle), 5.0e-13);
        EXPECT_NEAR(approx.second, std::cos(angle), 5.0e-13);
    }
}

} // namespace
