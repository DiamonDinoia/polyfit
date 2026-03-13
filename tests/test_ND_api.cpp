#include <gtest/gtest.h>

#include "test_ND_shared.hpp"

TEST(Eval, GenericFixedSizeContainerFitAndEval) {
    const auto approx = makeFixedApprox2();
    const Fixed2 x{{0.25, -0.5}};
    const auto y = approx(x);
    EXPECT_NEAR(y[0], std::cos(0.25) + std::sin(-0.5), 1e-8);
    EXPECT_NEAR(y[1], 0.25 * -0.5, 1e-8);
}

TEST(Eval, NdCanonicalConvenienceApisMatchReference) {
    const auto approx = makeArrayApprox2();
    const auto y = approx({0.25, -0.5});
    EXPECT_NEAR(y[0], std::cos(0.25) + std::sin(-0.5), 1e-8);
    EXPECT_NEAR(y[1], 0.25 * -0.5, 1e-8);

    const auto byArray = approx(Array2{0.25, -0.5});
    const auto byArgs = approx(0.25, -0.5);
    EXPECT_NEAR(byArgs[0], byArray[0], 1e-12);
    EXPECT_NEAR(byArgs[1], byArray[1], 1e-12);

    auto pts = samplePoints2<decltype(approx)::CanonicalInput>();
    std::vector<decltype(approx)::CanonicalOutput> out(pts.size());
    approx(pts.data(), out.data(), pts.size());
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const auto expected = approx(pts[i]);
        EXPECT_NEAR(out[i][0], expected[0], 1e-12);
        EXPECT_NEAR(out[i][1], expected[1], 1e-12);
    }

#if defined(__cpp_lib_span) && (__cpp_lib_span >= 202002L)
    approx(std::span<const decltype(approx)::CanonicalInput>(pts), std::span<decltype(approx)::CanonicalOutput>(out));
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const auto expected = approx(pts[i]);
        EXPECT_NEAR(out[i][0], expected[0], 1e-12);
        EXPECT_NEAR(out[i][1], expected[1], 1e-12);
    }
#endif
}

TEST(Eval, NdGenericConvenienceApisMatchReference) {
    const auto approx = makeFixedApprox2();
    auto pts = samplePoints2<Fixed2>();
    std::vector<Fixed2> out(pts.size());

    const auto y = approx(Fixed2{{0.25, -0.5}});
    EXPECT_NEAR(y[0], std::cos(0.25) + std::sin(-0.5), 1e-8);
    EXPECT_NEAR(y[1], 0.25 * -0.5, 1e-8);

    approx(pts, out);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const auto expected = approx(pts[i]);
        EXPECT_NEAR(out[i][0], expected[0], 1e-12);
        EXPECT_NEAR(out[i][1], expected[1], 1e-12);
    }
}

TEST(Eval, NdBraceInitArgumentsAreInferredFromFunction) {
    auto approx = poly_eval::fit([](const std::array<double, 2> &p) {
        return std::array<double, 2>{std::cos(p[0]) + std::sin(p[1]), p[0] * p[1]};
    }, 10, {-1.0, -1.0}, {1.0, 1.0});

    const auto y = approx({0.25, -0.5});
    EXPECT_NEAR(y[0], std::cos(0.25) + std::sin(-0.5), 1e-8);
    EXPECT_NEAR(y[1], 0.25 * -0.5, 1e-8);
}

#if __cplusplus >= 202002L
TEST(Eval, TemplateParameterNdFit) {
    constexpr std::array<double, 2> a{-1.0, -1.0};
    constexpr std::array<double, 2> b{1.0, 1.0};
    const auto approx = poly_eval::fit<8, a, b>([](const std::array<double, 2> &p) {
        return std::array<double, 2>{p[0] + p[1], p[0] * p[1]};
    });

    const auto y = approx({0.25, -0.5});
    EXPECT_NEAR(y[0], -0.25, 1e-8);
    EXPECT_NEAR(y[1], -0.125, 1e-8);
}
#endif
