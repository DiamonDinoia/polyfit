#include <gtest/gtest.h>
#include "polyfit/fast_eval.hpp"

using namespace poly_eval;

TEST(Utils, LinspaceScalar) {
    auto pts = detail::linspace<0>(0.0, 1.0, 5);
    ASSERT_EQ(pts.size(), 5u);
    EXPECT_DOUBLE_EQ(pts[0], 0.0);
    EXPECT_DOUBLE_EQ(pts.back(), 1.0);
    for (size_t i = 1; i < pts.size(); ++i) {
        double step = pts[i] - pts[i - 1];
        EXPECT_NEAR(step, 0.25, 1e-12);
    }
}

TEST(Utils, LinspaceArray) {
    using V = std::array<double, 2>;
    V a{0.0, -1.0};
    V b{1.0, 1.0};
    auto pts = detail::linspace<0>(a, b, 3);
    ASSERT_EQ(pts.size(), 3u);
    EXPECT_NEAR(pts[0][0], 0.0, 1e-12);
    EXPECT_NEAR(pts[0][1], -1.0, 1e-12);
    EXPECT_NEAR(pts[2][0], 1.0, 1e-12);
    EXPECT_NEAR(pts[2][1], 1.0, 1e-12);
    // middle point should be midpoint
    EXPECT_NEAR(pts[1][0], 0.5, 1e-12);
    EXPECT_NEAR(pts[1][1], 0.0, 1e-12);
}
