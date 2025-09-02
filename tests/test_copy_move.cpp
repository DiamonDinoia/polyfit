#include <array>
#include <cmath>
#include <gtest/gtest.h>

#include "polyfit/fast_eval.hpp"

using namespace poly_eval;

TEST(CopyMove, FuncEval1D) {
    auto f = [](double x) { return std::sin(x); };
    auto base = make_func_eval<5>(f, -1.0, 1.0);
    using FE = decltype(base);

    FE copy(base);
    EXPECT_NEAR(base(0.3), copy(0.3), 1e-12);

    FE copy_assign = make_func_eval<5>(f, -1.0, 1.0);
    copy_assign = base;
    EXPECT_NEAR(base(0.5), copy_assign(0.5), 1e-12);

    auto temp = make_func_eval<5>(f, -1.0, 1.0);
    auto ref = temp(0.7);
    FE move_ctor(std::move(temp));
    EXPECT_NEAR(ref, move_ctor(0.7), 1e-12);

    auto temp2 = make_func_eval<5>(f, -1.0, 1.0);
    auto ref2 = temp2(-0.4);
    FE move_assign = make_func_eval<5>(f, -1.0, 1.0);
    move_assign = std::move(temp2);
    EXPECT_NEAR(ref2, move_assign(-0.4), 1e-12);
}

TEST(CopyMove, FuncEvalMany) {
    auto f1 = [](double x) { return x * x; };
    auto f2 = [](double x) { return std::sin(x); };
    auto fe1 = make_func_eval<5>(f1, -1.0, 1.0);
    auto fe2 = make_func_eval<5>(f2, -1.0, 1.0);
    auto base = make_func_eval_many(fe1, fe2);
    using FEM = decltype(base);

    auto ref = base(0.3);
    FEM copy(base);
    auto vals = copy(0.3);
    EXPECT_NEAR(ref[0], vals[0], 1e-12);
    EXPECT_NEAR(ref[1], vals[1], 1e-12);

    FEM copy_assign = make_func_eval_many(fe1, fe2);
    copy_assign = base;
    vals = copy_assign(0.2);
    ref = base(0.2);
    EXPECT_NEAR(ref[0], vals[0], 1e-12);
    EXPECT_NEAR(ref[1], vals[1], 1e-12);

    auto temp = make_func_eval_many(fe1, fe2);
    ref = temp(0.5);
    FEM move_ctor(std::move(temp));
    vals = move_ctor(0.5);
    EXPECT_NEAR(ref[0], vals[0], 1e-12);
    EXPECT_NEAR(ref[1], vals[1], 1e-12);

    auto temp2 = make_func_eval_many(fe1, fe2);
    ref = temp2(-0.6);
    FEM move_assign = make_func_eval_many(fe1, fe2);
    move_assign = std::move(temp2);
    vals = move_assign(-0.6);
    EXPECT_NEAR(ref[0], vals[0], 1e-12);
    EXPECT_NEAR(ref[1], vals[1], 1e-12);
}

TEST(CopyMove, FuncEvalND) {
    using In = std::array<double, 2>;
    using Out = std::array<double, 2>;
    auto f = [](const In &pt) -> Out { return Out{pt[0] + pt[1], pt[0] - pt[1]}; };
    In low{-1.0, -1.0}, high{1.0, 1.0};
    auto base = make_func_eval<3>(f, low, high);
    using FEN = decltype(base);

    In x{0.1, -0.2};
    Out ref = base(x);
    FEN copy(base);
    auto val = copy(x);
    EXPECT_NEAR(ref[0], val[0], 1e-12);
    EXPECT_NEAR(ref[1], val[1], 1e-12);

    FEN copy_assign = make_func_eval<3>(f, low, high);
    copy_assign = base;
    val = copy_assign(x);
    EXPECT_NEAR(ref[0], val[0], 1e-12);
    EXPECT_NEAR(ref[1], val[1], 1e-12);

    auto temp = make_func_eval<3>(f, low, high);
    Out ref2 = temp(x);
    FEN move_ctor(std::move(temp));
    val = move_ctor(x);
    EXPECT_NEAR(ref2[0], val[0], 1e-12);
    EXPECT_NEAR(ref2[1], val[1], 1e-12);

    auto temp2 = make_func_eval<3>(f, low, high);
    Out ref3 = temp2(x);
    FEN move_assign = make_func_eval<3>(f, low, high);
    move_assign = std::move(temp2);
    val = move_assign(x);
    EXPECT_NEAR(ref3[0], val[0], 1e-12);
    EXPECT_NEAR(ref3[1], val[1], 1e-12);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

