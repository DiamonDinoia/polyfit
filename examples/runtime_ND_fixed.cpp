#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "polyfit/fast_eval.hpp"

// --- Helper to print results ---
template <typename Func, typename Poly, typename Input>
void test_and_print(const std::string &title, Func f, const Poly &poly, const Input &pt) {
    auto y_approx = poly(pt);
    auto y_true = f(pt);
    double error = poly_eval::detail::relative_l2_norm(y_approx, y_true);

    std::cout << title << std::endl;
    std::cout << "  - Approximation at point: " << y_approx[0] << ", ...\n";
    std::cout << "  - True value at point:    " << y_true[0] << ", ...\n";
    std::cout << "  - Relative L2 Error:      " << error << std::endl << std::endl;
}

int main() {
    // =========================================================================
    // Example 1: F: R^2 -> R^2 with runtime-specified degree
    // =========================================================================
    {
        using In = std::array<double, 2>;
        using Out = std::array<double, 2>;
        auto f = [](const In &p) { return Out{std::cos(p[0]) + std::sin(p[1]), std::exp(p[0] * p[1])}; };
        In a{-1.0, -1.0}, b{1.0, 1.0};
        int degree = 12;
        auto poly = poly_eval::make_func_eval(f, degree, a, b);
        test_and_print("F:R^2->R^2, Runtime Degree=12", f, poly, In{0.5, -0.2});
    }

    // =========================================================================
    // Example 2: F: R^2 -> R^3 with runtime-specified degree
    // =========================================================================
    {
        using In = std::array<double, 2>;
        using Out = std::array<double, 3>;
        auto f = [](const In &p) { return Out{std::cos(p[0]), std::sin(p[1]), p[0] + p[1]}; };
        In a{-1.0, -1.0}, b{1.0, 1.0};
        int degree = 10;
        auto poly = poly_eval::make_func_eval(f, degree, a, b);
        test_and_print("F:R^2->R^3, Runtime Degree=10", f, poly, In{0.1, 0.8});
    }

    // =========================================================================
    // Example 3: F: R^3 -> R^2 with runtime-specified degree
    // =========================================================================
    {
        using In = std::array<double, 3>;
        using Out = std::array<double, 2>;
        auto f = [](const In &p) { return Out{p[0] * p[1] + p[2], std::exp(p[0] + p[1] - p[2])}; };
        In a{-1.0, -1.0, -1.0}, b{1.0, 1.0, 1.0};
        int degree = 8;
        auto poly = poly_eval::make_func_eval(f, degree, a, b);
        test_and_print("F:R^3->R^2, Runtime Degree=8", f, poly, In{0.3, -0.4, 0.5});
    }

    return 0;
}