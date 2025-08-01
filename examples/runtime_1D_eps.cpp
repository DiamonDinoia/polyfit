#include <cmath>
#include <iostream>
#include <vector>

#include "polyfit/fast_eval.hpp"

int main() {
    // 1. Define the function, domain, and desired error tolerance.
    auto my_func = [](double x) { return std::cos(x); };
    double a = -1.0, b = 1.0;
    double epsilon = 1e-12;

    // 2. Create the approximation. The degree is chosen automatically
    // to meet the error tolerance 'epsilon'.
    auto poly = poly_eval::make_func_eval(my_func, epsilon, a, b);

    // 3. Verify the approximation at a random point.
    double x = 0.7;
    double y_approx = poly(x);
    double y_true = my_func(x);

    std::cout << "Error-driven evaluation for f(x) = cos(x) with eps = " << epsilon << std::endl;
    std::cout << "  - Approximation at x=0.7: " << y_approx << std::endl;
    std::cout << "  - True value:             " << y_true << std::endl;
    std::cout << "  - Error:                  " << std::abs(1 - y_approx / y_true) << std::endl;

    std::vector<double> xs = {0.1, 0.2, 0.3, 0.4};
    std::vector<double> ys(xs.size());
    poly(xs.data(), ys.data(), xs.size());

    std::cout << "\nBatch evaluation:" << std::endl;
    for (size_t i = 0; i < xs.size(); ++i) {
        std::cout << "  - f(" << xs[i] << ") ≈ " << ys[i] << std::endl;
    }

    return 0;
}