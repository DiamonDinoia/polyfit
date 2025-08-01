#include <iostream>
#include <vector>
#include <cmath>
#include "polyfit/fast_eval.hpp"

int main() {
    // 1. Define the function to approximate and its domain.
    auto my_func = [](double x) { return std::cos(x); };
    double a = -1.0, b = 1.0;
    int degree = 16;

    // 2. Create the polynomial approximation object.
    // The degree is a runtime variable.
    auto poly = poly_eval::make_func_eval(my_func, degree, a, b);

    // 3. Evaluate at a single point.
    double x = 0.5;
    double y_approx = poly(x);
    double y_true = my_func(x);

    std::cout << "Single point evaluation at x = " << x << std::endl;
    std::cout << "  - Approximation: " << y_approx << std::endl;
    std::cout << "  - True value:    " << y_true << std::endl;
    std::cout << "  - Error:         " << std::abs(1 - y_approx / y_true) << std::endl;

    // 4. Evaluate a batch of points.
    std::vector<double> xs = {0.1, 0.2, 0.3, 0.4};
    std::vector<double> ys(xs.size());
    poly(xs.data(), ys.data(), xs.size());

    std::cout << "\nBatch evaluation:" << std::endl;
    for (size_t i = 0; i < xs.size(); ++i) {
        std::cout << "  - f(" << xs[i] << ") ≈ " << ys[i] << std::endl;
    }

    return 0;
}