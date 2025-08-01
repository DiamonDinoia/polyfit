#include <cmath>
#include <iostream>
#include <vector>

#include "polyfit/fast_eval.hpp"

int main() {
    // 1. Define the function and domain.
    constexpr auto my_func = [](double x) { return 2.0 * x * x * x - 3.0 * x + 1.0; };
    double a = -1.0, b = 1.0;

    // 2. Create the approximation with a compile-time fixed degree of 16.
    constexpr size_t degree = 16;
    auto poly = poly_eval::make_func_eval<degree>(my_func, a, b);

    // 3. Evaluate at a single point.
    double x = 0.5;
    double y_approx = poly(x);
    double y_true = my_func(x);

    std::cout << "Compile-time degree evaluation (N=" << degree << ")" << std::endl;
    std::cout << "  - Approximation at x=0.5: " << y_approx << std::endl;
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