#include <iostream>

#include "polyfit/polyfit.hpp"

int main() {
    constexpr auto f = [](double x) { return x * x + 2.0; };
    constexpr double eps = 1e-10;
    constexpr auto poly = poly_eval::make_func_eval<eps, -1.0, 1.0>(f);

    std::cout << poly(0.5) << '\n';
}
