#include <cmath>
#include <iostream>

#include "polyfit/polyfit.hpp"

int main() {
    const auto f = [](double x) { return std::cos(x); };
    const double eps = 1e-12;
    const auto poly = poly_eval::make_func_eval(f, eps, -1.0, 1.0);
    const double x = 0.5;
    std::cout << "f(" << x << ") ~= " << poly(x) << '\n';
}
