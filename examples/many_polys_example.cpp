#include <cmath>
#include <iostream>

#include "polyfit/polyfit.hpp"

int main() {
    const auto sinPoly = poly_eval::make_func_eval<8>([](double x) { return std::sin(x); }, -1.0, 1.0);
    const auto cosPoly = poly_eval::make_func_eval<8>([](double x) { return std::cos(x); }, -1.0, 1.0);
    const auto packed = poly_eval::make_func_eval_many(sinPoly, cosPoly);
    const auto y = packed(0.5);

    std::cout << y[0] << ", " << y[1] << '\n';
}
