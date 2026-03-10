#include <iostream>

#include "polyfit/polyfit.hpp"

int main() {
    constexpr auto f = [](double x) { return x * x + 2.0; };
    constexpr auto poly = poly_eval::make_func_eval<3>(f, -1.0, 1.0);

    std::cout << poly(0.5) << '\n';
}
