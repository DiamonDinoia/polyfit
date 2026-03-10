#include <array>
#include <cmath>
#include <iostream>

#include "polyfit/polyfit.hpp"

int main() {
    using In = std::array<double, 2>;
    using Out = std::array<double, 2>;

    const auto f = [](const In &p) {
        return Out{std::cos(p[0]) + std::sin(p[1]), p[0] * p[1]};
    };

    const auto poly = poly_eval::make_func_eval(f, 1e-6, In{-1.0, -1.0}, In{1.0, 1.0});
    const auto y = poly(In{0.25, -0.5});

    std::cout << y[0] << ", " << y[1] << '\n';
}
