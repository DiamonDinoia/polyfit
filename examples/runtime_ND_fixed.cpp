#include <array>
#include <cmath>
#include <iostream>

#include "polyfit/polyfit.hpp"

int main() {
    const auto approx = poly_eval::fit([](const std::array<double, 2> &p) {
        return std::array<double, 2>{std::cos(p[0]) + std::sin(p[1]), p[0] * p[1]};
    }, 10, {-1.0, -1.0}, {1.0, 1.0});
    const auto y = approx({0.25, -0.5});

    std::cout << y[0] << ", " << y[1] << '\n';
}
