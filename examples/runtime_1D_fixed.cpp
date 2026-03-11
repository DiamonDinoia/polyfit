#include <cmath>
#include <iostream>

#include "polyfit/polyfit.hpp"

int main() {
    const auto approx = poly_eval::fit([](double x) { return std::cos(x); }, 12, -1.0, 1.0);
    const double x = 0.5;
    std::cout << "cos(" << x << ") ~= " << approx(x) << '\n';
}
