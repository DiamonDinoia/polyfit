#include <iostream>

#include "polyfit/polyfit.hpp"

int main() {
    constexpr auto approx = poly_eval::fit<3>([](double x) { return x * x + 2.0; }, -1.0, 1.0);
    std::cout << approx(0.5) << '\n';
}
