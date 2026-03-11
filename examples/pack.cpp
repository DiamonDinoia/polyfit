#include <cmath>
#include <iostream>

#include "polyfit/polyfit.hpp"

int main() {
    const auto sinApprox = poly_eval::fit<8>([](double x) { return std::sin(x); }, -1.0, 1.0);
    const auto cosApprox = poly_eval::fit<8>([](double x) { return std::cos(x); }, -1.0, 1.0);
    const auto packed = poly_eval::pack(sinApprox, cosApprox);
    const auto y = packed(0.5);

    std::cout << y[0] << ", " << y[1] << '\n';
}
