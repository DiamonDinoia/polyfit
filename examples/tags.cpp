#include <cmath>
#include <iostream>

#include "polyfit/polyfit.hpp"

int main() {
    const auto approx = poly_eval::fit(
        [](double x) { return std::sin(x); },
        1e-12,
        -1.0,
        1.0,
        poly_eval::MaxCoeffs<48>{},
        poly_eval::EvalPts<200>{},
        poly_eval::Iters<2>{},
        poly_eval::FuseNever{});

    std::cout << approx(0.25) << '\n';
}
