#include <cmath>
#include <iostream>

#include "polyfit/polyfit.hpp"

int main() {
#if PF_HAS_CONSTEXPR_EPS_OVERLOAD
    constexpr auto approx = poly_eval::fit<1e-12, -1.0, 1.0>([](double x) constexpr {
        return std::sin(x) + x * x;
    });
    static_assert(approx(0.25) > 0.0);
    std::cout << approx(0.25) << '\n';
#else
    std::cout << "compile-time eps fit is unavailable on this compiler\n";
#endif
}
