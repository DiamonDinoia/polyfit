#include <iostream>
#include "polyfit/fast_eval.hpp"

int main() {
    using namespace poly_eval;

    auto f1 = [](double x) { return std::sin(x); };
    auto f2 = [](double x) { return std::cos(2 * x); };

    auto fe1 = make_func_eval<8>(f1, -1.0, 1.0);
    auto fe2 = make_func_eval<8>(f2, -1.0, 1.0);

    auto packed = make_func_eval_many(fe1, fe2);

    double x = 0.5;
    auto out = packed(x);
    std::cout << "f1(" << x << ") approx= " << out[0] << "\n";
    std::cout << "f2(" << x << ") approx= " << out[1] << "\n";

    // Bulk evaluation example
    std::array<double, 2> xs = {0.1, 0.2};
    auto out2 = packed(xs);
    std::cout << "packed(xs)[0] = " << out2[0] << "\n";
    std::cout << "packed(xs)[1] = " << out2[1] << "\n";
    return 0;
}
