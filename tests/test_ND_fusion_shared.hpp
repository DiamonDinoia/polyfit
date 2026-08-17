#pragma once

#include <array>
#include <cmath>
#include <cstddef>

namespace polyfit_test_nd_fusion {

constexpr std::size_t kNumRandomTests = 4000;

using Arr2 = std::array<double, 2>;
using Out1 = std::array<double, 1>;
using Out2 = std::array<double, 2>;

inline auto smooth_2d = [](const Arr2 &x) -> Out1 {
    return {std::sin(1.3 * x[0]) * std::cos(0.7 * x[1])};
};

} // namespace polyfit_test_nd_fusion
