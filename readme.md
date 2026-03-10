# polyfit

`polyfit` is a header-only C++ library for fitting fast polynomial approximations of functions on a bounded domain.

It builds an interpolating polynomial on Chebyshev nodes, computes Newton divided differences with Björck-Pereyra, converts to monomial coefficients, and evaluates with Horner's method.

## Requirements

- C++17 for runtime fitting and evaluation
- C++20 for `constexpr` fitting
- CMake 3.14+

## Install with CMake

### `FetchContent`

```cmake
include(FetchContent)

FetchContent_Declare(
    polyfit
    GIT_REPOSITORY https://github.com/DiamonDinoia/polyfit.git
    GIT_TAG main
)
FetchContent_MakeAvailable(polyfit)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE polyfit::polyfit)
```

### `find_package`

```cmake
find_package(polyfit REQUIRED)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE polyfit::polyfit)
```

## Public API

Use these headers:

- `#include <polyfit/polyfit.hpp>`
- `#include <polyfit/polyeval.hpp>`

`polyfit.hpp` is the umbrella header. `polyeval.hpp` contains the public API.

The main entry point is `poly_eval::make_func_eval(...)`.

## Examples

### Runtime 1D fit with a fixed coefficient count

```cpp
#include <cmath>
#include <iostream>
#include <polyfit/polyfit.hpp>

int main() {
    auto f = [](double x) { return std::cos(x); };
    double a = -1.0;
    double b = 1.0;
    int nCoeffs = 12;

    auto poly = poly_eval::make_func_eval(f, nCoeffs, a, b);
    std::cout << poly(0.5) << '\n';
}
```

### Runtime 1D fit with an error target

```cpp
#include <cmath>
#include <iostream>
#include <polyfit/polyfit.hpp>

int main() {
    auto f = [](double x) { return std::cos(x); };
    auto poly = poly_eval::make_func_eval(f, 1e-12, -1.0, 1.0);
    std::cout << poly(0.5) << '\n';
}
```

### Compile-time 1D fit

```cpp
#include <iostream>
#include <polyfit/polyfit.hpp>

int main() {
    constexpr auto f = [](double x) { return x * x + 2.0; };
    constexpr auto poly = poly_eval::make_func_eval<3>(f, -1.0, 1.0);
    std::cout << poly(0.5) << '\n';
}
```

### Runtime N-D fit

```cpp
#include <array>
#include <iostream>
#include <polyfit/polyfit.hpp>

int main() {
    using In = std::array<double, 2>;
    using Out = std::array<double, 2>;

    auto f = [](const In& p) {
        return Out{std::cos(p[0]) + std::sin(p[1]), p[0] * p[1]};
    };

    auto poly = poly_eval::make_func_eval(f, 10, In{-1.0, -1.0}, In{1.0, 1.0});
    auto y = poly(In{0.25, -0.5});

    std::cout << y[0] << ", " << y[1] << '\n';
}
```

### Pack several 1D evaluators

```cpp
#include <cmath>
#include <iostream>
#include <polyfit/polyfit.hpp>

int main() {
    auto sinPoly = poly_eval::make_func_eval<8>([](double x) { return std::sin(x); }, -1.0, 1.0);
    auto cosPoly = poly_eval::make_func_eval<8>([](double x) { return std::cos(x); }, -1.0, 1.0);

    auto packed = poly_eval::make_func_eval_many(sinPoly, cosPoly);
    auto y = packed(0.5);

    std::cout << y[0] << ", " << y[1] << '\n';
}
```

More examples: `examples/`

Full API reference: [docs/API.md](docs/API.md)

## Notes

- Runtime `make_func_eval(func, eps, a, b)` searches from coefficient count `1` up to `MaxN` and returns the first fit that satisfies the checked error bound.
- Coefficients are stored in Horner order: highest degree first.
- Degenerate domains are rejected.
- Compensated arithmetic is used in the main fitting kernels, not everywhere.

## Build examples and tests

```bash
cmake -S . -B build -DPOLYFIT_BUILD_TESTS=ON
cmake --build build
ctest --test-dir build --output-on-failure
```
