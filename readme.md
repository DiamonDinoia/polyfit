# polyfit

`polyfit` is a lightweight, header-only C++ library for creating fast, polynomial approximations of functions.

In mathematical terms, to approximate a function $f(x)$ over an interval $(a, b)$ with a polynomial $P(x)$ of degree $d$, the library constructs an interpolating polynomial. The process is as follows:

1.  **Node Selection**: The function $f(x)$ is sampled at $d+1$ Chebyshev nodes within the interval $[a, b]$. These nodes are chosen to minimize the maximum approximation error, preventing the large oscillations typical of interpolation with equally spaced points (Runge's phenomenon).

2.  **Newton Divided Differences (Björck–Pereyra)**: The samples are used to compute Newton divided differences via the Björck–Pereyra algorithm. This produces the polynomial in Newton form, which is numerically stable and avoids forming the ill-conditioned Vandermonde matrix.

3.  **Monomial Conversion**: The Newton-form coefficients are converted to the standard monomial basis $P(x) = \sum_{j=0}^{d} c_j x^j$ via a stable recurrence.

4.  **Optional Refinement**: Iterative refinement reruns the fitting step on the residual to further reduce numerical error. The number of refinement iterations is controlled by the `Iters_compile_time` template parameter (default: 1).

The library supports approximating both scalar ($f: \mathbb{R} \to \mathbb{R}$) and vector-valued ($f: \mathbb{R}^N \to \mathbb{R}^M$) functions.

## Requirements

*   A C++17 compatible compiler (GCC, Clang, or MSVC).
*   A C++20 compatible compiler is required for compile-time approximation features.
*   CMake version 3.14 or higher.
*   **MSVC note**: C++17 is fully supported. For C++20 `constexpr` fitting with MSVC, add `/Zc:__cplusplus` to your compiler flags.

## How to use with CMake

`polyfit` is designed for easy integration with CMake. You can link against the `polyfit::polyfit` target.

### Using `FetchContent`

This is the recommended approach. Add the following to your `CMakeLists.txt`:

```cmake
include(FetchContent)
FetchContent_Declare(
    polyfit
    GIT_REPOSITORY https://github.com/DiamonDinoia/polyfit.git
    GIT_TAG        main
)
FetchContent_MakeAvailable(polyfit)

# ... later in your CMakeLists.txt
add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE polyfit::polyfit)
```

### Using `find_package`

If you have installed `polyfit` on your system, you can use `find_package`:

```cmake
find_package(polyfit REQUIRED)

# ... later in your CMakeLists.txt
add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE polyfit::polyfit)
```

## API & Examples

The core of the library is the `poly_eval::make_func_eval` factory function, which creates a callable polynomial object that approximates a given function.

Full API reference: [docs/API.md](docs/API.md)

### 1D Function Approximation (Runtime)

Create an approximation for a scalar function $f: \mathbb{R} \to \mathbb{R}$ by specifying either a fixed degree or an error tolerance at runtime.

#### Fixed Degree

```cpp
#include <polyfit/fast_eval.hpp>
#include <cmath>
#include <iostream>

int main() {
    auto my_func = [](double x) { return std::cos(x); };
    double a = -1.0, b = 1.0;
    int degree = 16;

    auto poly = poly_eval::make_func_eval(my_func, degree, a, b);
    std::cout << "my_func(0.5) ≈ " << poly(0.5) << std::endl;
}
```

#### Error Tolerance

```cpp
#include <polyfit/fast_eval.hpp>
#include <cmath>
#include <iostream>

int main() {
    auto my_func = [](double x) { return std::cos(x); };
    double a = -1.0, b = 1.0;
    double epsilon = 1e-12;

    auto poly = poly_eval::make_func_eval(my_func, epsilon, a, b);
    std::cout << "my_func(0.7) ≈ " << poly(0.7) << std::endl;
}
```

### 1D Function Approximation (Compile-time, C++20)

For maximum performance, the polynomial can be generated at compile time if the function, domain, and degree/epsilon are `constexpr`.

#### Fixed Degree

```cpp
#include <polyfit/fast_eval.hpp>
#include <cmath>
#include <iostream>

int main() {
    constexpr auto my_func = [](double x) { return x*x+2; };
    constexpr double a = 0.0, b = 3.14159;
    constexpr size_t degree = 12;

    constexpr auto poly = poly_eval::make_func_eval<degree>(my_func, a, b);
    // if c++20 is not available, use:
    // auto poly = poly_eval::make_func_eval<degree>(my_func, a, b);
    // which will do the fitting at runtime. But generates a faster evaluator than the runtime version.

    std::cout << "my_func(1.0) ≈ " << poly(1.0) << std::endl;
}
```

#### Error Tolerance

```cpp
#include <polyfit/fast_eval.hpp>
#include <cmath>
#include <iostream>

int main() {
    constexpr auto my_func = [](double x) { return  x*x+2; };
    constexpr double a = 0.0, b = 3.14159;
    constexpr double eps = 1e-8;

    constexpr auto poly = poly_eval::make_func_eval<eps, a, b>(my_func);
    std::cout << "my_func(1.0) ≈ " << poly(1.0) << std::endl;
}
```

### N-D Function Approximation (Runtime)

The same API can be used to approximate vector-valued functions $f: \mathbb{R}^N \to \mathbb{R}^M$.

#### Fixed Degree

```cpp
#include <polyfit/fast_eval.hpp>
#include <array>
#include <cmath>
#include <iostream>

int main() {
    using In = std::array<double, 2>;
    using Out = std::array<double, 2>;

    auto f = [](const In& p) {
        return Out{std::cos(p[0]) + std::sin(p[1]), std::exp(p[0] * p[1])};
    };
    In a{-1.0, -1.0}, b{1.0, 1.0};
    int degree = 12;

    auto poly = poly_eval::make_func_eval(f, degree, a, b);

    // if the degree is known at compile time, this version generates a faster evaluator:
    // auto poly = poly_eval::make_func_eval<degree>(f, a, b);

    In point = {0.5, -0.2};
    Out result = poly(point);
    std::cout << "f(0.5, -0.2) ≈ {" << result[0] << ", " << result[1] << "}" << std::endl;
}
```

#### Error Tolerance

```cpp
#include <polyfit/fast_eval.hpp>
#include <array>
#include <cmath>
#include <iostream>

int main() {
    using In = std::array<double, 2>;
    using Out = std::array<double, 2>;

    auto f = [](const In& p) {
        return Out{std::cos(p[0]) + std::sin(p[1]), std::exp(p[0] * p[1])};
    };
    In a{-1.0, -1.0}, b{1.0, 1.0};
    double eps = 1e-3;

    auto poly = poly_eval::make_func_eval(f, eps, a, b);
    // if the eps is known at compile time and c++20 is availble, this version generates a faster evaluator:
    // auto poly = poly_eval::make_func_eval<eps>(f, a, b);

    In point = {0.5, -0.2};
    Out result = poly(point);
    std::cout << "f(0.5, -0.2) ≈ {" << result[0] << ", " << result[1] << "}" << std::endl;
}
```

## Building the Examples

The project includes several examples. You can build and run them all using the `run_examples` CMake target.

```bash
# Clone the repository and navigate to the project directory then

# Configure the project
cmake -S . -B build

# Build and run the examples
cmake --build build --target run_examples
```

## Accuracy & Performance

### Compensated Arithmetic

Both `bjorck_pereyra` (Newton divided differences) and `newton_to_monomial` (monomial conversion) use **error-free transformations (EFT) / compensated summation**. This yields near-machine-precision fitting accuracy even at high polynomial degrees where naive floating-point accumulation would lose significant digits.

### Refinement Iterations

The `Iters_compile_time` template parameter on `FuncEval` controls the number of iterative refinement steps applied after the initial fit (default: 1). Each iteration reruns the fitting pipeline on the current residual and corrects accumulated floating-point error. For high-degree polynomials or stringent accuracy requirements, increasing `Iters_compile_time` to 2–3 can meaningfully reduce the approximation error. Refinement uses a compensated Horner scheme internally to avoid divergence at high degrees.

### Truncation

When using the error-tolerance overload of `make_func_eval`, the evaluator may drop trailing near-zero coefficients via `truncate(eps)` after fitting converges. This reduces the polynomial degree actually evaluated, lowering per-call cost without exceeding the requested error bound.

### SIMD Bulk Evaluation

Bulk evaluation paths use SIMD intrinsics (via xsimd) automatically. For best throughput, prefer the bulk overload:

```cpp
poly(pts, out, num_points);  // operator()(const T* pts, T* out, std::size_t n)
```

For evaluating multiple polynomials over the same inputs simultaneously, use `FuncEvalMany` (see API reference).

## Notes on Coefficients and Mapping

- Coefficient order for Horner’s method is reversed: pass `[c_N, …, c_0]` to the low-level 1D APIs.
- Evaluators internally map the user domain `[a,b]` to the Chebyshev domain for sampling and back for evaluation.

## Testing Across C++ Standards

- The test suite builds and runs for C++17 and C++20.
- Enable tests and build:
  - `cmake -S . -B build -DPOLYFIT_BUILD_TESTS=ON`
  - `cmake --build build`
  - `ctest --test-dir build --output-on-failure`
- Per-standard test targets are generated (e.g., `test_horner_cxx17`, `test_horner_cxx20`).
- Benchmarks follow the same pattern (e.g., `bench_horner_cxx17`, `bench_horner_cxx20`).
