# polyfit

`polyfit` is a header-only C++ library for fitting polynomial approximations on bounded domains.

It exposes two main entry points:

- `poly_eval::fit(...)` builds a polynomial evaluator.
- `poly_eval::pack(...)` packs several 1D evaluators for SIMD-friendly batch evaluation.

## Requirements

- C++17 for runtime fitting and evaluation
- C++20 for constexpr fitting
- CMake 3.14+

## Install with CMake

```cmake
include(FetchContent)

FetchContent_Declare(
    polyfit
    GIT_REPOSITORY https://github.com/DiamonDinoia/polyfit.git
    GIT_TAG main
)
FetchContent_MakeAvailable(polyfit)

target_link_libraries(my_target PRIVATE polyfit::polyfit)
```

## Quick Start

```cpp
#include <cmath>
#include <polyfit/polyfit.hpp>

int main() {
    auto approx = poly_eval::fit([](double x) { return std::cos(x); }, 12, -1.0, 1.0);
    return approx(0.5) > 0.0 ? 0 : 1;
}
```

Compile-time coefficient count:

```cpp
constexpr auto approx = poly_eval::fit<8>([](double x) { return x * x + 1.0; }, -1.0, 1.0);
static_assert(approx(0.5) > 1.0);
```

Compile-time fit with a compile-time error target:

```cpp
constexpr auto approx = poly_eval::fit<1e-12, -1.0, 1.0>([](double x) constexpr {
    return std::sin(x) + x * x;
});
static_assert(approx(0.25) > 0.0);
```

This overload is available only when `PF_HAS_CONSTEXPR_EPS_OVERLOAD` is enabled.
Today that means the supported C++20 GCC path. It is 1D-only and uses template
parameters instead of runtime tags:

- `EPS`
- `a`
- `b`
- optional `MAX_NCOEFFS`, `EVAL_POINTS`, `ITERS`

Adaptive fit with options:

```cpp
auto approx = poly_eval::fit(
    [](double x) { return std::sin(x); },
    1e-12,
    -1.0,
    1.0,
    poly_eval::MaxCoeffs<48>{},
    poly_eval::EvalPts<200>{},
    poly_eval::Iters<2>{},
    poly_eval::FuseNever{});
```

The tags can appear in any order. You can pass only the ones you need.

ND fit:

```cpp
auto approx = poly_eval::fit([](const std::array<double, 2> &p) {
    return std::array<double, 2>{std::cos(p[0]) + std::sin(p[1]), p[0] * p[1]};
}, 10, {-1.0, -1.0}, {1.0, 1.0});

auto y = approx({0.25, -0.5});
```

Pack several 1D evaluators:

```cpp
auto sinApprox = poly_eval::fit<8>([](double x) { return std::sin(x); }, -1.0, 1.0);
auto cosApprox = poly_eval::fit<8>([](double x) { return std::cos(x); }, -1.0, 1.0);
auto packed = poly_eval::pack(sinApprox, cosApprox);
auto y = packed(0.5);
```

## Notes

- `fit(f, eps, a, b, ...)` searches from 1 coefficient upward and throws if the error target is not met before `MaxCoeffs`.
- 1D coefficients are stored in Horner order: highest degree first.
- Degenerate domains throw `std::invalid_argument`.
- `pack(...)` is for 1D `FuncEval` objects, not `FuncEvalND`.

## Tags

- Tags are optional.
- Tags can appear in any order.
- Repeating the same tag kind is a compile-time error.
- `MaxCoeffs<N>{}` sets the search cap for adaptive `fit(f, eps, a, b, ...)`.
- `EvalPts<N>{}` sets how many points are checked when validating an adaptive fit.
- `Iters<N>{}` controls refinement iterations for 1D fitting.
- `FuseAuto{}`, `FuseAlways{}`, and `FuseNever{}` control domain fusion in 1D fits.
- `MaxCoeffs` and `EvalPts` do not affect fixed-count or ND runtime fits.
- The compile-time epsilon overload does not use tags; it takes `MAX_NCOEFFS`,
  `EVAL_POINTS`, and `ITERS` as template parameters.

Example with the same tags in a different order:

```cpp
auto approx = poly_eval::fit(
    [](double x) { return std::sin(x); },
    1e-12,
    -1.0,
    1.0,
    poly_eval::FuseNever{},
    poly_eval::Iters<2>{},
    poly_eval::EvalPts<200>{},
    poly_eval::MaxCoeffs<48>{});
```

## Build Examples and Tests

```bash
cmake -S . -B build -DPOLYFIT_BUILD_TESTS=ON -DPOLYFIT_BUILD_EXAMPLES=ON
cmake --build build
cmake --build build --target run_examples
ctest --test-dir build --output-on-failure
```

More detail: [`docs/API.md`](docs/API.md)

## Performance

Benchmarks run automatically on pushes to `main` across `gcc-14`, `gcc-15`, `llvm-20`, and `llvm-21` on Ubuntu 24.04 with `-mavx2 -mfma`.
Published charts and raw summaries live on the [`benchmark-results`](https://github.com/DiamonDinoia/polyfit/tree/benchmark-results) branch.

### Horner evaluation

![Horner evaluation performance](https://raw.githubusercontent.com/DiamonDinoia/polyfit/benchmark-results/horner_performance.svg)

### 1D fitting

![1D fitting performance](https://raw.githubusercontent.com/DiamonDinoia/polyfit/benchmark-results/fitting_performance.svg)

Raw numbers: [`summary.md`](https://raw.githubusercontent.com/DiamonDinoia/polyfit/benchmark-results/summary.md)
