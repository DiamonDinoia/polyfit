# polyfit

`polyfit` is a header-only C++ library for fitting polynomial approximations on bounded domains.

It exposes two main entry points:

- `poly_eval::fit(...)` builds a polynomial evaluator.
- `poly_eval::pack(...)` packs several 1D evaluators for SIMD-friendly batch evaluation.

## Requirements

- C++17 for runtime fitting and evaluation
- C++20 for constexpr fitting
- Optional C++23/C++26 improvements are enabled automatically when the compiler and standard library support them
- CMake 3.14+

## Support Matrix

| API shape | Dimensions | Runtime | `constexpr` | Notes |
| --- | --- | --- | --- | --- |
| `fit(f, nCoeffs, a, b, tags...)` | 1D | Yes | No | Fixed coefficient count |
| `fit(f, eps, a, b, tags...)` | 1D | Yes | No | Adaptive search |
| `fit<NCOEFFS>(f, a, b, tags...)` | 1D | Yes | Yes, in C++20 | Fixed coefficient count |
| `fit(f, nCoeffs, a, b)` | ND | Yes | No | Runtime-sized ND fit |
| `fit<NCOEFFS, a, b>(f)` | ND | Yes | Yes, in C++20 | Fixed-count ND fit with template bounds |
| `fit<EPS, a, b, MAX_NCOEFFS, EVAL_POINTS, ITERS>(f)` | 1D | No | Yes, when `PF_HAS_CONSTEXPR_EPS_OVERLOAD` is enabled | Compile-time epsilon search |

For the full overload and tag reference, see [`docs/API.md`](docs/API.md).
For usage-focused guides, see [`docs/guides/fit.rst`](docs/guides/fit.rst) and [`docs/guides/pack.rst`](docs/guides/pack.rst).

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
In the current CI matrix, that path is exercised on the Linux GCC jobs: `gcc`,
`gcc-13`, and `gcc-14`. It is 1D-only.

Full form:

```cpp
constexpr auto approx = poly_eval::fit<1e-12, -1.0, 1.0, 48, 200, 2>([](double x) constexpr {
    return std::sin(x) + x * x;
});
```

The template parameters are:

- `EPS`: target relative error
- `a`, `b`: domain bounds
- `MAX_NCOEFFS` (default `32`): maximum coefficient count to try
- `EVAL_POINTS` (default `100`): number of points used to validate each candidate
- `ITERS` (default `1`): refinement passes after the initial fit

## Forward-Compatible Feature Gates

`polyfit` keeps its current baseline, but newer language/library features are
used automatically when available:

- C++23: native `if consteval`, clearer compile-time dispatch, and constexpr local statics in supported compilers
- C++26-capable standard libraries: native constexpr `<cmath>` paths replace internal fallback logic automatically

These are internal compatibility improvements. Consumer-facing API requirements
do not change.

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
For tag defaults, applicability, and ordering rules, use the API reference.

ND fit:

```cpp
auto approx = poly_eval::fit([](const std::array<double, 2> &p) {
    return std::array<double, 2>{std::cos(p[0]) + std::sin(p[1]), p[0] * p[1]};
}, 10, {-1.0, -1.0}, {1.0, 1.0});

auto y = approx({0.25, -0.5});
```

For fixed-count ND fits, the template-bounds form supports constexpr
construction in C++20 when the callable is `constexpr`:

```cpp
constexpr std::array<double, 2> a{-1.0, -1.0};
constexpr std::array<double, 2> b{1.0, 1.0};
constexpr auto approx = poly_eval::fit<4, a, b>([](const std::array<double, 2> &p) constexpr {
    return std::array<double, 2>{p[0] + p[1], p[0] - p[1]};
});
static_assert(approx.nCoeffsPerAxis() == 4);
```

Runtime-sized ND fits remain runtime-only.

Pack several 1D evaluators:

```cpp
auto sinApprox = poly_eval::fit<8>([](double x) { return std::sin(x); }, -1.0, 1.0);
auto cosApprox = poly_eval::fit<8>([](double x) { return std::cos(x); }, -1.0, 1.0);
auto packed = poly_eval::pack(sinApprox, cosApprox);
auto y = packed(0.5);
```

Use `pack(...)` when you have several independent 1D evaluators of the same value type and you want to evaluate them together:

- `packed(x)`: every polynomial sees the same input
- `packed(xs)`: each polynomial sees its own input
- `packed(pts, out, count)`: evaluate many inputs into a flat output buffer

It is usually worthwhile when you evaluate the same bundle of 1D approximations repeatedly. It is not for ND evaluators and it does not help a single polynomial on its own.

## Notes

- `fit(f, eps, a, b, ...)` searches from 1 coefficient upward and throws if the error target is not met before `MaxCoeffs`.
- 1D coefficients are stored in Horner order: highest degree first.
- Degenerate domains throw `std::invalid_argument`.
- `pack(...)` is for 1D `FuncEval` objects, not `FuncEvalND`.

## Build Examples, Tests, and Docs

```bash
cmake -S . -B build -DPOLYFIT_BUILD_TESTS=ON -DPOLYFIT_BUILD_EXAMPLES=ON
cmake --build build
cmake --build build --target run_examples
ctest --test-dir build --output-on-failure
```

To build the documentation:

```bash
python3 -m pip install -r docs/requirements.txt
cmake -S . -B build-docs-check
cmake --build build-docs-check --target docs
```

The generated HTML is written under `build-docs-check/docs/html/`.

Reference: [`docs/API.md`](docs/API.md)

## Performance

Benchmarks run automatically on pushes to `main` across `gcc-14`, `gcc-15`, `llvm-20`, and `llvm-21` on Ubuntu 24.04 with `-mavx2 -mfma`.
Published charts and raw summaries live on the [`benchmark-results`](https://github.com/DiamonDinoia/polyfit/tree/benchmark-results) branch.

To run the local benchmark suite:

```bash
cmake -S . -B build-bench -DPOLYFIT_BUILD_TESTS=ON -DPOLYFIT_BUILD_EXAMPLES=OFF -DPOLYFIT_ENABLE_BENCH_NATIVE_TUNING=ON
cmake --build build-bench --target run_benchmarks
```

For stable numbers, benchmark on a tuned system:

- use a `performance` CPU governor
- run `pyperf system tune` when available
- avoid sanitizer-enabled builds for timing runs

### Horner evaluation

![Horner evaluation performance](https://raw.githubusercontent.com/DiamonDinoia/polyfit/benchmark-results/horner_performance.svg)

### 1D fitting

![1D fitting performance](https://raw.githubusercontent.com/DiamonDinoia/polyfit/benchmark-results/fitting_performance.svg)

Raw numbers: [`summary.md`](https://raw.githubusercontent.com/DiamonDinoia/polyfit/benchmark-results/summary.md)
