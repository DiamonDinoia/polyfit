# polyfit API

## Headers

- `#include <polyfit/polyfit.hpp>`: umbrella header
- `#include <polyfit/polyeval.hpp>`: public API header

## Main Entry Points

### `poly_eval::fit(...)`

Supported overloads:

- `poly_eval::fit(f, nCoeffs, a, b, tags...)`
- `poly_eval::fit(f, errorTarget, a, b, tags...)`
- `poly_eval::fit<NCOEFFS>(f, a, b, tags...)`
- `poly_eval::fit<NCOEFFS, a, b>(f)` for ND fitting with template-parameter bounds
- `poly_eval::fit<EPS, a, b, MAX_NCOEFFS, EVAL_POINTS, ITERS>(f)` for compile-time epsilon-driven 1D fitting when `PF_HAS_CONSTEXPR_EPS_OVERLOAD` is enabled

Returned type:

- `FuncEval` for 1D scalar or complex callables
- `FuncEvalND` for `std::array`-like inputs

### `poly_eval::pack(...)`

Packs several 1D `FuncEval` objects into a `FuncEvalMany`.

## `fit(...)` Overloads

### Runtime fixed coefficient count

```cpp
auto approx = poly_eval::fit(
    [](double x) { return std::cos(x); },
    12,
    -1.0,
    1.0);
```

Notes:

- works for 1D and ND callables
- `nCoeffs` must be positive
- 1D overloads accept the `Iters` and fusion tags
- ND runtime fixed-count fits ignore the 1D-only tuning tags

### Runtime adaptive fit

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

Notes:

- 1D only
- searches upward from `1` coefficient
- throws if the requested error is not met before `MaxCoeffs`

### Compile-time fixed coefficient count

```cpp
constexpr auto approx = poly_eval::fit<8>(
    [](double x) { return x * x + 1.0; },
    -1.0,
    1.0,
    poly_eval::Iters<2>{});
```

Notes:

- 1D overloads accept `Iters` and fusion tags
- `NCOEFFS` must be positive

### ND fit with template-parameter bounds

```cpp
constexpr std::array<double, 2> a{-1.0, -1.0};
constexpr std::array<double, 2> b{1.0, 1.0};

auto approx = poly_eval::fit<8, a, b>([](const std::array<double, 2> &p) {
    return std::array<double, 2>{p[0] + p[1], p[0] * p[1]};
});
```

Notes:

- C++20 only
- for `std::array`-like inputs
- uses template parameters for the domain bounds and coefficient count

### Compile-time epsilon fit

```cpp
constexpr auto approx = poly_eval::fit<1e-12, -1.0, 1.0>([](double x) constexpr {
    return std::sin(x) + x * x;
});
```

Full form:

```cpp
constexpr auto approx = poly_eval::fit<1e-12, -1.0, 1.0, 48, 200, 2>([](double x) constexpr {
    return std::sin(x) + x * x;
});
```

Defaults:

- `MAX_NCOEFFS = 32`
- `EVAL_POINTS = 100`
- `ITERS = 1`

Notes:

- C++20 only
- available only when `PF_HAS_CONSTEXPR_EPS_OVERLOAD` is enabled
- currently intended for the supported GCC constexpr path
- 1D only
- uses template parameters, not runtime tags

## Tags

Tags apply to the 1D runtime fixed-count, runtime adaptive, and compile-time fixed-count overloads.

Rules:

- tags are optional
- tags can appear in any order
- duplicate tag kinds are rejected at compile time

Available tags:

- `poly_eval::Iters<N>`
- `poly_eval::MaxCoeffs<N>`
- `poly_eval::EvalPts<N>`
- `poly_eval::FuseAuto`
- `poly_eval::FuseAlways`
- `poly_eval::FuseNever`

Defaults:

- `Iters<1>`
- `MaxCoeffs<32>`
- `EvalPts<100>`
- `FuseAuto{}`

Meaning:

- `Iters<N>`: refinement passes for 1D fitting
- `MaxCoeffs<N>`: search cap for adaptive `fit(f, eps, a, b, ...)`
- `EvalPts<N>`: validation sample count for adaptive `fit(f, eps, a, b, ...)`
- fusion tags: control whether the domain mapping is baked into the polynomial coefficients

Not used by:

- runtime ND fixed-count fits
- compile-time epsilon fits

Equivalent tag orderings:

```cpp
auto p1 = poly_eval::fit(
    f, 1e-12, -1.0, 1.0,
    poly_eval::MaxCoeffs<48>{},
    poly_eval::EvalPts<200>{},
    poly_eval::Iters<2>{},
    poly_eval::FuseNever{});

auto p2 = poly_eval::fit(
    f, 1e-12, -1.0, 1.0,
    poly_eval::FuseNever{},
    poly_eval::Iters<2>{},
    poly_eval::EvalPts<200>{},
    poly_eval::MaxCoeffs<48>{});
```

## Main Types

### `FuncEval`

1D evaluator returned by the 1D `fit(...)` overloads.

Useful members:

- `operator()(x)` evaluates one point
- `operator()(pts, out, count)` evaluates many points
- `coeffs()` returns coefficients in Horner order
- `nCoeffs()` returns the active coefficient count
- `truncate(eps)` removes leading small coefficients for runtime-sized fits

Useful constants:

- `FuncEval::NCOEFFS`
- `FuncEval::ITERS`

### `FuncEvalND`

ND evaluator returned by ND `fit(...)` overloads.

Useful members:

- `operator()(x)` evaluates one point
- `nCoeffsPerAxis()` returns the active coefficient count used on each axis

### `FuncEvalMany`

Packed evaluator returned by `pack(...)`.

Useful members:

- `size()` returns the number of packed evaluators
- `nCoeffs()` returns the shared active coefficient count
- `operator()(x)` evaluates all polynomials at the same point
- `operator()(xs)` evaluates each polynomial at its corresponding point
- `operator()(first, rest...)` evaluates from separate scalar arguments
- `operator()(tuple)` evaluates from a tuple of inputs
- `operator()(pts, out, count)` evaluates many points into a flat output buffer
- `truncate(eps)` trims trailing coefficient rows

Useful constants:

- `FuncEvalMany::COUNT`
- `FuncEvalMany::SIMD_WIDTH`
- `FuncEvalMany::PADDED_COUNT`
- `FuncEvalMany::MAX_NCOEFFS`

## `pack(...)` Constraints

- `pack(...)` is for 1D `FuncEval` objects
- all packed evaluators must use the same sizing mode
- mixing runtime-sized and compile-time-sized evaluators is unsupported
- runtime-sized evaluators can be packed together only if they use the same active coefficient count

Basic example:

```cpp
auto sinApprox = poly_eval::fit<8>([](double x) { return std::sin(x); }, -1.0, 1.0);
auto cosApprox = poly_eval::fit<8>([](double x) { return std::cos(x); }, -1.0, 1.0);
auto packed = poly_eval::pack(sinApprox, cosApprox);

auto samePoint = packed(0.5);
auto perPoly = packed(std::array<double, 2>{0.25, 0.75});
```
