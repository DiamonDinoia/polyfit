# polyfit API

## Headers

- `#include <polyfit/polyfit.hpp>`: umbrella header
- `#include <polyfit/polyeval.hpp>`: public API header

## Main Entry Points

For usage-first examples, see [`guides/fit.rst`](guides/fit.rst) and [`guides/pack.rst`](guides/pack.rst).

### `poly_eval::fit(...)`

Supported overloads:

- `poly_eval::fit(f, nCoeffs, a, b, tags...)`
- `poly_eval::fit(f, errorTarget, a, b, tags...)`
- `poly_eval::fit<NCOEFFS>(f, a, b, tags...)`
- `poly_eval::fit<NCOEFFS, a, b>(f)` for ND fitting with template-parameter bounds
- `poly_eval::fit<EPS, a, b, MAX_NCOEFFS, EVAL_POINTS, ITERS>(f)` for compile-time epsilon-driven 1D fitting when `PF_HAS_CONSTEXPR_EPS_OVERLOAD` is enabled

Returned type:

- `FuncEval` for 1D scalar or complex callables
- `FuncEvalND` for fixed-size indexable ND inputs and outputs

Support matrix:

| Form | Dimensions | Runtime | `constexpr` | Notes |
| --- | --- | --- | --- | --- |
| `fit(f, nCoeffs, a, b, tags...)` | 1D | Yes | No | Fixed-count |
| `fit(f, eps, a, b, tags...)` | 1D | Yes | No | Adaptive search |
| `fit<NCOEFFS>(f, a, b, tags...)` | 1D | Yes | Yes, in C++20 | Fixed-count |
| `fit(f, nCoeffs, a, b)` | ND | Yes | No | Runtime-sized ND |
| `fit<NCOEFFS, a, b>(f)` | ND | Yes | Yes, in C++20 | Fixed-count ND with template bounds |
| `fit<EPS, a, b, MAX_NCOEFFS, EVAL_POINTS, ITERS>(f)` | 1D | No | Yes, when `PF_HAS_CONSTEXPR_EPS_OVERLOAD` is enabled | Compile-time epsilon search |

Tested support:

| Area | Current CI coverage |
| --- | --- |
| Runtime 1D, runtime ND, `pack(...)` | Linux `gcc`, `gcc-13`, `gcc-14`, `llvm`, `llvm-18`, `llvm-21`; macOS `apple-clang`; Windows `msvc` |
| Fixed-count `constexpr` 1D and ND | Built and tested in the C++20/C++23 Linux matrix |
| Compile-time epsilon 1D | Exercised when `PF_HAS_CONSTEXPR_EPS_OVERLOAD` is true; in current CI that means the Linux GCC jobs |

### `poly_eval::pack(...)`

Packs several 1D `FuncEval` objects into a `FuncEvalMany`.

Best use case:

- several independent 1D evaluators
- same input/output scalar family
- evaluated together often enough that a grouped path matters

Not for:

- `FuncEvalND`
- a single evaluator by itself

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
- ND runtime fits are not intended to be constant-evaluated

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
- `MaxCoeffs<N>` sets the search cap
- `EvalPts<N>` sets how many validation points are checked per candidate
- `Iters<N>` sets refinement passes after the initial fit

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
- can be used at runtime or in constant evaluation in C++20

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
- for fixed-size indexable inputs
- uses template parameters for the domain bounds and coefficient count
- can be constant-evaluated when the callable is `constexpr`
- runtime-sized ND fits remain runtime-only

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

Notes:

- C++20 only
- available only when `PF_HAS_CONSTEXPR_EPS_OVERLOAD` is enabled
- currently CI-tested on the Linux GCC jobs where `PF_HAS_CONSTEXPR_EPS_OVERLOAD` is true
- 1D only
- uses template parameters, not runtime tags

Template parameters:

- `EPS`: target relative error
- `a`, `b`: domain bounds
- `MAX_NCOEFFS` (default `32`): maximum coefficient count to try
- `EVAL_POINTS` (default `100`): number of validation points per candidate
- `ITERS` (default `1`): refinement passes after the initial fit

This overload does not take runtime tags. The template parameters replace
`MaxCoeffs`, `EvalPts`, and `Iters`.

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
- `operator()(x0, x1, ..., xN)` evaluates from separate coordinates
- `operator()(pts, out, count)` evaluates a batch of canonical ND points
- `operator()(span_pts, span_out)` evaluates a batch through `std::span` when available
- `operator()(points, out)` evaluates batches stored in `data()`-backed outer containers
- `nCoeffsPerAxis()` returns the active coefficient count used on each axis

Input/output model:

- the callable's declared ND point and return types must be fixed-size containers
- accepted point types expose `size()` and `operator[]`
- the container-pair `operator()` also requires `data()` on the outer point/output containers

### `FuncEvalMany`

Packed evaluator returned by `pack(...)`.

When it helps:

- repeated grouped evaluation of several 1D approximations
- especially when the same bundle is evaluated on the same input or on batches of inputs

Useful members:

- `size()` returns the number of packed evaluators
- `nCoeffs()` returns the shared active coefficient count
- `operator()(x)` evaluates all polynomials at the same point
- `operator()(xs)` evaluates each polynomial at its corresponding point
- `operator()(first, rest...)` evaluates from separate scalar arguments
- `operator()(tuple)` evaluates from a tuple of inputs
- `operator()(pts, out, count)` evaluates many points into a flat output buffer
- `truncate(eps)` trims trailing coefficient rows

Useful constant:

- `FuncEvalMany::COUNT`

## `pack(...)` Constraints

- `pack(...)` is for 1D `FuncEval` objects
- all packed evaluators must use the same input and output types
- all packed evaluators must use the same sizing mode
- mixing runtime-sized and compile-time-sized evaluators is unsupported
- runtime-sized evaluators can be packed together only if they use the same active coefficient count

Common call forms:

- `packed(x)`: same input for every polynomial
- `packed(xs)`: one input per polynomial
- `packed(pts, out, count)`: bulk evaluation into a flat output buffer

Basic example:

```cpp
auto sinApprox = poly_eval::fit<8>([](double x) { return std::sin(x); }, -1.0, 1.0);
auto cosApprox = poly_eval::fit<8>([](double x) { return std::cos(x); }, -1.0, 1.0);
auto packed = poly_eval::pack(sinApprox, cosApprox);

auto samePoint = packed(0.5);
auto perPoly = packed(std::array<double, 2>{0.25, 0.75});
```
