# polyfit API

## Main Entry Points

### `poly_eval::fit(...)`

Builds a polynomial evaluator from a callable and a domain.

Supported forms:

- `poly_eval::fit(f, nCoeffs, a, b, ...)`
- `poly_eval::fit(f, errorTarget, a, b, ...)`
- `poly_eval::fit<NCOEFFS>(f, a, b, ...)`
- `poly_eval::fit<eps, a, b>(f)` in the supported constexpr GCC path

The returned type is:

- `FuncEval` for 1D scalar or complex input/output
- `FuncEvalND` for `std::array`-like inputs and outputs

### `poly_eval::pack(...)`

Packs several 1D `FuncEval` objects into a `FuncEvalMany`.

All packed evaluators must use the same sizing mode:

- all compile-time coefficient counts, or
- all runtime coefficient counts

Mixing runtime-sized and compile-time-sized evaluators is unsupported.

## Main Types

### `FuncEval`

1D evaluator returned by `fit(...)`.

Key members:

- `operator()(x)` evaluates one point
- `operator()(pts, out, count)` evaluates many points
- `coeffs()` returns Horner-order coefficients
- `nCoeffs()` returns the current coefficient count
- `truncate(eps)` removes leading small coefficients for runtime-sized fits

### `FuncEvalND`

N-dimensional evaluator returned when the input type is `std::array`-like.

Key member:

- `nCoeffsPerAxis()` returns the coefficient count used on each axis

### `FuncEvalMany`

Packed evaluator returned by `pack(...)`.

Key members:

- `operator()(x)` evaluates all packed polynomials at one point
- `operator()(xs)` evaluates each packed polynomial at its corresponding point
- `operator()(pts, out, count)` evaluates many points
- `nCoeffs()` returns the shared active coefficient count
- `truncate(eps)` trims trailing coefficient rows

## 1D Fit Options

Tags are order-independent and apply to the 1D `fit(...)` overloads.

### `poly_eval::Iters<N>`

Refinement passes after the initial fit.

Default: `Iters<1>`

Use larger values for higher-degree fits if residual accuracy matters.

### `poly_eval::MaxCoeffs<N>`

Upper bound for the adaptive `fit(f, eps, a, b, ...)` overload.

Default: `MaxCoeffs<32>`

If no fit satisfies the requested error before `N`, the overload throws.

### `poly_eval::EvalPts<N>`

Number of evaluation points used by the adaptive overload when checking the error target.

Default: `EvalPts<100>`

### `poly_eval::FuseAuto`, `poly_eval::FuseAlways`, `poly_eval::FuseNever`

Control whether the domain mapping is baked into the coefficients.

- `FuseAuto` picks a safe default
- `FuseAlways` always fuses
- `FuseNever` keeps the mapping separate

## Examples

Fixed-size runtime fit:

```cpp
auto approx = poly_eval::fit([](double x) { return std::cos(x); }, 12, -1.0, 1.0);
```

Adaptive fit:

```cpp
auto approx = poly_eval::fit(
    [](double x) { return std::sin(x); },
    1e-12,
    -1.0,
    1.0,
    poly_eval::MaxCoeffs<48>{},
    poly_eval::EvalPts<200>{},
    poly_eval::Iters<2>{});
```

Compile-time fit:

```cpp
constexpr auto approx = poly_eval::fit<8>([](double x) { return x * x + 1.0; }, -1.0, 1.0);
```

Packed evaluation:

```cpp
auto sinApprox = poly_eval::fit<8>([](double x) { return std::sin(x); }, -1.0, 1.0);
auto cosApprox = poly_eval::fit<8>([](double x) { return std::cos(x); }, -1.0, 1.0);
auto packed = poly_eval::pack(sinApprox, cosApprox);
```
