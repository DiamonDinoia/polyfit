# polyfit API Reference

This document summarizes the public API of polyfit and the responsibilities of the main types.

## Primary entry points

- make_func_eval(...)
  - Overloads support:
    - Compile-time coefficient count: make_func_eval<N>(func, a, b)
    - Runtime fixed coefficient count: make_func_eval(func, nCoeffs, a, b)
    - Runtime adaptive (epsilon): make_func_eval(func, eps, a, b)
    - ND / vector-valued functions use std::array-like inputs/outputs
  - Behavior:
    - Samples f on Chebyshev nodes, fits via Björck–Pereyra Newton interpolation → monomial conversion. No Vandermonde matrix is formed.
    - Returns a callable evaluator (FuncEval or FuncEvalND) that evaluates with Horner's method.
    - For C++20 constexpr overloads, fitting may occur at compile time.
  - Key template parameters on the returned `FuncEval`:
    - compile-time coefficient count: number of monomial coefficients known at compile time (`0` = runtime coefficient count)
    - refinement iterations: number of iterative refinement passes (default `1`); increase to `2-3` when residual error matters
  - Runtime `eps` overload:
    - Searches coefficient counts from `1` to `MaxN` and returns the first evaluator whose checked error is within `eps`.

## Main types

- poly_eval::FuncEval<Func, NCoeffsCt = 0, RefineIters = 1>
  - 1D polynomial evaluator, supports:
    - `NCoeffsCt > 0`: fixed compile-time coefficient count
    - `NCoeffsCt == 0`: runtime-sized coefficient storage (`nCoeffs` passed to the constructor)
  - Key members:
    - operator()(InputType) — scalar evaluation
    - operator()(const InputType* pts, OutputType* out, std::size_t num) — bulk SIMD evaluation
    - coeffs() — access monomial coefficients (monomial basis, highest-order term first)
    - nCoeffs() — current coefficient count
    - truncate(eps) — remove leading coefficients whose absolute value is ≤ eps, reducing evaluation work

- poly_eval::FuncEvalND<Func, NCoeffsCt = 0>
  - N‑D variant using mdspan extents for storage.
  - Use for functions taking tuple-like / std::array inputs.
  - `nCoeffsPerAxis()` reports the coefficient count used along each axis.

- poly_eval::FuncEvalMany<EvalTypes...>
  - Packs several FuncEval instances for SIMD-friendly batched evaluation.
  - All packed evaluators must use the same coefficient-count mode:
    - all compile-time coefficient count, or
    - all runtime coefficient count.
    Mixing runtime and compile-time coefficient count evaluators in one pack is unsupported.
  - Key members:
    - operator()(InputType x) — evaluates all packed polynomials at x; returns an array of outputs
    - operator()(InputType first, Ts... rest) — variadic overload; evaluates each packed polynomial at its corresponding argument
    - operator()(const std::tuple<Ts...>&) — tuple overload; elements of the tuple are forwarded to the corresponding polynomial
    - operator()(const InputType* pts, OutputType* out, std::size_t n) — bulk evaluation across many points using SIMD
    - nCoeffs() — current shared coefficient count
    - truncate(eps) — delegates to each packed evaluator's truncate(eps)

## Config / macros

Macros used by the implementation are internal. See include/polyfit/internal/macros.h for internal definitions (these are not part of the public API):

- PF_ALWAYS_INLINE, PF_NO_INLINE, PF_RESTRICT, PF_C20CONSTEXPR, PF_C23STATIC, PF_ASSUME

## Mapping semantics

- Domain mapping:
  - Internal representation maps [a,b] to the canonical Chebyshev domain used for sampling.
  - Evaluator maps back to original domain when evaluating points.

## Coefficient Order

- 1D evaluators and helpers expect monomial coefficients in Horner order (highest-order term first).
- Internal fitting converts Newton-form coefficients to monomial form and reverses the order before evaluation.
- When interacting with low-level APIs such as `horner`, pass coefficients as `[c_N, c_{N-1}, …, c_0]`.

## Accuracy

### Compensated Fitting Pipeline

- **Björck–Pereyra** (Newton divided differences) and **newton_to_monomial** (monomial conversion) use compensated arithmetic for arithmetic and complex coefficient types.
- **Compensated Horner** is used in the iterative refinement step to prevent divergence when evaluating large fits during residual correction.
- Compensation is selective, not universal: compile-time fitting uses reduced compensation where constexpr restrictions apply, and unsupported coefficient types fall back to ordinary arithmetic.

### Truncation

- `FuncEval::truncate(eps)` removes leading coefficients with absolute value ≤ eps.
- The adaptive (`epsilon`) overload finds the first coefficient count that satisfies the error bound; it does not rely on automatic truncation afterward.

## MSVC Compatibility

- **C++17**: fully supported. Internal constants use `detail::constants::pi` instead of `M_PI`; `detail::cos()` uses a Cody-Waite minimax approximation (consistent across standards) to produce numerically stable Chebyshev nodes at high polynomial degrees.
- **C++20 constexpr fitting**: requires `/Zc:__cplusplus` compiler flag so that MSVC correctly reports the language standard version to the library's feature-detection macros. Without this flag, constexpr fitting silently falls back to runtime fitting.

## C++ Standards

- C++17: fully supported for runtime fitting and evaluation (tests build as `_cxx17` targets).
- C++20: enables additional constexpr paths and cleaner syntax (tests build as `_cxx20` targets).
- The public headers avoid requiring C++20 features where not necessary; internal implementations use feature-detection macros.

## Running Tests and Benchmarks

- Configure with tests enabled:
  - `cmake -S . -B build -DPOLYFIT_BUILD_TESTS=ON`
- Build and run tests (creates both C++17 and C++20 variants):
  - `cmake --build build`
  - `ctest --test-dir build --output-on-failure`
- Test targets are generated per standard, e.g. `test_horner_cxx17`, `test_horner_cxx20`.
- Benchmark targets mirror the matrix as `bench_*_cxx17` and `bench_*_cxx20`.

## Notes & recommendations

- Use compile-time overloads when possible to generate the fastest evaluators.
- For adaptive fitting, provide a sensible `maxNCoeffs` limit and evaluation grid size when customizing.
- For multidimensional functions, prefer std::array-based inputs/outputs for clarity.

## Examples

See examples/ for runnable examples demonstrating the overloads:

- examples/runtime_1D_fixed.cpp
- examples/runtime_1D_eps.cpp
- examples/runtime_ND_fixed.cpp

### FuncEvalMany (packing example)

`FuncEvalMany` packs multiple `FuncEval` instances so that coefficients are
laid out for efficient SIMD/batched evaluation. This is useful when you want
to evaluate several related polynomials for the same input values.

Minimal usage:

```cpp
auto f1 = make_func_eval<8>([](double x){ return std::sin(x); }, 0.0, 1.0);
auto f2 = make_func_eval<8>([](double x){ return std::cos(x); }, 0.0, 1.0);
auto packed = make_func_eval_many(f1, f2);
double out0 = packed(0.5)[0]; // value for first polynomial
double out1 = packed(0.5)[1]; // value for second polynomial
```

Padding and SIMD notes:

- If the number of packed polynomials is not a multiple of the SIMD width,
  the implementation pads coefficients to the next vector width (`kFPad`).
- Use the bulk `operator()(const InputType *x, OutputType *out, std::size_t)`
  overload to drive efficient vectorized evaluation across many points.
