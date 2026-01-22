# polyfit API Reference

This document summarizes the public API of polyfit and the responsibilities of the main types.

## Primary entry points

- make_func_eval(...)
  - Overloads support:
    - Compile-time degree: make_func_eval<N>(func, a, b)
    - Runtime fixed degree: make_func_eval(func, n, a, b)
    - Runtime adaptive (epsilon): make_func_eval(func, eps, a, b)
    - ND / vector-valued functions use std::array-like inputs/outputs
  - Behavior:
    - Samples f on Chebyshev nodes, fits via Newton→monomial conversion.
    - Returns a callable evaluator (FuncEval or FuncEvalND) that evaluates with Horner's method.
    - For C++20 constexpr overloads, fitting may occur at compile time.

## Main types

- poly_eval::FuncEval<Func, N_compile_time = 0, Iters_compile_time = 1>
  - 1D polynomial evaluator, supports:
    - N_compile_time > 0: fixed compile-time degree
    - N_compile_time == 0: runtime degree (n passed to constructor)
  - Key members:
    - operator()(InputType) — scalar evaluation
    - operator()(const InputType* pts, OutputType* out, std::size_t num) — bulk eval
    - coeffs() — access monomial coefficients (monomial basis)

- poly_eval::FuncEvalND<Func, N_compile_time = 0>
  - N‑D variant using mdspan extents for storage.
  - Use for functions taking tuple-like / std::array inputs.

- poly_eval::FuncEvalMany<EvalTypes...>
  - Packs several FuncEval instances for SIMD-friendly batched evaluation.

## Config / macros

Macros used by the implementation are internal. See include/polyfit/internal/macros.h for internal definitions (these are not part of the public API):
- PF_ALWAYS_INLINE, PF_NO_INLINE, PF_RESTRICT, PF_C20CONSTEXPR, PF_C23STATIC, PF_ASSUME

## Mapping semantics

- Domain mapping:
  - Internal representation maps [a,b] to the canonical Chebyshev domain used for sampling.
  - Evaluator maps back to original domain when evaluating points.

## Coefficient Order

- 1D evaluators and helpers expect monomial coefficients in reversed order for Horner's method (highest degree first).
- Internal fitting converts Newton-form coefficients to monomial form and reverses the order before evaluation.
- When interacting with low-level APIs such as `horner`, pass coefficients as `[c_N, c_{N-1}, …, c_0]`.

## C++ Standards

- C++17: fully supported for runtime fitting and evaluation (tests build as `_cxx17` targets).
- C++20: enables additional constexpr paths and cleaner syntax (tests build as `_cxx20` targets).
- The public headers avoid requiring C++20 features where not necessary; internal implementations may use feature-detection macros.

## Running Tests and Benchmarks

- Configure with tests enabled:
  - `cmake -S . -B build -DMONOFIT_BUILD_TESTS=ON`
- Build and run tests (creates both C++17 and C++20 variants):
  - `cmake --build build`
  - `ctest --test-dir build --output-on-failure`
- Test targets are generated per standard, e.g. `test_horner_cxx17`, `test_horner_cxx20`.
- Benchmark targets mirror the matrix as `bench_*_cxx17` and `bench_*_cxx20`.

## Notes & recommendations

- Use compile-time overloads when possible to generate the fastest evaluators.
- For adaptive fitting, provide a sensible MaxN (default in code) and NumEvalPoints when customizing.
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
  the implementation pads coefficients to the next vector width (`kF_pad`).
- Use the bulk `operator()(const InputType *x, OutputType *out, std::size_t)`
  overload to drive efficient vectorized evaluation across many points.

