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

## Configuration tags

All `make_func_eval` overloads accept an optional, order-independent pack of
tag arguments after the domain endpoints. Tags customise fitting behaviour
without changing the function signature.

```cpp
// Any subset, any order:
auto poly = poly_eval::make_func_eval(f, nCoeffs, a, b,
    poly_eval::iters<2>{},
    poly_eval::fuse_never{});
```

### `iters<N>` — iterative refinement passes

**Default:** `iters<1>`

Each pass re-evaluates the current polynomial at the Chebyshev nodes using
compensated Horner, computes residuals `f(xᵢ) − p(xᵢ)`, fits a correction
polynomial to those residuals, and adds it to the coefficients. This converges
the fit toward machine precision when the initial Björck–Pereyra solve leaves a
non-trivial residual (typically at high polynomial degree or near-degenerate
domains).

For `nCoeffs > 32` two extra passes are added automatically regardless of `N`,
because compensated residuals are needed to prevent divergence at high degree.

| Coefficient count | Recommended `N` |
|---|---|
| ≤ 40 | 1 (default) |
| 40 – 80 | 2 |
| 80 – 128 | 3 |

```cpp
// High-degree fit: more refinement passes
auto poly = poly_eval::make_func_eval(f, 64, -1.0, 1.0,
    poly_eval::iters<3>{});
```

---

### `maxNCoeffs<N>` — search limit for the epsilon overload

**Default:** `maxNCoeffs<32>`
**Applies to:** `make_func_eval(f, eps, a, b, ...)` only.

The epsilon overload tries coefficient counts 1, 2, … up to `N` and returns
the first evaluator whose maximum error on the evaluation grid is within `eps`.
If no count satisfies the bound, it returns the evaluator with the smallest
observed error.

Increase `N` when the default ceiling of 32 is not enough to reach the target
accuracy. Decrease it to cap runtime cost when an approximate fit is acceptable.

```cpp
// Allow up to 64 coefficients to hit 1e-14 error
auto poly = poly_eval::make_func_eval(f, 1e-14, -1.0, 1.0,
    poly_eval::maxNCoeffs<64>{});
```

---

### `evalPts<N>` — error-checking grid size for the epsilon overload

**Default:** `evalPts<100>`
**Applies to:** `make_func_eval(f, eps, a, b, ...)` only.

The epsilon overload checks the fitted polynomial against `N` equispaced points
in `[a, b]` to estimate the approximation error. The default of 100 is
sufficient for smooth functions. Increase it for functions with sharp features
that might be missed by a coarse grid.

```cpp
// Dense grid to catch a near-discontinuity
auto poly = poly_eval::make_func_eval(f, 1e-8, -1.0, 1.0,
    poly_eval::evalPts<500>{},
    poly_eval::maxNCoeffs<48>{});
```

---

### `fuse_auto` / `fuse_always` / `fuse_never` — domain fusion

**Default:** `fuse_auto`

Domain fusion bakes the linear mapping `[a, b] → [−1, 1]` directly into the
polynomial coefficients at construction time. The evaluator stores the
transformed polynomial `q(x) = p(α·x + β)` so that evaluating at a point `x`
requires no per-point division or affine operation — just the Horner recurrence
on `x` directly.

| Tag | Behaviour |
|---|---|
| `fuse_auto` | Fuse when the domain's condition number is small enough to avoid coefficient blow-up (see below). |
| `fuse_always` | Always fuse, regardless of numerical properties. |
| `fuse_never` | Never fuse; apply the domain mapping at every evaluation point. |

**When `fuse_auto` decides to fuse** — it checks:

```
cond = |α| + |β| + 1   (where α = 2/(b−a), β = −(b+a)/(b−a))
fuse  iff  (nCoeffs − 1) · log₁₀(cond) < digits₁₀(T) − 3
```

For double this threshold is roughly `log₁₀(cond) < 12 / (nCoeffs − 1)`.
Narrow domains near `[−1, 1]` almost always fuse; wide or heavily offset
domains (e.g. `[1000, 2000]`) may not fuse at high degree.

**When to override:**

Use `fuse_never` when:
- The domain is wide or offset and you intend to call `truncate(eps)` after
  fitting — fusion changes coefficient magnitudes and makes the truncation
  threshold meaningless relative to the original polynomial scale.
- You need reproducible coefficients independent of domain scaling.

Use `fuse_always` when:
- The domain is narrow (e.g. `[−0.1, 0.1]`) and you want to guarantee zero
  per-point overhead even for small `nCoeffs` where `auto` might not fuse.

```cpp
// Wide offset domain: let auto decide (will likely skip fusion)
auto poly1 = poly_eval::make_func_eval(f, 12, 1000.0, 2000.0);

// Force no fusion — preserve coefficient scale for truncation
auto poly2 = poly_eval::make_func_eval(f, 12, 1000.0, 2000.0,
    poly_eval::fuse_never{});
poly2.truncate(1e-10);

// Force fusion on a narrow domain where auto might hedge
auto poly3 = poly_eval::make_func_eval<16>(f, -0.1, 0.1,
    poly_eval::fuse_always{});
```

---

### Combining tags

All tags are order-independent and composable:

```cpp
auto poly = poly_eval::make_func_eval(f, 1e-12, -1.0, 1.0,
    poly_eval::maxNCoeffs<64>{},
    poly_eval::evalPts<200>{},
    poly_eval::iters<2>{},
    poly_eval::fuse_never{});
```

Tags also work with the compile-time overload:

```cpp
auto poly = poly_eval::make_func_eval<32>(f, -1.0, 1.0,
    poly_eval::iters<2>{},
    poly_eval::fuse_always{});
```

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
