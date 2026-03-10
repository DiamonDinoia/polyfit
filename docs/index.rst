polyfit
=======

**polyfit** is a fast, header-only C++17/20 library for polynomial approximation of arbitrary
functions. Given a function and an interval, it fits a polynomial using Chebyshev nodes and
Björck–Pereyra Newton interpolation, converting to monomial form for efficient evaluation via
Horner's method with optional SIMD acceleration.

Key properties:

- **No Vandermonde matrix** — uses Newton divided differences directly.
- **Compensated arithmetic** — near-machine-precision fitting via EFT in the fitting pipeline.
- **SIMD-accelerated evaluation** — bulk ``operator()`` uses xsimd for vectorized Horner.
- **C++17 and C++20** — constexpr fitting paths available in C++20.
- **Header-only** — drop into your project with CMake FetchContent or CPM.

Algorithm overview
------------------

1. Sample ``f`` on ``nCoeffs`` Chebyshev nodes on ``[a, b]``.
2. Compute Newton divided differences via Björck–Pereyra (with EFT compensation).
3. Convert Newton form to monomial basis via ``newton_to_monomial`` (with EFT compensation).
4. Optionally apply iterative refinement to correct residual error.
5. Evaluate with Horner's method; SIMD batch evaluation for bulk points.

Quick start
-----------

**Via CMake (recommended):**

.. code-block:: cmake

   include(FetchContent)
   FetchContent_Declare(polyfit
     GIT_REPOSITORY https://github.com/DiamonDinoia/polyfit.git
     GIT_TAG        main
   )
   FetchContent_MakeAvailable(polyfit)
   target_link_libraries(my_target PRIVATE polyfit::polyfit)

**Public headers** (see :doc:`install` for integration details):

.. code-block:: cpp

   #include <polyfit/polyfit.hpp>

   // Runtime coefficient count
   auto eval = poly_eval::make_func_eval(
       [](double x){ return std::sin(x); }, 16, 0.0, 3.14159);
   double y = eval(1.5);

   // Compile-time coefficient count
   auto eval2 = poly_eval::make_func_eval<16>(
       [](double x){ return std::exp(x); }, 0.0, 1.0);

Benchmark results
-----------------

Charts from the latest ``main`` branch benchmarks (GCC 14, GCC 15, LLVM 20, LLVM 21 on AVX2):

.. image:: https://raw.githubusercontent.com/DiamonDinoia/polyfit/benchmark-results/horner_performance.svg
   :alt: Horner evaluation performance

.. image:: https://raw.githubusercontent.com/DiamonDinoia/polyfit/benchmark-results/fitting_performance.svg
   :alt: 1D fitting performance

.. image:: https://raw.githubusercontent.com/DiamonDinoia/polyfit/benchmark-results/cross_compiler_overview.svg
   :alt: Cross-compiler overview

.. toctree::
   :maxdepth: 2
   :caption: Contents

   install
   guides/make_func_eval
   guides/func_eval_many
   guides/accuracy
   api/library_root
