make_func_eval
==============

``make_func_eval`` is the primary factory function for creating polynomial approximations.
All overloads sample the function on Chebyshev nodes, fit via Björck–Pereyra Newton
interpolation, convert to monomial form, and return a callable evaluator.

Include
-------

.. code-block:: cpp

   #include <polyfit/polyfit.hpp>
   using namespace poly_eval;

Overloads
---------

Runtime fixed coefficient count
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   auto eval = make_func_eval(f, nCoeffs, a, b);
   // f: callable (double) -> double
   // nCoeffs: coefficient count (int)
   // a, b: interval endpoints (double)
   // Returns: FuncEval<Func, 0, 1>

Fits a polynomial with ``nCoeffs`` monomial coefficients on ``[a, b]``.

.. code-block:: cpp

   auto sin_approx = make_func_eval(
       [](double x){ return std::sin(x); }, 16, 0.0, 3.14159265);
   double y = sin_approx(1.5);  // ≈ sin(1.5)

Runtime adaptive (epsilon)
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   auto eval = make_func_eval(f, eps, a, b);
   // eps: double, target uniform error bound
   // Returns: FuncEval<Func, 0, 1>

Searches coefficient counts from ``1`` through ``MaxN`` and returns the first evaluator whose checked
error on the evaluation grid is ≤ ``eps``.

.. code-block:: cpp

   auto cos_approx = make_func_eval(
       [](double x){ return std::cos(x); }, 1e-10, 0.0, 6.28);
   // coefficient count chosen automatically

Compile-time coefficient count
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   auto eval = make_func_eval<N>(f, a, b);
   // N: compile-time coefficient count (int, template parameter)
   // Returns: FuncEval<Func, N, 1>

The coefficient count is encoded in the type, allowing the compiler to fully unroll Horner evaluation.

.. code-block:: cpp

   auto eval = make_func_eval<16>(
       [](double x){ return std::exp(x); }, 0.0, 1.0);

Compile-time coefficient count with custom refinement iterations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   auto eval = make_func_eval<N, Iters>(f, a, b);
   // Iters: number of iterative refinement passes (default 1)

Use ``Iters = 2`` or ``3`` for large fits (roughly 40+ coefficients) where a single refinement
pass may not fully correct residual error.

.. code-block:: cpp

   auto eval = make_func_eval<64, 2>(
       [](double x){ return std::tanh(x); }, -5.0, 5.0);

C++20 compile-time epsilon (NTTP)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   #if __cplusplus >= 202002L
   auto eval = make_func_eval<eps, a, b>(f);
   // eps, a, b: non-type template parameters (double)
   // Fitting may occur at compile time (constexpr path)
   #endif

.. code-block:: cpp

   constexpr double A = 0.0, B = 1.0;
   auto eval = make_func_eval<1e-12, A, B>(
       [](double x){ return std::log1p(x); });

N-D functions
^^^^^^^^^^^^^

For functions whose input or output is a ``std::array``, use the N-D overloads which return
a ``FuncEvalND``:

.. code-block:: cpp

   // Vector-valued: R -> R^2
   auto eval = make_func_eval(
       [](double t) -> std::array<double,2> { return {std::cos(t), std::sin(t)}; },
       32, 0.0, 6.28);

   auto [cx, cy] = eval(1.0);

RefineIters parameter
-----------------------------

The ``RefineIters`` template parameter on ``FuncEval`` controls how many iterative
refinement passes are applied after the initial fit:

- ``Iters = 1`` (default): one refinement pass, sufficient for fits up to roughly 40 coefficients.
- ``Iters = 2`` or ``3``: recommended for larger fits where residual from Newton-form
  conversion accumulates.

Refinement uses compensated Horner evaluation to avoid introducing new rounding errors
during the correction step.
