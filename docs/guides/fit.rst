fit
===

``poly_eval::fit(...)`` builds the library's evaluator types.

This guide is for choosing the right overload quickly. For the complete overload
and tag reference, see the API reference page.

Returned type
-------------

- 1D scalar or complex callables return ``FuncEval``
- ``std::array``-like inputs return ``FuncEvalND``

Quick overload guide
--------------------

===============================  ==========  =======  ==========================
Form                             Dimensions  Runtime  Constant evaluation
===============================  ==========  =======  ==========================
``fit(f, nCoeffs, a, b, ...)``   1D          Yes      No
``fit(f, eps, a, b, ...)``       1D          Yes      No
``fit<NCOEFFS>(f, a, b, ...)``   1D          Yes      Yes, in C++20
``fit(f, nCoeffs, a, b)``        ND          Yes      No
``fit<NCOEFFS, a, b>(f)``        ND          Yes      Yes, in C++20
``fit<EPS, a, b, ...>(f)``       1D          No       Yes, when enabled
===============================  ==========  =======  ==========================

Runtime 1D fixed coefficient count
----------------------------------

.. code-block:: cpp

   auto approx = poly_eval::fit(
       [](double x) { return std::cos(x); },
       12,
       -1.0,
       1.0);

Runtime 1D adaptive fit
-----------------------

.. code-block:: cpp

   auto approx = poly_eval::fit(
       [](double x) { return std::sin(x); },
       1e-12,
       -1.0,
       1.0,
       poly_eval::MaxCoeffs<48>{},
       poly_eval::EvalPts<200>{},
       poly_eval::Iters<2>{},
       poly_eval::FuseNever{});

The adaptive overload searches upward from ``1`` coefficient and throws if the
requested error is not met before ``MaxCoeffs``.

Tag meanings for the adaptive overload:

- ``MaxCoeffs<N>``: maximum coefficient count to try
- ``EvalPts<N>``: validation point count for each candidate polynomial
- ``Iters<N>``: refinement passes after the initial fit
- ``FuseAuto``, ``FuseAlways``, ``FuseNever``: control whether the domain map is folded into the coefficients

Compile-time 1D fixed coefficient count
---------------------------------------

.. code-block:: cpp

   constexpr auto approx = poly_eval::fit<8>(
       [](double x) { return x * x + 1.0; },
       -1.0,
       1.0,
       poly_eval::Iters<2>{});

This overload is still a 1D fit. The coefficient count is fixed at compile
time, but the evaluator can still be used at runtime.

Runtime ND fixed coefficient count
----------------------------------

.. code-block:: cpp

   auto approx = poly_eval::fit([](const std::array<double, 2> &p) {
       return std::array<double, 2>{std::cos(p[0]) + std::sin(p[1]), p[0] * p[1]};
   }, 10, {-1.0, -1.0}, {1.0, 1.0});

   auto y = approx({0.25, -0.5});

This is the ND runtime path. It uses a runtime coefficient count per axis and
does not use the 1D tuning tags.

ND fit with template-parameter bounds
-------------------------------------

.. code-block:: cpp

   constexpr std::array<double, 2> a{-1.0, -1.0};
   constexpr std::array<double, 2> b{1.0, 1.0};

   auto approx = poly_eval::fit<8, a, b>([](const std::array<double, 2> &p) {
       return std::array<double, 2>{p[0] + p[1], p[0] * p[1]};
   });

This overload is C++20-only. It can be constant-evaluated when the callable is
``constexpr``. Runtime-sized ND fits remain runtime-only.

Compile-time epsilon fit
------------------------

.. code-block:: cpp

   constexpr auto approx = poly_eval::fit<1e-12, -1.0, 1.0>([](double x) constexpr {
       return std::sin(x) + x * x;
   });

This overload:

- is C++20 only
- is enabled only when ``PF_HAS_CONSTEXPR_EPS_OVERLOAD`` is true
- is currently CI-tested on the Linux GCC jobs where ``PF_HAS_CONSTEXPR_EPS_OVERLOAD`` is true
- is 1D only
- uses template parameters instead of runtime tags

Full form:

.. code-block:: cpp

   constexpr auto approx = poly_eval::fit<1e-12, -1.0, 1.0, 48, 200, 2>([](double x) constexpr {
       return std::sin(x) + x * x;
   });

Template parameters:

- ``EPS``: target relative error
- ``a``, ``b``: domain bounds
- ``MAX_NCOEFFS`` (default ``32``): maximum coefficient count to try
- ``EVAL_POINTS`` (default ``100``): number of validation points per candidate
- ``ITERS`` (default ``1``): refinement passes after the initial fit

Those last three template parameters serve the same role as ``MaxCoeffs``,
``EvalPts``, and ``Iters`` in the runtime adaptive API.

Tags
----

The 1D runtime fixed-count, runtime adaptive, and compile-time fixed-count
overloads accept optional tuning tags. See the API reference for the full tag
list and defaults.

Rules:

- tags are optional
- tags can appear in any order
- duplicate tag kinds are a compile-time error

Applicability:

- ``MaxCoeffs`` and ``EvalPts`` only affect adaptive ``fit(f, eps, a, b, ...)``
- ``Iters`` and fusion tags affect the 1D fixed-count overloads
- runtime ND fits ignore these 1D-only tuning tags
- compile-time epsilon fits use template parameters, not tags
