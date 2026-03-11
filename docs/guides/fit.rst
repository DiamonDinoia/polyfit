fit
===

``poly_eval::fit(...)`` builds the library's evaluator types.

Returned type
-------------

- 1D scalar or complex callables return ``FuncEval``
- ``std::array``-like inputs return ``FuncEvalND``

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

Compile-time 1D fixed coefficient count
---------------------------------------

.. code-block:: cpp

   constexpr auto approx = poly_eval::fit<8>(
       [](double x) { return x * x + 1.0; },
       -1.0,
       1.0,
       poly_eval::Iters<2>{});

Runtime ND fixed coefficient count
----------------------------------

.. code-block:: cpp

   auto approx = poly_eval::fit([](const std::array<double, 2> &p) {
       return std::array<double, 2>{std::cos(p[0]) + std::sin(p[1]), p[0] * p[1]};
   }, 10, {-1.0, -1.0}, {1.0, 1.0});

   auto y = approx({0.25, -0.5});

ND fit with template-parameter bounds
-------------------------------------

.. code-block:: cpp

   constexpr std::array<double, 2> a{-1.0, -1.0};
   constexpr std::array<double, 2> b{1.0, 1.0};

   auto approx = poly_eval::fit<8, a, b>([](const std::array<double, 2> &p) {
       return std::array<double, 2>{p[0] + p[1], p[0] * p[1]};
   });

Compile-time epsilon fit
------------------------

.. code-block:: cpp

   constexpr auto approx = poly_eval::fit<1e-12, -1.0, 1.0>([](double x) constexpr {
       return std::sin(x) + x * x;
   });

This overload:

- is C++20 only
- is enabled only when ``PF_HAS_CONSTEXPR_EPS_OVERLOAD`` is true
- is currently intended for the supported GCC constexpr path
- is 1D only
- uses template parameters instead of runtime tags

Full form:

.. code-block:: cpp

   constexpr auto approx = poly_eval::fit<1e-12, -1.0, 1.0, 48, 200, 2>([](double x) constexpr {
       return std::sin(x) + x * x;
   });

Tags
----

The 1D runtime fixed-count, runtime adaptive, and compile-time fixed-count
overloads accept these tags:

- ``poly_eval::Iters<N>``
- ``poly_eval::MaxCoeffs<N>``
- ``poly_eval::EvalPts<N>``
- ``poly_eval::FuseAuto``
- ``poly_eval::FuseAlways``
- ``poly_eval::FuseNever``

Rules:

- tags are optional
- tags can appear in any order
- duplicate tag kinds are a compile-time error

Defaults:

- ``Iters<1>``
- ``MaxCoeffs<32>``
- ``EvalPts<100>``
- ``FuseAuto{}``

Applicability:

- ``MaxCoeffs`` and ``EvalPts`` only affect adaptive ``fit(f, eps, a, b, ...)``
- ``Iters`` and fusion tags affect the 1D fixed-count overloads
- runtime ND fits ignore these 1D-only tuning tags
- compile-time epsilon fits use template parameters, not tags
