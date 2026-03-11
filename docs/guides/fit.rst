fit
===

``poly_eval::fit(...)`` is the main factory for polynomial evaluators.

1D fixed coefficient count
--------------------------

.. code-block:: cpp

   auto approx = poly_eval::fit(
       [](double x) { return std::cos(x); },
       12,
       -1.0,
       1.0);

1D adaptive fit
---------------

.. code-block:: cpp

   auto approx = poly_eval::fit(
       [](double x) { return std::sin(x); },
       1e-12,
       -1.0,
       1.0,
       poly_eval::MaxCoeffs<48>{},
       poly_eval::EvalPts<200>{},
       poly_eval::Iters<2>{});

The adaptive overload searches from 1 coefficient upward and throws if the requested error is not met before ``MaxCoeffs``.

Compile-time coefficient count
------------------------------

.. code-block:: cpp

   constexpr auto approx = poly_eval::fit<8>(
       [](double x) { return x * x + 1.0; },
       -1.0,
       1.0);

N-D fit
-------

.. code-block:: cpp

   using In = std::array<double, 2>;
   using Out = std::array<double, 2>;

   auto approx = poly_eval::fit(
       [](const In& p) {
           return Out{std::cos(p[0]) + std::sin(p[1]), p[0] * p[1]};
       },
       10,
       In{-1.0, -1.0},
       In{1.0, 1.0});

Tags
----

The 1D overloads accept these order-independent tags:

- ``poly_eval::Iters<N>``
- ``poly_eval::MaxCoeffs<N>``
- ``poly_eval::EvalPts<N>``
- ``poly_eval::FuseAuto``
- ``poly_eval::FuseAlways``
- ``poly_eval::FuseNever``
