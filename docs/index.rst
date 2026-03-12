polyfit
=======

**polyfit** fits polynomial approximations on bounded domains and evaluates them efficiently.

Quick start
-----------

.. code-block:: cpp

   #include <polyfit/polyfit.hpp>

   auto approx = poly_eval::fit(
       [](double x) { return std::cos(x); },
       12,
       -1.0,
       1.0);

   double y = approx(0.5);

Core API
--------

- ``poly_eval::fit(...)`` builds ``FuncEval`` or ``FuncEvalND``.
- ``poly_eval::pack(...)`` combines several 1D evaluators into ``FuncEvalMany``.
- 1D fitting options use tags such as ``Iters<N>``, ``MaxCoeffs<N>``, ``EvalPts<N>``, and ``FuseAuto``.

Quick support guide
-------------------

===============================  ==========  =======  ==========================
Form                             Dimensions  Runtime  Constant evaluation
===============================  ==========  =======  ==========================
``fit(f, nCoeffs, a, b, ...)``   1D          Yes      No
``fit(f, eps, a, b, ...)``       1D          Yes      No
``fit<NCOEFFS>(f, a, b, ...)``   1D          Yes      Yes, in C++20
``fit(f, nCoeffs, a, b)``        ND          Yes      No
``fit<NCOEFFS, a, b>(f)``        ND          Yes      Yes, in C++20
``fit<EPS, a, b, ...>(f)``       1D          No       Yes, when enabled
``pack(e1, e2, ...)``            1D bundle   Yes      See ``guides/pack``
===============================  ==========  =======  ==========================

If you are choosing between overloads, start with :doc:`guides/fit`. If you
already have several 1D evaluators and want to evaluate them together, go to
:doc:`guides/pack`.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   install
   guides/fit
   guides/pack
   guides/accuracy
   api/library_root
