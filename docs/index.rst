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

.. toctree::
   :maxdepth: 2
   :caption: Contents

   install
   guides/fit
   guides/pack
   guides/accuracy
   api/library_root
