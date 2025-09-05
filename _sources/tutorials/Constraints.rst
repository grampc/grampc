.. _sec:Tut:Constraints:

Constraint Tuning
-----------------

This section shows how the options of GRAMPC can be adjusted in such a
way that the performance in terms of computation time and optimality is
improved. To this end, a specific OCP problem is considered. However,
the approach can serve as template for different problems.

.. _problem-formulation-4:

Problem formulation
~~~~~~~~~~~~~~~~~~~

The system at hand is a double integrator with one inequality constraint
and a terminal equality constraint for each state
:footcite:`BrysonHo:book:1969`. The OCP problem is then defined as

.. math::
    :label: eq:BrysonExample

   	\min_{u} \quad & J(\mb{x}, u) = \int_0^T r u^2 \, {\rm d}t

   	\textrm{s.t.} \quad & {\dot x_1}(t) = x_2(t) \,, \qquad x_1(0) = x_{1,0}

    & {\dot x_2}(t) = u(t) \,, \phantom{_2} \qquad x_2(0) = x_{2,0}

    & h(\mb{x}(t)) = x_1(t) - 0.1 \leq 0

    & g_{\mathrm{T},1}(\mb{x}(T)) = x_1(T) = 0

    & g_{\mathrm{T},2}(\mb{x}(T)) = x_2(T) + 1 = 0,

with the weight :math:`r= 0.5` and the initial state
:math:`\mb{x}(0) = [0, 1]^\mathsf{T}`. The target of the problem
is to steer the double integrator states to the terminal state
:math:`\mb{x}(T) = [0, -1]^\mathsf{T}` without violating the
inequality constraint :math:`h(\mb{x}(t))` . 
The time in which the set point change should be executed is set to
:math:`T=1\,\mathrm{s}`.

Tuning approach
~~~~~~~~~~~~~~~

.. figure:: ../img/tikz/BrysonIntegrator1.*
    :name: fig:BrysonIntegrator1

    The prediction plot of GRAMPC for the double integrator example in :math:numref:`eq:BrysonExample`.

In the listing below, the step by step approach of tuning the parameters
of the augmented Lagrangian algorithm are detailed. The computation time
as well as the number of outer and inner iterations are shown in
:numref:`tab:TuningCompTime` along with the corresponding option
that was added or changed in each step.

.. list-table::  Computation time and outer / inner iteration count for the different settings. The last column shows which additional parameter was set different from the initial values in each step.
    :name: tab:TuningCompTime
    :widths: auto
    :header-rows: 1

    * - Step
      - Time
      - MultIter
      - GradIter (mean)
      - Additional option
    * - 0.
      - 343 ms
      - 3883
      - 10.3
      - Default settings
    * - 1.
      - 129 ms
      - 637
      - 22.4
      - ``grampc_estim_penmin``
    * - 2.
      - 59 ms
      - 171
      - 42.0
      - :math:`\text{PenaltyIncreaseFactor}=1.5`
    * - 3.
      - 36 ms
      - 98
      - 46.6
      - :math:`\text{PenaltyIncreaseThreshold}=0.75`
    * - 4.
      - 22 ms
      - 42
      - 61.2
      - :math:`\text{LineSearchInit}=5\text{e-7}`
    * - 4.
      - 22 ms
      - 42
      - 61.2
      - :math:`\text{PenaltyMin}=2\text{e4}`

#. The choice of the initial penalty parameter is crucial for numerical
   conditioning and therefore convergence. A value that is too high will
   initially put an unnecessary amount of weight on constraint
   satisfaction and mostly ignore the optimality. If the value is too
   low, the opposite will happen, i.e. at first the cost function is
   decreased at the cost of constraint violation. This will be
   especially detrimental in MPC applications. If no prior knowledge is
   available, it is recommended to use the function
   ``grampc_estim_penmin`` (or the corresponding Cmex interface ``grampc_estim_penmin_Cmex``). 
   For the example at hand, this reduces
   the computation time by approximately a third.
   :numref:`fig:BrysonIntegrator1` shows the simulation results
   using the estimated value for ``PenaltyMin``.

#. While the convergence speed is significantly increased, the penalty
   parameters at :math:`0.3\,\mathrm{s}` and :math:`0.7\,\mathrm{s}`
   are several magnitudes greater than the initial value. Since this
   huge penalty parameter occurs only at two points during the
   simulation, it is advisable to set the ``PenaltyIncreaseFactor`` to a bigger value. This again
   reduces the computation time by more than half, since fewer increases
   of the penalty parameter during the outer iterations are necessary.
   The value should not be chosen too big, as this will have an adverse
   effect on the numerical conditioning and convergence.

#. In accordance with the previous step, the threshold to increase the
   penalty parameter, i.e. ``PenaltyIncreaseThreshold`` is lowered in order to increase the penalty
   parameter more aggressively. This step almost doubles the convergence
   speed. Note that this step and the previous step are interchangeable.

#. Another common tuning possibility is the initial step size of the
   line search, especially if one of the explicit methods is used,
   cf. :ref:`sec:AlgOpt:LineSearchExplicit`. Note that the
   initial value ``LineSearchInit`` is used in the case that the explicit formula results
   in a negative step size. One approach is to use the ``LineSearchExpAutoFallback`` option. However,
   problem specific tuning (mostly trial and error) can result in a
   significant performance boost. In the example at hand, this results
   in approximately 40% faster convergence speed.

#. To further optimize the parameters, the estimation function for the
   minimal penalty parameter can be deactivated again and a better value
   for ``PenaltyMin`` be used (note that the parameter is increased until there is no
   further improvement or an decrease in performance). This results in
   an additional 15%` decrease of computation time.

.. footbibliography::