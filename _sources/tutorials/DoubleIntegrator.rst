.. _sec:Tut:DoubleIntegrator:

Optimal control of a double integrator
--------------------------------------

.. figure:: ../img/tikz/DoubleInt.*
    :name: fig:DoubleInt

    Numerical OCP solution of the double integrator problem with fixed end time.

This section describes how GRAMPC can be used to solve OCPs by
considering the example of a double integrator problem. The problem
formulation includes equality and inequality constraints. Both the
control variable :math:`u` and the end time :math:`T` serve as
optimization variables. In addition to the problem formulation of the
OCP, one focus of the following discussion will be the appropriate
convergence check of the augmented Lagrangian and gradient method, and
the optimization of the end time.

.. _problem-formulation-1:

Problem formulation
~~~~~~~~~~~~~~~~~~~

The double integrator problems reads as

.. math::
    :label: eq:OCP_doubleInt

    \min_{u, T} \quad & J(u, T;\mb{x}_0) = T + \int_0^T q_1u^2(t) \, {\rm d}t \\
    \textrm{s.t.} \quad & {\dot x_1}(t) = x_2(t) \,, \quad x_1(0) = x_{1,0} \\
    & {\dot x_2}(t) = u(t) \,, \quad x_2(0) = x_{2,0} \\
    & \mb{g}_{T}(\mb{x}(T))  = [x_1(T)\,, x_2(T)]^\mathsf{T}= \mb{0} \\
    & h(\mb{x} (t)) = x_2(t) - 0.5 \leq 0 \\
    & u(t) \in \left[u_{\min}, u_{\max}\right] \,,\quad T \in \left[T_{\min}, T_{\max}\right]

including two terminal equality constraints :math:`\mb{g}_{T}(\mb{x}(T))`
and one general inequality constraint :math:`h(\mb{x} (t))`.
The control task consists in a setpoint transition from the initial
state :math:`x_{1,0}=x_{2,0}=-1` to the origin :math:`x_1(T)=x_2(T)=0`.
The cost functional :math:`J(u, T;\mb{x}_0)` represents a
trade-off between time and energy optimality depending on the weight :math:`q_1=0.1`.

The problem is formulated in GRAMPC using the C file
``probfct_DOUBLE_INTEGRATOR``, which can be found in ``<grampc_root>/examples/Double_Integrator``. In particular,
the terminal equality constraints are formulated in GRAMPC via the functions ``gTfct`` , ``dgTdx_vec`` and ``dgTdT_vec``. 
The number of terminal equality constraints is set to ``NgT`` =2 in the function ``ocp_dim``. Similarly,
the inequality constraint is formulated by means of the functions ``hfct`` , ``dhdx_vec`` and ``dhdu_vec`` and setting to ``Nh`` =1 in ``ocp_dim``.
More details on implementing the OCP can be found in
:ref:`sec:ProblemImplementation` and in the example provided in the folder ``<grampc root>/examples/Double Integrator``.

The options ``OptimControl`` and ``OptimTime`` are activated to optimize not only the control
variable :math:`u` but also the end time :math:`T`. The lower and upper
bounds of the optimzation variables are set to :math:`u\in[-1,1]` and
:math:`T\in[1,10]` using the parameters ``umin`` , ``umax`` , ``Tmin`` and ``Tmax`` ,
cf. :ref:`chap:ProblemFormulation`. Note that the option ``ShiftControl`` is
deactivated, as the shifting of the control trajectory typically only
applies to MPC problems.

The option ``ConvergenceCheck`` is activated to terminate the augmented Lagrangian algorithm
as soon as the state constraints are fulfilled with sufficient accuracy
and the optimization variable has converged to an optimal value. To this
end, the convergence
criteria :math:numref:`eq:AlgOpt:ConvGradient` and
:math:numref:`eq:AlgOpt:ConvConstraints` are
evaluated after each gradient and augmented Lagrangian step using the
thresholds
:math:`\mb{\varepsilon_{\mb{g}_T}}=[1e-6,1e-6]^\mathsf{T}`
and :math:`\varepsilon_{\mb{h}}=1e-6`, cf. the option ``ConstraintsAbsTol`` .
The threshold for checking convergence of the optimization variable
:math:`\varepsilon_\text{rel,c}` is set to :math:`1e-9` via the
option ``ConvergenceGradientRelTol`` . Note that the activation of the convergence check using the
option ``ConvergenceCheck`` causes the gradient method to be aborted when the convergence
condition :math:`\eta^{i|j+1} \leq \varepsilon_\text{rel,c}` defined by
:math:numref:`eq:AlgOpt:ConvConstraints` is reached.

Optimization with fixed end time
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: ../img/tikz/DoubleInt_constr.*
    :name: fig:DoubleInt_constr

    Numerical OCP solution of the double integrator problem with fixed end time and state constraint.

In a first scenario, the OCP :math:numref:`eq:OCP_doubleInt` is numerically solved with the fixed end time :math:`T=4\,\mathrm{s}`,
i.e., the option ``OptimTime`` is set ``off`` , and without the inequality
constraint, i.e., the option ``InequalityConstraints`` is also set ``off`` . The time interval :math:`[0,T]` is
discretized using ``Nhor=50`` discretization points. During the
augmented Lagrangian iterations, a decrease of the penalty parameters is
prevented by setting the adaptation factor :math:`\beta_\text{de}=1`
(see option ``PenaltyDecreaseFactor``). The increase factor :math:`\beta_\text{in}` of the
penalty parameter update :math:numref:`eq:AlgOpt:UpdatePenh` is set to the value
:math:`1.25`. These settings ensure a fast convergence, since the
very low tolerances :math:`\mb{\varepsilon_{\mb{g}_T}}`
require high penalty parameters. However, starting with high penalties
can lead to instabilities, as high penalties distort the optimization
problem.

As shown in :numref:`fig:DoubleInt`, the state variables
:math:`\mb{x}(t)` are transferred to the origin as specified by
the terminal state constraints.
Note that only 5 augmented Lagrangian steps are required in total for
numerically solving the state constrained optimization problem in
accordance with the formulated convergence criterion for the terminal
equality constraints and the optimization variable :math:`u`. The number
of gradient iterations varies in each augmented Lagrangian step as shown
in :numref:`fig:DoubleInt`. The violation of the formulated
terminal equality constraints continuously decreases below the specified thresholds
:math:`\mb{\varepsilon_{\mb{g}_T}}=[1e-6,1e-6]^\mathsf{T}`.
As a result, the augmented cost functional and the original cost
functional converge to the same value, i.e. the so-called duality gap
is zero.

In a second scenario, :numref:`fig:DoubleInt_constr` shows the
optimal solution of OCP :math:numref:`eq:OCP_doubleInt` with
activated inequality constraint using the fixed end time :math:`T=5.25\,\mathrm{s}`. Further problem
settings are identical to the first simulation scenario. Again, the
terminal equality constraints are satisfied by the optimal solution and the state variables
:math:`\mb{x}(t)` are transferred to the origin. The control
variable :math:`u` is slightly adapted compared to the first simulation
scenario in :numref:`fig:DoubleInt` in order to comply with the inequality constraint.
In view of the additional inequality constraint, 17 augmented Lagrangian
steps are required to be able to solve the optimization problem with
sufficient accuracy.

The number of gradient iterations varies in each augmented Lagrangian
step as shown in :numref:`fig:DoubleInt_constr`. The
violation of the state constraints is almost continuously decreased
below the specified thresholds
:math:`\mb{\varepsilon_{\mb{g}_T}}=[1e-6,1e-6]^\mathsf{T}`
and :math:`\varepsilon_{\mb{h}}=1e-6`, respectively. As
before, the augmented cost functional and the original cost functional
converge to the same value. The computation time for solving the problem
on a Windows 10 machine with an Intel(R) Core(TM) i5-5300U CPU running
at 2.3GHz using the Microsoft Visual C++ 2013 Professional (C) compiler
amounts to 1.1ms and 14.6ms, respectively.

Optimization with free end time
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: ../img/tikz/DoubleInt_Timeopt_constr.*
    :name: fig:DoubleInt_Timeopt_constr

    Numerical OCP solution of the double integrator problem with free end time and state constraint.

In a third simulation scenario, the
OCP :math:numref:`eq:OCP_doubleInt` is numerically solved in
the free end time setting. The initial end time is set to
:math:`T=5.25\,\mathrm{s}` as given before. To weight the update of
the end time :math:`T` against the control update when updating the
optimization variables according to :math:numref:`eq:AlgOpt:LineSearchProblem`,
the adaptation factor :math:`\gamma_{T}` is set using the option ``OptimTimeLineSearchFactor``. 
The value :math:`\gamma_{T}=1.75` is used in the scenario. Note, however,
that values :math:`\gamma_{T}<1` typically increase the algorithmic
stability at the expense of the calculation time and vice versa.

The numerical results for the free end time case are shown in
:numref:`fig:DoubleInt_Timeopt_constr`. In contrast to the first
two simulation scenarios, the reduction of the end time below
:math:`4.6\,\mathrm{s}` allows one to carry out the setpoint
transition with a significantly more aggressive control
trajectory :math:`u`. The free end time optimization is a more
challenging problem than before and is accompanied by a higher
computational effort. This can be observed both in the larger number of
augmented Lagrangian steps and gradient steps as well as in terms of the
slower improvement of the violation of the state constraints.

Nevertheless, the state constraints are fulfilled at the last augmented
Lagrangian step in accordance with the thresholds
:math:`\mb{\varepsilon_{\mb{g}_T}}=[1e-6,1e-6]^\mathsf{T}`
and :math:`\varepsilon_{\mb{h}}=1e-6`. The improvement of
the control performance when optimizing the end time compared to a fixed
end time can be specified by the lower value of the cost functional,
cf. :numref:`fig:DoubleInt_constr` and :numref:`fig:DoubleInt_Timeopt_constr`. 
However, this results in a slightly increased computation time of 21.96ms compared to 14.66ms
with a fixed end time on the same Windows 10 machine.

.. footbibliography::