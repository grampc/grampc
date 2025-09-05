.. _sec:AlgOpt:ConvCheck:

Convergence criterion
---------------------

While the usage of fixed iteration counts :math:`i_\text{max}` and
:math:`j_\text{max}` for the outer and inner loops is typical in
real-time MPC or MHE applications, GRAMPC also provides an optional
convergence check that is useful for solving optimal control or
parameter optimization problems, or if the computation time in MPC is of
minor importance.

The inner gradient loop in :prf:ref:`alg:AlgOpt:GrampcAlgorithm`
evaluates the maximum relative gradient

.. math::
   :label: eq:AlgOpt:RelGrad

   \eta^{i|j+1} = 
   \max\left\{
   \frac{ \|\mb{u}^{i|j+1} - \mb{u}^{i|j}\|_{L_2} }
        { \|\mb{u}^{i|j+1}\|_{L_2} } \,,
   \frac{ \left\| \mb{p}^{i|j+1} - \mb{p}^{i|j} \right\|_2 }
        { \left\| \mb{p}^{i|j+1} \right\|_2 } \,,
   \frac{ |T^{i|j+1} - T^{i|j}| }
        { T^{i|j+1} } 
   \right\} 

in each iteration :math:`i|j` and terminates if

.. math::
    :label: eq:AlgOpt:ConvGradient

    \eta^{i|j} \leq \varepsilon_\text{rel,c} \,. 

Otherwise, the inner loop is continued until the maximum number of
iterations :math:`j_\text{max}` is reached. The last value
:math:`\eta^{i} = \eta^{i|j+1}` is returned to the outer loop and used
for the convergence check as well as for the update of multipliers and
penalties.

The augmented Lagrangian loop is terminated in iteration :math:`i` if
the inner loop is converged, that is
:math:`\eta^{i} \leq \varepsilon_\text{rel,c}`, and all constraints are
sufficiently satisfied, i.e.

.. math::
    :label: eq:AlgOpt:ConvConstraints

    \begin{bmatrix}
    |\mb{g}^{i}(t)| \\ \mb{\max}\{\mb{0}, \mb{\bar h}^{i}(t)\}
    \end{bmatrix} 
    \leq
    \begin{bmatrix}
    \mb{\varepsilon}_g \\ \mb{\varepsilon}_h
    \end{bmatrix}
    \, \forall t \in [0, T]
    \quad \land \quad
    \begin{bmatrix}
    |\mb{g}_T^{i}| \\ \mb{\max}\{\mb{0}, \mb{\bar h}_T^{i}\}
    \end{bmatrix} 
    \leq
    \begin{bmatrix}
    \mb{\varepsilon}_{g_T} \\ \mb{\varepsilon}_{h_T}
    \end{bmatrix} \,,

whereby the notation :math:`|\cdot|` denotes the
component-wise absolute value. Otherwise, the outer loop is continued
until the maximum number of iterations :math:`i_\text{max}` is reached.
The thresholds
:math:`\mb{\varepsilon}_g, \mb{\varepsilon}_h, \mb{\varepsilon}_{g_T}, \mb{\varepsilon}_{h_T}`
are vector-valued in order to rate each constraint individually.

The following options can be used to adjust the convergence criterion:

-  ``ConvergenceCheck``: This option activates the convergence criterion. Otherwise, the
   inner and outer loops always perform the maximum number of
   iterations, see the options ``MaxGradIter`` and ``MaxMultIter``.

-  ``ConvergenceGradientRelTol``: This option sets the threshold :math:`\varepsilon_\text{rel,c}` for
   the maximum relative gradient of the inner minimization problem that
   is used in the convergence criterion. Note that this threshold is
   different from the one that is used in the update of multipliers and
   penalties.

-  ``ConstraintsAbsTol``: Thresholds
   :math:`(\mb{\varepsilon_{\mb{g}}}, \mb{\varepsilon_{\mb{h}}}, \mb{\varepsilon_{\mb{g_T}}}, \mb{\varepsilon_{\mb{h_T}}}) \in \mathbb{R}^{N_{c}}`
   for the equality, inequality, terminal equality, and terminal
   inequality constraints.

.. _sec:AlgOpt:Scaling:

Scaling
-------

Scaling is recommended for improving the numerical conditioning when the
states :math:`\mb{x}` and the optimization variables
:math:`(\mb{u}, \mb{p}, T)` of the given optimization
problem differ in several orders of magnitude. Although GRAMPC allows
one to scale a specific problem automatically using the option
``ScaleProblem=1``, it should be noted that this typically increases the
computational load due to the additional multiplications in the
algorithm, cf. the tutorial on controlling a permanent magnet
synchronous machine in :ref:`sec:TUT:PMSM`. This issue can
be avoided by directly formulating the scaled problem within the C file
template ``probfct_TEMPLATE.c`` included in the folder
``examples/TEMPLATE``, also see :ref:`sec:ProblemImplementation`.

The scaling in GRAMPC is performed according to

.. math::
    :label: eq:AlgOpt:Scaling

    \bar{\mb{x}}(t) & = (\mb{x}(t) - \mb{x}_{\text{offset}}) \,./\, \mb{x}_{\text{scale}} \\
    \bar{\mb{u}}(t) & = (\mb{u}(t) - \mb{u}_{\text{offset}}) \,./\, \mb{u}_{\text{scale}} \\ 
    \bar{\mb{p}} & = (\mb{p} - \mb{p}_{\text{offset}}) \,./\, \mb{p}_{\text{scale}} \\
    \bar T & = \frac{T - T_{\text{offset}}}{T_{\text{scale}}} \,,

where :math:`\mb{x}_{\text{offset}} \in \mathbb{R}^x`,
:math:`\mb{u}_{\text{offset}} \in \mathbb{R}^u`,
:math:`\mb{p}_{\text{offset}} \in \mathbb{R}^p` and
:math:`\mb{T}_\text{offset} \in \mathbb{R}` denote offset values
and :math:`\mb{x}_{\text{scale}} \in \mathbb{R}^x`,
:math:`\mb{u}_{\text{scale}} \in \mathbb{R}^u`,
:math:`\mb{p}_{\text{scale}} \in \mathbb{R}^p` and
:math:`\mb{T}_\text{scale} \in \mathbb{R}` are scaling values.
The symbol :math:`./` in :math:numref:`eq:AlgOpt:Scaling`
denotes element-wise division by the scaling vectors.

Furthermore, GRAMPC provides a scaling factor :math:`J_\text{scale}`
for the cost functional as well as scaling factors
:math:`\mb{c}_\text{scale} = [\mb{c}_{\text{scale},\mb{g}}, \mb{c}_{\text{scale},\mb{h}}, \mb{c}_{\text{scale},\mb{g}_T}, \mb{c}_{\text{scale},\mb{h}_T}] \in \mathbb{R}^{N_c}`
for the constraints. The scaling of the cost functional is relevant as
the constraints are adjoined to the cost functional by means of
Lagrangian multipliers and penalty parameters and the original cost
functional should be of the same order of magnitude as these additional
terms.

The following options can be used to adjust the scaling:

-  ``ScaleProblem``: Activates or deactivates scaling. Note that GRAMPC requires more
   computation time if scaling is active.

-  ``xScale``, ``xOffset``: Scaling factors :math:`\mb{x}_\text{scale}` and offsets
   :math:`\mb{x}_\text{offset}` for each state variable.

-  ``uScale``, ``uOffset``: Scaling factors :math:`\mb{u}_\text{scale}` and offsets
   :math:`\mb{u}_\text{offset}` for each control variable.

-  ``pScale``, ``pOffset``: Scaling factors :math:`\mb{p}_\text{scale}` and offsets
   :math:`\mb{p}_\text{offset}` for each parameter.

-  ``TScale``, ``TOffset``: Scaling factor :math:`T_\text{scale}` and offset
   :math:`T_\text{offset}` for the horizon length.

-  ``JScale``: Scaling factor :math:`J_\text{scale}` for the cost functional.

-  ``cScale``: Scaling factors :math:`\mb{c}_\text{scale}` for each state
   constraint. The elements of the vector refer to the equality,
   inequality, terminal equality and terminal inequality constraints.

.. _sec:AlgOpt:ControlShift:

Control shift
-------------

The principle of optimality for an infinite horizon MPC problem
motivates to shift the control trajectory :math:`\mb{u}(t)`,
:math:`t\in[0,T]` from the previous MPC step :math:`k-1` by the sampling
time :math:`\Delta t` before the first GRAMPC iteration in the current
MPC step :math:`k`,
cf. :prf:ref:`alg:AlgOpt:GrampcAlgorithm`.
The last time segment of the shifted trajectory is hold on the last
value of the trajectory. Shifting the control can lead to a faster
convergence behavior of the gradient algorithm for many MPC problems.

If the control shift is activated for a problem with free end time
:math:`T`, GRAMPC assumes a shrinking horizon problem, because time
optimization is unusual in classical model predictive control. The
principle of optimality then motivates to subtract the sampling time
from the horizon :math:`T` after each MPC step, which corresponds to a
control shift for the end time.

-  ``ShiftControl``: Activates or deactivates the shifting of the control trajectory and
   the adaptation of :math:`T` in case of a free end time, i.e., if ``OptimTime`` is
   active.

.. _sec:AlgOpt:StatusFlags:

Status flags
------------

Several status flags are set in the solution structure
``grampc.sol.status`` during the execution of
:prf:ref:`alg:AlgOpt:GrampcAlgorithm`.
These flags can be printed as short messages by the function
``grampc_printstatus`` for the levels error, warn, info and debug.

The following status flags are printed on the level
``STATUS_LEVEL_ERROR`` and require immediate action:

-  ``STATUS_INTEGRATOR_INPUT_NOT_CONSISTENT``: This flag is set by the integrator ``rodas`` if the input values
   are not consistent. See :footcite:`Rodas:Webpage:2018` for
   further details.

-  ``STATUS_INTEGRATOR_MAXSTEPS``: This flag is set by the integrators ``ruku45`` or ``rodas`` if too
   many steps are required.

-  ``STATUS_INTEGRATOR_STEPS_TOO_SMALL``: This flag is set by the integrator ``rodas`` if the step size
   becomes too small.

-  ``STATUS_INTEGRATOR_MATRIX_IS_SINGULAR``: This flag is set by the integrator ``rodas`` if a singular Jacobian
   :math:`\frac{\partial \mb{f}}{\partial \mb{x}}` or
   :math:`\left(\frac{\partial \mb{f}}{\partial \mb{x}}\right)^\mathsf{T}`
   is detected. See :footcite:`Rodas:Webpage:2018` for further
   details.

-  ``STATUS_INTEGRATOR_H_MIN``: This flag is set by the integrator ``ruku45`` if a smaller step
   size than the minimal allowed value is required.

The following flags are printed in addition to the previous ones on the
level ``STATUS_LEVEL_WARN``:

-  ``STATUS_MULTIPLIER_MAX``: This flag is set if one of the multipliers
   :math:`\mb{\bar \mu}` reaches the upper limit
   :math:`\mu_\text{max}` or the lower limit :math:`-\mu_\text{max}`.
   The situation may occur for example if the problem is infeasible or
   ill-conditioned or if the penalty parameters are too high.

-  ``STATUS_PENALTY_MAX``: This flag is set if one of the penalty parameters
   :math:`\mb{\bar c}` reaches the upper limit
   :math:`c_\text{max}`. The situation may occur for example if the
   problem is infeasible or ill-conditioned or if the penalty increase
   factor :math:`\beta_{\mathrm{in}}` is too high.

-  ``STATUS_INFEASIBLE``: This flag is set if the constraints are not satisfied and one run
   of :prf:ref:`alg:AlgOpt:GrampcAlgorithm`
   does not reduce the norm of the constraints. The situation may occur
   in single runs, if few iterations :math:`i_\text{max}` and
   :math:`j_\text{max}` are used for a suboptimal solution. However, if
   the flag is set in multiple successive runs, it is a strong indicator
   for an infeasible optimization problem.

The following flags are printed in addition to the previous ones on the
level ``STATUS_LEVEL_INFO``:

-  ``STATUS_GRADIENT_CONVERGED``: This flag is set if the convergence check is activated and the
   relative tolerance :math:`\varepsilon_\text{rel,c}` is satisfied for
   the controls :math:`\mb{u}`, the parameters
   :math:`\mb{p}`, and the end time :math:`T`.

-  ``STATUS_CONSTRAINTS_CONVERGED``: This flag is set if the convergence check is activated and the
   absolute tolerances
   :math:`\mb{\varepsilon}_g, \mb{\varepsilon}_h, \mb{\varepsilon}_{g_T}, \mb{\varepsilon}_{h_T}`
   are satisfied for all constraints.

-  ``STATUS_LINESEARCH_INIT``: This flag is set if the gradient algorithm uses the initial step
   size :math:`\alpha_\text{init}` as fallback for the explicit line
   search strategy in one iteration.

The following flags are printed in addition to the previous ones on the
level ``STATUS_LEVEL_DEBUG``:

-  ``STATUS_LINESEARCH_MAX``: This flag is set if the gradient algorithm uses the maximum step
   size :math:`\alpha_\text{max}` in one iteration.

-  ``STATUS_LINESEARCH_MIN``: This flag is set if the gradient algorithm uses the minimum step
   size :math:`\alpha_\text{min}` in one iteration.

-  ``STATUS_MULTIPLIER_UPDATE``: This flag is set if the relative tolerance
   :math:`\varepsilon_\text{rel,u}` is satisfied for the controls
   :math:`\mb{u}`, the parameters :math:`\mb{p}` and the
   end time :math:`T` and therefore the update of the multipliers
   :math:`\mb{\bar \mu}` and the penalty parameters
   :math:`\mb{\bar c}` is performed, cf. :ref:`sec:AlgOpt:UpdateMultPen`.

.. footbibliography::