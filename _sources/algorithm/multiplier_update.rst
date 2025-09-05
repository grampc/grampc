.. _sec:AlgOpt:UpdateMultPen:

Update of multipliers and penalties
-----------------------------------

GRAMPC handles general nonlinear constraints using an augmented
Lagrangian approach or, alternatively, using an external penalty method.
The key to the efficient solution of constrained problems using these
approaches are the updates of the multipliers in the outer loop for
:math:`i = 1,\, \dots,\, i_\text{max}`.

.. _sec:AlgOpt:UpdateMult:

Update of Lagrangian multipliers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The outer loop of the GRAMPC algorithm in :prf:ref:`alg:AlgOpt:GrampcAlgorithm`
maximizes the augmented Lagrangian function with respect to the
multipliers :math:`\mb{\bar \mu}`. This update is carried out in
the direction of steepest ascent that is given by the constraint
residual. The penalty parameter is used as step size, as it is typically
done in augmented Lagrangian methods.

For an arbitrary equality constraint
:math:`g^{i} = g(\mb{x}^{i}, \dots)` with multiplier
:math:`\mu_g^{i}`, penalty :math:`c_g^i`, and tolerance
:math:`\varepsilon_g`, the update is defined by

.. math::

   \mu_g^{i+1} = \zeta_{g}(\mu_g^{i}, c_g^{i}, g^{i}, \varepsilon_g) 
   = 
   \begin{cases}
   \mu_g^{i} + (1-\rho) c_g^{i} g^{i}
   & \text{if } \left|g^{i}\right| >  \varepsilon_g
   \, \land \,
   \eta^{i} \leq \varepsilon_\text{rel,u}
   \\
   {\mu}_g^{i}
   & \text{else} \,.
   \end{cases}

The update is not performed if the constraint is satisfied within its
tolerance :math:`\varepsilon` or if the inner minimization is not
sufficiently converged, which is checked by the maximum relative
gradient :math:`\eta^i` (see :math:numref:`eq:AlgOpt:RelGrad` for the definition) and the
threshold :math:`\varepsilon_\text{rel,u}`. Similarly, for an inequality
constraint :math:`\bar h^{i} = \bar h(\mb{x}^{i}, \dots)` with multiplier
:math:`\mu_h^{i}`, penalty :math:`c_h^{i}`, and tolerance
:math:`\varepsilon_h`, the update is defined by

.. math::

   \mu_h^{i+1} = \zeta_{h}(\mu_h^{i}, c_h^{i}, \bar h^{i}, \varepsilon_h) 
   = 
   \begin{cases}
   \mu_h^{i} + (1-\rho) c_h^{i} \bar h^{i}
   & \text{if } \left(\bar h^{i} > \varepsilon_h
   \, \land \, 
   \eta^{i} \leq \varepsilon_\text{rel,u} \right) 
   \, \lor \, 
   \bar h^{i} < 0
   \\
   \mu_h^{i}
   & \text{else} \,.
   \end{cases}

Similar update rules are used for the terminal equality and terminal
inequality constraints.

GRAMPC provides several means to increase the robustness of the
multiplier update, which may be required if few iterations
:math:`j_\text{max}` are used for the suboptimal solution of the inner
minimization problem. The damping factor :math:`\rho \in [0, 1)` can be
used to scale the step size of the steepest ascent and the tolerance
:math:`\varepsilon_\text{rel,u}` can be used to skip the multiplier
update in case that the minimization is not sufficiently converged.
Furthermore, the multipliers are limited by lower and upper bounds
:math:`\mu_g \in [-\mu_\text{max}, \mu_\text{max}]` for equalities and
:math:`\mu_h \leq \mu_\text{max}` for inequalities, respectively, to
avoid unlimited growth. A status flag is set if one of the multipliers
reaches this bound and the user should check the problem formulation as
this case indicates an ill-posed or even infeasible optimization
problem.

The following options can be used to adjust the update of the Lagrangian
multipliers:

-  ``MultiplierMax``: Upper bound :math:`\mu_\text{max}` and lower bound
   :math:`-\mu_\text{max}` for the Lagrangian multpliers.

-  ``MultiplierDampingFactor``: Damping factor :math:`\rho \in [0,1)` for the multiplier update.

-  ``AugLagUpdateGradientRelTol``: Threshold :math:`\varepsilon_\text{rel,u}` for the maximum relative
   gradient of the inner minimization problem.

-  ``ConstraintsAbsTol``: Thresholds
   :math:`(\mb{\varepsilon_{\mb{g}}}, \mb{\varepsilon_{\mb{h}}}, \mb{\varepsilon_{\mb{g_T}}}, \mb{\varepsilon_{\mb{h_T}}}) \in \mathbb{R}^{N_{c}}`
   for the equality, inequality, terminal equality, and terminal
   inequality constraints.

.. _sec:AlgOpt:UpdatePen:

Update of penalty parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The penalty parameters :math:`\mb{\bar c}` are adapted in each
outer iteration according to a heuristic rule that is motivated by the
LANCELOT package :footcite:`Conn2013,Nocedal2006`. A carefully
tuned adaptation of the penalties can speed-up the convergence
significantly and is therefore highly recommended (also see
:ref:`sec:AlgOpt:EstimPenMin` and the tutorial in
:ref:`sec:TUT:PMSM`). Note that the penalty parameters are
also updated if external penalties are used instead of the augmented
Lagrangian method, i.e. ``ConstraintsHandling`` is set to ``extpen``. In order to keep the
penalty parameters at the initial value ``PenaltyMin``, the options ``PenaltyIncreaseFactor`` and ``PenaltyDecreaseFactor`` can be set to
:math:`1.0`, which basically deactivates the penalty update.

For an arbitrary equality constraint
:math:`g^{i} = g(\mb{x}^{i}, \dots)` with penalty :math:`c_g^i`
and tolerance :math:`\varepsilon_g`, the update is defined by

.. math::
    :label: eq:AlgOpt:UpdatePeng

    c_g^{i+1} = \xi_g(c_g^{i}, g^{i}, g^{i-1}, \varepsilon_g) 
    =
    \begin{cases}
    \beta_{\mathrm{in}} \, c_g^{i}
    &\text{if } \left|g^{i}\right| \geq \gamma_{\mathrm{in}} \left|g^{i-1}\right| 
    \, \land \,
    \left|g^{i}\right| > \varepsilon_g
    \, \land \,
    \eta^{i} \leq \varepsilon_\text{rel,u}
    \\
    \beta_{\mathrm{de}} \, c_g^{i}
    &\text{else if } \left|g^{i}\right| \leq \gamma_{\mathrm{de}} \, \varepsilon_g
    \\
    c_g^{i}
    & \textrm{else} \,.
    \end{cases}

The penalty :math:`c_g^{i}` is increased by the factor
:math:`\beta_\text{in}> 1` if the (sub-optimal) solution of the inner minimization problem
does not generate sufficient progress in the constraint, which is rated
by the factor :math:`\gamma_\text{in} > 0` and compared to the previous
iteration :math:`i-1`. This update is skipped if the inner minimization
is not sufficiently converged, which is checked by the maximum relative
gradient :math:`\eta^i` and the threshold
:math:`\varepsilon_\text{rel,u}`. The penalty :math:`c_g^{i}` is
decreased by the factor :math:`\beta_\text{de} < 1` if the constraint
:math:`g^{i}` is sufficiently satisfied within its tolerance, whereby
currently the constant factor :math:`\gamma_\text{de} = 0.1` is used.
The setting :math:`\beta_\text{in} = \beta_\text{de} = 1` can be used to
keep the penalty constant, i.e., to deactivate the penalty adaptation.
Similarly, for an inequality constraint
:math:`\bar h^{i} = \bar h(\mb{x}^{i}, \dots)` with penalty
:math:`c_h^i` and tolerance :math:`\varepsilon_h`, the update is defined
by

.. math::
    :label: eq:AlgOpt:UpdatePenh

    c_h^{i+1} = \xi_h(c_h^{i}, \bar h^{i}, \bar h^{i-1}, \varepsilon_h) 
    = 
    \begin{cases}
    \beta_{\mathrm{in}} \, c_h^{i}
    &\text{if } \bar h^{i} \geq \gamma_{\mathrm{in}} \bar h^{i-1}
    \, \land \,
    \bar h^{i} > \varepsilon_h
    \, \land \,
    \eta^{i} \leq \varepsilon_\text{rel,u}
    \\
    \beta_{\mathrm{de}} \, c_h^{i}
    & \text{else if } \bar h^{i} \leq \gamma_{\mathrm{de}} \, \varepsilon_h
    \\
    c_h^{i}
    & \textrm{else} \,.
    \end{cases}

Similar update rules are used for the terminal equality and inequality
constraints. In analogy to the multiplier update, the penalty parameters
are restricted to upper and lower bounds :math:`c_\text{max} \gg
c_\text{min} > 0` in order to avoid unlimited growth as well as
negligible values.

The following options can be used to adjust the update of the penalty
parameters:

-  ``PenaltyMax``: This option sets the upper bound :math:`c_\text{max}` of the
   penalty parameters.

-  ``PenaltyMin``: This option sets the lower bound :math:`c_\text{min}` of the
   penalty parameters.

-  ``PenaltyIncreaseFactor``: This option sets the factor :math:`\beta_\text{in}` by which
   penalties are increased.

-  ``PenaltyDecreaseFactor``: This option sets the factor :math:`\beta_\text{de}` by which
   penalties are decreased.

-  ``PenaltyIncreaseThreshold``: This option sets the factor :math:`\gamma_\text{in}` that rates the
   progress in the constraints between the last two iterates.

-  ``AugLagUpdateGradientRelTol``: Threshold :math:`\varepsilon_\text{rel,u}` for the maximum relative
   gradient of the inner minimization problem.

-  ``ConstraintsAbsTol``: Thresholds
   :math:`(\mb{\varepsilon_{\mb{g}}}, \mb{\varepsilon_{\mb{h}}}, \mb{\varepsilon_{\mb{g_T}}}, \mb{\varepsilon_{\mb{h_T}}}) \in \mathbb{R}^{N_{c}}`
   for the equality, inequality, terminal equality, and terminal
   inequality constraints.

.. _sec:AlgOpt:EstimPenMin:

Estimation of minimal penalty parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In real-time or embedded MPC applications, where only a limited number
of iterations per step is computed, it is crucial that the penalty
parameter is not decreased below a certain threshold
:math:`c_\text{min}`. This lower bound should be large enough that an
inactive constraint that becomes active is still sufficiently penalized
in the augmented Lagrangian cost functional. However, it should not be
chosen too high to prevent ill-conditioning. A suitable value of
:math:`c_\text{min}` tailored to the given MPC problem therefore is of
importance to ensure a high performance of GRAMPC.

In order to support the user, GRAMPC offers the routine ``grampc_estim_penmin`` to compute a
problem-specific estimate of :math:`c_\text{min}`. The basic idea behind
this estimation is to determine :math:`c_\text{min}` such that the
actual costs :math:`J` are in the same order of magnitude as the squared
constraints multiplied by :math:`c_\text{min}`, see
equation :math:numref:`eq:AlgOpt:AugLag`. This approach
requires initial values for the states :math:`\mb{x}`, controls
:math:`\mb{u}`, and cost :math:`J`. If the GRAMPC structure
includes only default values, i.e. zeros (cold start), the estimation
function ``grampc_estim_penmin`` can be called with the argument ``rungrampc=1`` to perform one optimization
or MPC step, where the possible maximum numbers of gradient iterations ``MaxGradIter``
and augmented Lagrangian iterations ``MaxMultIter`` are limited to 20. Afterwards, the
estimated value for ``PenaltyMin`` is set as detailed below and, if ``rungrampc=1``, the initial
states :math:`\mb{x}`, controls :math:`\mb{u}` and costs
:math:`J` are reset.

Based on the initial values, a first estimate of the minimal penalty
parameter is computed according to

.. math::
    :label: eq:AlgOpt:C1

    \hat{c}_\text{min}^\text{I} 
    %= \frac{\tfrac{1}{2}\,|J|\,(N_{\vm{g}}+N_{\vm{h}})}
    %{\|\vm g(\vm x(t), \vm u(t), \vm p, t)\|_{L_1}^2+\|\vm h(\vm x(t), \vm u(t), \vm p, t)\|_{L_1}^2 }
    %+ \frac{\tfrac{1}{2}\,|J|\,(N_{\vm{g}_T}+N_{\vm{h}_T})} {\|\vm g_T(\vm x(T), \vm p, T)\|_{1}^2+\|\vm h_T(\vm x(T), \vm p, T)\|_{1}^2}
    = \frac{2\,|J| %\,(N_{\vm{g}} + N_{\vm{h}} + N_{\vm{g}_T} + N_{\vm{h}_T})
    }
    {\|\mb{g}(\mb{x}(t), \mb{u}(t), \mb{p}, t)\|_{L_2}^2 + \|\mb{h}(\mb{x}(t), \mb{u}(t), \mb{p}, t)\|_{L_2}^2 + \|\mb{g}_T(\mb{x}(T), \mb{p}, T)\|_{2}^2 + \|\mb{h}_T(\mb{x}(T), \mb{p}, T)\|_{2}^2}

However, if the inequality constraints are initially inactive and are
far away from their bounds (i.e. large negative values are returned),
the estimate :math:`\hat{c}_\text{min}^\text{I}` may be too small.

To deal with these cases, a second estimate for the minimal penalty
parameter

.. math::
    :label: eq:AlgOpt:C2

    \hat{c}_\text{min}^\text{II} 
    %= \kappa \left(\frac{\tfrac{1}{2}\,|J|}
    %{T \max\{\|\vm{\varepsilon}_{\vm{g}}\|_1,\,\|\vm{\varepsilon}_{\vm{h}}\|_1\}^2}
    %+ \frac{\tfrac{1}{2}\,|J|} 
    %{\max\{\|\vm{\varepsilon}_{\vm{g}_T}\|_1,\,\|\vm{\varepsilon}_{\vm{h}_T}\|_1\}^2}\right)
    = \frac{2\,|J|}
    {
    T \left(\|\mb{\varepsilon}_{\mb{g}}\|_2^2 +
    \|\mb{\varepsilon}_{\mb{h}}\|_2^2\right) +
    \|\mb{\varepsilon}_{\mb{g}_T}\|_2^2 +
    \|\mb{\varepsilon}_{\mb{h}_T}\|_2^2}

is computed in the same spirit using the constraint tolerances (see ``ConstraintsAbsTol``,
:math:`\mb{\varepsilon}_{\mb{g}}`,
:math:`\mb{\varepsilon}_{\mb{h}}`,
:math:`\mb{\varepsilon}_{\mb{g}_T}` and
:math:`\mb{\varepsilon}_{\mb{h}_T}`) instead of the
constraint values [2]_. Since the norms of the tolerances are summed,
more conservative values for :math:`c_\text{min}` are estimated and
therefore instabilities can be avoided. Note that it is recommended to
scale all constraints so that they are in the same order of magnitude,
see e.g. the PMSM example in :ref:`sec:TUT:PMSM`.

Finally, the minimal penalty parameter

.. math::

   \hat{c}_\text{min} = \min\left\{
   \max\left\{\hat{c}_\text{min}^\text{I} ,\,
   \kappa \, \hat{c}_\text{min}^\text{II} \right\} ,\, 
   \frac{c_\text{max}}{500}\right\}

is chosen as the maximum of :math:numref:`eq:AlgOpt:C1` and
:math:numref:`eq:AlgOpt:C2` and additionally limited to
:math:`0.2\%` of the maximum penalty parameter :math:`c_\text{max}`.
This limitation ensures reasonable values even with very small
constraint tolerances. The relation factor :math:`\kappa` has been
determined to :math:`10^{-6}` on the basis of various example systems.
Please note that this estimation is intended to assist the user in
making an initial guess. Problem-specific tuning of ``PenaltyMin`` can lead to further
performance improvements and is therefore recommended. All MPC example
problems in ``<grampc_root>/examples`` contain an initial call of ``grampc_estim_penmin`` to estimate
:math:`\hat{c}_\text{min}` and, as alternative, manually tuned values
that can further enhance the performance of GRAMPC for fixed numbers
of iterations.

.. rubric:: Footnotes

.. [2] The integration behind the :math:`L^2`-norm can be replaced by a multiplication by the horizon length :math:`T`, as the constraint tolerances :math:`\epsilon_g` and :math:`\epsilon_h` are no functions of time.

.. footbibliography::