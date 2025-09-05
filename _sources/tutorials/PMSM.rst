.. _sec:TUT:PMSM:

Model predictive control of a PMSM
----------------------------------

The torque or current control of a permanent magnet synchronous machine
(PMSM) is a challenging example for nonlinear constrained model
predictive control. The following subsections illustrate the problem
formulation as well as useful options of GRAMPC to improve the control
performance. Corresponding m and C files can be found in the folder ``<grampc_root>/examples/PMSM``.

Problem formulation
~~~~~~~~~~~~~~~~~~~


.. figure:: ../img/tikz/PMSM_constr_unscaled.*
    :name: fig:PMSM_constr_unscaled

    Simulated MPC trajectories for the PMSM example with default settings.

The system dynamics of a PMSM :footcite:`Englert:CEP:2018`

.. math::
    :label: eq:tut:PMSM:model

    L_\text{d}\tfrac{{\rm d}}{{\rm d}t}i_\text{d} &= -R i_\text{d}+L_\text{q}\omega i_\text{q}+u_\text{d}\\ 
    L_\text{q}\tfrac{{\rm d}}{{\rm d}t}i_\text{q} &= -R i_\text{q}-L_\text{d}\omega i_\text{d}-\omega\psi_\text{p}+u_\text{q} \\
    J\tfrac{{\rm d}}{{\rm d}t}\omega &= \left(\tfrac{3}{2}z_\text{p}\left(\psi_\text{p}i_\text{q} + 
    i_\text{d}i_\text{q}(L_\text{d} - L_\text{q})\right)-\tfrac{\mu_\text{f}}{z_\text{p}}\omega-T_\text{L}\right)z_\text{p}\\
    \tfrac{{\rm d}}{{\rm d}t}\phi &= \omega

is given in the dq-coordinates. The system state
:math:`\mb{x} = [i_\text{d},i_\text{q},\omega,\phi]^\mathsf{T}`
comprises the dq-currents, the electrical rotor speed as well as the
electrical angle. The dq-voltages serve as controls
:math:`\mb{u} = [u_\text{d},u_\text{q}]^\mathsf{T}`. Further
system parameters are the stator resistance :math:`R = 3.5\,\Omega`,
the number of pole-pairs :math:`z_\text{p} = 3`, the permanent magnet
flux :math:`\psi_p = 0.17\,\mathrm{V\,s}`, the dq-inductivities
:math:`L_\text{d}=L_\text{q}=17.5\,\mathrm{mH}`, the moment of
inertia :math:`J= 0.9\,\mathrm{g\,m^2}` as well as the friction
coefficient :math:`\mu_\text{f}=0.4\,\mathrm{mN\,m\,s}`.

The magnitude of the dq-currents is limited by the maximum current
:math:`I_\text{max}=10\,A`, i.e.

.. math::
    :label: eq:tut:PMSM:stateconstr

    i_\text{abs} =i_\text{d}^2+i_\text{q}^2 \leq I_\text{max}^2\,,

in order prevent damage of the electrical components. Another constraint
concerns the dq-voltages. Through the modulation stage between the
controller and the voltage source inverter, the dq-voltages are limited
inside the circle

.. math::
    :label: eq:tut:PMSM:inputconstr

    u_\text{abs} = u_\text{d}^2+u_\text{q}^2 \leq U_\text{max}^2

with the maximum voltage :math:`U_\text{max}=323\,\mathrm{V}`.

The optimal control problem

.. math::

   	\min_{\mb{u}} \quad & J(\mb{u};\mb{x}_k) = \int_0^T l(\mb{x}(\tau), \mb{u}(\tau)) \, {\rm d}\tau

   	\textrm{s.t.} \quad & \mb{\dot x}(\tau) = \mb{f}(\mb{x}(\tau), \mb{u}(\tau)) \,, \quad \mb{x}(0) = \mb{x}_k

   	& h_1(\mb{x} (\tau)) = x_1(\tau)^2 +  x_2(\tau)^2 - I_\text{max}^2 \leq 0

   	& h_2(\mb{u} (\tau)) = u_1(\tau)^2 +  u_2(\tau)^2 - U_\text{max}^2 \leq 0

   	& \mb{u}(\tau) \in \left[\mb{u}_{\min}, \mb{u}_{\max}\right]

is subject to the system dynamics given by :math:numref:`eq:tut:PMSM:model` and the
constraints given by :math:numref:`eq:tut:PMSM:stateconstr` and
:math:numref:`eq:tut:PMSM:inputconstr`. 
The control constraints :math:numref:`eq:tut:PMSM:inputconstr` are
nonlinear and are therefore handled by the augmented Lagrangian
framework and not by the projection gradient method itself. It is
therefore reasonable to add the box constraints for :math:`\mb{u}` with :math:`\mb{u}_\text{min} = [-U_\text{max},-U_\text{max}]^\mathsf{T}` and
:math:`\mb{u} _\text{max} = [U_\text{max},U_\text{max}]^\mathsf{T}` to the OCP formulation to enhance the overall robustness of the
algorithm. The cost functional consists of the integral part

.. math::

   l(\mb{x}, \mb{u}) = q_1 (\mb{i}_\text{d} -  i_\text{d,des} )^2 + q_2 (\mb{i}_\text{q} -  i_\text{q,des} )^2 + (\mb{u} - \mb{u}_\text{des} ) ^\mathsf{T}\mb{R}  (\mb{u} - \mb{u}_\text{des} )\,,

with the setpoints for the states :math:`i_\text{d,des}`,
:math:`i_\text{q,des}` and controls
:math:`\mb{u}_\mathrm{des} \in \mathbb{R}^2` respectively. The
weights are set to :math:`q_1 = 8\,\mathrm{A^{-2}}`,
:math:`q_2=200\,\mathrm{A^{-2}}` and :math:`\mb{R} ={\rm diag}(0.001\,\mathrm{V^{-2}},0.001\,\mathrm{V^{-2}})`. The example considers a
startup of the motor from standstill by defining the setpoints
:math:`i_\text{d,des}=0\,\mathrm{A}` and
:math:`i_\text{q,des}=10\,\mathrm{A}`, corresponding to a constant
torque demand of :math:`7.65\,\mathrm{N\,m}`. The desired controls are set to
:math:`\mb{u}_\text{d,des}=[0\,\mathrm{V}, 0\,\mathrm{V}]^\mathsf{T}`.

The resulting OCP is solved by GRAMPC with the sampling time
:math:`\Delta t = 125\,\mathrm{\mu s}` (parameter ``dt``) and the horizon
:math:`T = 5\,\mathrm{ms}` (parameter ``Thor``) using standard options
almost exclusively. Only the number of discretization points ``Nhor=11`` , the
number of gradient iterations ``MaxGradIter=3``  and the number of augmented
Lagrangian iterations ``MaxMultIter=3``  are adapted to the problem. In addition, the
constraints tolerances ``ConstraintsAbsTol`` are set to 0.1% of the respective limit, i.e.
:math:`0.1\,\mathrm{A^2}` and :math:`104.5\,\mathrm{V^2}`. ``PenaltyMin`` is set to 2.5 x 10\ :sup:`-7` by
the estimation method of GRAMPC, see :ref:`sec:AlgOpt:EstimPenMin`.

:numref:`fig:PMSM_constr_unscaled` illustrates the simulation
results. The setpoints are reached very fast and are stabilized almost
exactly. However, an overshoot can be observed, which also leads to a
small violation of the dq-current constraint by :math:`0.45\,\mathrm{A}`. With increasing
rotor speed, the voltage also increases until the voltage constraint
becomes active. While the voltage constraint is almost exactly hold, the
dq-current constraint is clearly violated. The figure shows that at the
end of the simulation the dq-current constraint violation is more than
:math:`1\,\mathrm{A}` or 10%. Though a larger number of iterations might be used to reduce
the constraint violation, the main reason for this deviation is that the
nonlinear voltage and current constraints differ in several orders of
magnitude. The next section therefore shows how to scale the problem in
GRAMPC.

Also note that the increase of the cost functional does not indicate
instability, but can be explained by the increasing speed, which affects
the control term in the cost with the control setpoints
:math:`(\mb{u}_\text{d,des}=[0\,\mathrm{V}, 0\,\mathrm{V}]^\mathsf{T})`.
Moreover, the current setpoints
:math:`(i_\text{d,des}=0\,\mathrm{A},\,i_\text{q,des}=10\,\mathrm{A})`
cannot be hold due to the
constraints :math:numref:`eq:tut:PMSM:stateconstr` and
:math:numref:`eq:tut:PMSM:inputconstr`, which leads to
an additional cost increase.

Constraints scaling
~~~~~~~~~~~~~~~~~~~

.. figure:: ../img/tikz/PMSM_constr_scaled.*
    :name: fig:PMSM_constr_scaled

    Simulated MPC trajectories for the PMSM example with scaled constraints.

The two spherical constraints :math:numref:`eq:tut:PMSM:stateconstr` and
:math:numref:`eq:tut:PMSM:inputconstr` lie in very
different orders of magnitude, i.e. :math:`I_\text{max}^2 =100\,\mathrm{A^2}` 
and :math:`U_\text{max}^2 =104329\,\mathrm{V^2}`. 
Consequently, the two constraints should be scaled by the maximum value

.. math::

   \frac{i_\text{d}^2+i_\text{q}^2}{I_\text{max}^2}-1 \leq 0\,,\qquad  \frac{u_\text{d}^2+u_\text{q}^2}{U_\text{max}^2}-1 \leq 0\,.

This scaling can either be done by hand directly in the problem function
or by activating the option ``ScaleProblem`` and setting ``cScale`` :math:`=[I_\text{max}^2, U_\text{max}^2]`. 
The scaling option of GRAMPC, however, causes
additional computing effort (approx. 45% for the PMSM problem). Hence,
this option is suitable for testing the scaling, but eventually should
be done manually in the problem formulation to achieve the highest
computational efficiency. In accordance with the scaling of the
constraints, the tolerances are also adapted to 1 x 10\ :sup:`-3`
corresponding to 0.1% of the scaled constraints limits.

Besides the scaling and ``PenaltyMin`` that is set to 2 x 10\ :sup:`3` by the
estimation routine of GRAMPC, all parameters and options are the same
as in the last subsection. Please note that the estimation method for ``PenaltyMin``
strongly depends on reasonable constraint tolerances. In general, the
method returns rather conservative values, which may lead to constraint
violations if the order of magnitude of the constraints is very
different.

:numref:`fig:PMSM_constr_scaled` shows a clear improvement in
terms of the dq-current constraint that is now fully exploited. The only
violation results from the overshoot at the beginning, which is in the
same range as in the unscaled case (approx. 0.35 A). Further
improvements, e.g. reduction of the overshoot, can be achieved by
optimizing the penalty update as described in the next subsection.

Optimization of the penalty update
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: ../img/tikz/PMSM_penalty_opt.*
    :name: fig:PMSM_penalty_opt

    Simulated MPC trajectories for the PMSM example with scaled constraints and optimized penalty update.

In order to improve compliance with the dq-current constraint at the
beginning of the simulation, the number of augmented Lagrangian updates
is increased. To this end ``AugLagUpdateGradientRelTol`` is raised to 1, which means that in every
outer iteration an update of the multipliers and penalties is performed,
even if the inner minimization is not converged. Furthermore ``PenaltyMin`` , is raised
to 1 x 10\ :sup:`4` compared to the estimated value of 2 x 10\ :sup:`3`.
In addition, the plot of the step size, see the plot functions described
in :ref:`sec:Plotfunctions`, shows that the maximum value
:math:`\alpha_\text{max}` is often used. Consequently, setting the
maximum step size ``LineSearchMax`` to 10 allows larger optimization steps, especially at
the beginning and at the end of the simulation. All other parameters and
options, in particular the scaling options, are the same as in the
previous subsection.

:numref:`fig:PMSM_penalty_opt` illustrates the simulation result
with the optimized penalty update. The initial dq-current overshoot is
further reduced and the constraint is only violated by less than 0.07 A.
Furthermore, no oscillations occur in the costs and the augmented and
original cost are almost the same, which indicates that GRAMPC is well
tuned.

The computation time on a Windows 10 machine with Intel(R) Core(TM)
i5-5300U CPU running at 2.3 GHz using the Microsoft Visual C++ 2013
Professional (C) compiler amounts to 0.032 ms. On the dSpace real-time
hardware DS1202, the computation time is 0.13 ms.

.. footbibliography::