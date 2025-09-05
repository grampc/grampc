.. _sec:AlgOpt:BasicAlgorithm:

Optimization algorithm
----------------------

The optimization algorithm of GRAMPC is based on an augmented
Lagrangian formulation with an inner projected gradient method as
minimization step and an outer multiplier and penalty update. This
section gives a brief sketch of the algorithm. Note that a more detailed
description is given in :footcite:`Englert:OE:2019`.

.. _sec:AlgOpt:AugLag:

Augmented Lagrangian method
~~~~~~~~~~~~~~~~~~~~~~~~~~~

GRAMPC implements the augmented Lagrangian approach to handle the equality and inequality constraints of the OCP. 
The constraints are adjoined to the integral cost function using the time-dependent multipliers
:math:`\mb{\mu}= [\mb{\mu}_{\mb{g}}^\mathsf{T}, \mb{\mu}_{\mb{h}}^\mathsf{T}]^\mathsf{T}`
and penalties :math:`\mb{c}= [\mb{c}_{\mb{g}}^\mathsf{T},\mb{c}_{\mb{h}}^\mathsf{T}]^\mathsf{T}`. 
Similarly, multipliers
:math:`\mb{\mu}_T = [\mb{\mu}_{\mb{g}_T}^\mathsf{T}, \mb{\mu}_{\mb{h}_T}^\mathsf{T}]^\mathsf{T}`
and penalties
:math:`\mb{c}_T = [\mb{c}_{\mb{g}_T}^\mathsf{T}, \mb{c}_{\mb{h}_T}^\mathsf{T}]^\mathsf{T}`
are used for the terminal constraints. 
Where appropriate, the syntax
:math:`\mb{\bar \mu} = (\mb{\mu}_{\mb{g}}, \mb{\mu}_{\mb{h}}, \mb{\mu}_{\mb{g}_T}, \mb{\mu}_{\mb{h}_T})`
and
:math:`\mb{\bar c} = (\mb{c}_{\mb{g}}, \mb{c}_{\mb{h}}, \mb{c}_{\mb{g}_T}, \mb{c}_{\mb{h}_T})`
is used to denote all multipliers and penalties. 
The algorithm requires a reformulation of the inequality constraints that leads to the transformed functions (see :footcite:`Englert:OE:2019` for details)

.. math::

   \mb{\bar h}(\mb{x}, \mb{u}, \mb{p}, t, \mb{\mu}_{\mb{h}}, \mb{c}_{\mb{h}}) &=
   \mb{\max}\left\{ \mb{h}(\mb{x}, \mb{u}, \mb{p}, t), -\mb{C}_{\mb{h}}^{-1} \mb{\mu}_{\mb{h}} \right\}
   \\
   \mb{\bar h}_T(\mb{x}, \mb{p}, T, \mb{\mu}_{\mb{h}_T}, \mb{c}_{\mb{h}_T}) &= 
   \mb{\max}\left\{ \mb{h}_T(\mb{x}, \mb{p}, T), -\mb{C}_{\mb{h}_T}^{-1} \mb{\mu}_{\mb{h}_T} \right\}

with the component-wise **max**-function and the diagonal matrix syntax :math:`\mb{C}= {\rm diag}(\mb{c})`. 
The augmented Lagrangian function is defined as

.. math::
    :label: eq:AlgOpt:AugLag

    \bar J(\mb{u}, \mb{p}, T, \mb{\bar \mu}, \mb{\bar c};\mb{x}_0) = 
    \bar V(\mb{x}, \mb{p}, T, \mb{\mu}_T, \mb{c}_T) 
    + \int_0^T \bar l(\mb{x}, \mb{u}, \mb{p}, t, \mb{\mu}, \mb{c}) \, \mathrm dt

with the augmented terminal cost term

.. math::

   \bar V(\mb{x}, \mb{p}, T, \mb{\mu}_T, \mb{c}_T) =\,& V(\mb{x}, \mb{p}, T) 
   + \mb{\mu}_{\mb{g}_T}^\mathsf{T}\, \mb{g}_T(\mb{x}, \mb{p}, T)
   + \frac12 \| \mb{g}_T(\mb{x}, \mb{p}, T) \|^2_{\mb{C}_{\mb{g}_T}} 
   \nonumber\\
   &+ \mb{\mu}_{\mb{h}_T}^\mathsf{T}\, \mb{\bar h}_T(\mb{x}, \mb{p}, T, \mb{\mu}_{\mb{h}_T},\mb{c}_{\mb{h}_T})
   + \frac12 \| \mb{\bar h}_T(\mb{x}, \mb{p}, T, \mb{\mu}_{\mb{h}_T},\mb{c}_{\mb{h}_T}) \|^2_{\mb{C}_{\mb{h}_T}}

and the augmented integral cost term

.. math::

   \bar l(\mb{x}, \mb{u}, \mb{p}, t, \mb{\mu}, \mb{c}) =\,& l(\mb{x}, \mb{u}, \mb{p}, t)
   + \mb{\mu}_{\mb{g}}^\mathsf{T}\, \mb{g}(\mb{x}, \mb{u}, \mb{p}, t)
   + \frac12 \| \mb{g}(\mb{x}, \mb{u}, \mb{p}, t) \|^2_{\mb{C}_{\mb{g}}}
   \nonumber\\
   &+ \mb{\mu}_{\mb{h}}^\mathsf{T}\, \mb{\bar h}(\mb{x}, \mb{u}, \mb{p}, t, \mb{\mu_h},\mb{c_h})
   + \frac12 \| \mb{\bar h}(\mb{x}, \mb{u}, \mb{p}, t, \mb{\mu_h},\mb{c_h}) \|^2_{\mb{C}_{\mb{h}}} \,.

Instead of solving the original problem, the algorithm solves the max-min-problem

.. math::
    :label: eq:MaxMin

    \max_{\mb{\bar\mu}} \, \min_{\mb{u}, \mb{p}, T} \quad& 
    \bar{ J}(\mb{u}, \mb{p}, T, \mb{\bar\mu}, \mb{\bar c}; \mb{x}_0) 
    \\ 
    \textrm{s.t.} \quad& \mb{M} \mb{\dot x}(t) = \mb{f}(\mb{x}, \mb{u}, \mb{p}, t) 
    \,,\quad 
    \mb{x}(0) = \mb{x}_0
    \\
    & \mb{u}(t) \in [\mb{u}_{\min}, \mb{u}_{\max}] %\,,\quad t \in [0, T]
    \,,\quad 
    t\in[0,T]
    \\
    & \mb{p} \in [\mb{p}_{\min}, \mb{p}_{\max}]
    \,,\quad 
    T \in [T_{\min}, T_{\max}] \,,

whereby the augmented Lagrangian function :math:`\bar{J}` is maximized with respect to the multipliers :math:`\mb{\bar \mu}` and minimized with respect to the controls :math:`\mb{u}`, the parameters :math:`\mb{p}` and the end time :math:`\mb{T}`. 
Note that the full set of optimization variables :math:`(\mb{u},\mb{p},T)` is considered in what follows for the sake of completeness. 
The max-min-problem :math:numref:`eq:MaxMin` corresponds to the dual problem of :math:numref:`eq:OCP` in the case of :math:`\mb{\bar c} = \mb{0}`. 
The maximization step is performed by steepest ascent using the constraint residual as direction and the penalty parameter as step size. 
See :footcite:`Englert:OE:2019` for a detailed description of the augmented Lagrangian algorithm.

.. _sec:AlgOpt:ProjGrad:

Projected gradient method
~~~~~~~~~~~~~~~~~~~~~~~~~

GRAMPC uses a projected gradient method to solve the inner
minimization problem subject to the dynamics
:math:`\mb{f}(\mb{x}(t), \mb{u}(t), \mb{p}, t)` as well as the box constraints
:math:`\mb{u}(t) \in \left[\mb{u}_{\min}, \mb{u}_{\max}\right]` and
:math:`\mb{p} \in \left[\mb{p}_{\min}, \mb{p}_{\max}\right]`. 
The algorithm is based on the first-order optimality conditions that can be compactly stated using the Hamiltonian

.. math::

   H(\mb{x}, \mb{u}, \mb{p}, \mb{\lambda}, t, \mb{\mu}, \mb{c}) = 
   \bar{ l}(\mb{x}, \mb{u}, \mb{p}, t, \mb{\mu}, \mb{c}) + \mb{\lambda}^\mathsf{T}\mb{f}(\mb{x}, \mb{u}, \mb{p}, t)


with the adjoint states :math:`\mb{\lambda}`. 
The canonical equations are then given by

.. math::
    :label: eq:AlgOpt:OptCondLambda

    \mb{M}\mb{\dot x} &= \mb{f} (\mb{x}, \mb{u}, \mb{p}, t) \,, &\mb{x}(0) &= \mb{x}_0 \,,

    \mb{M}^\mathsf{T}\mb{\dot \lambda} &= -H_{\mb{x}}(\mb{x}, \mb{u}, \mb{p}, \mb{\lambda}, t, \mb{\mu}, \mb{c}) \,, &\mb{M}^\mathsf{T}\mb{\lambda}(T) &= \bar{ V}_{\mb{x}}(\mb{x}(T), \mb{p}, T, \mb{\mu}_T, \mb{c}_T)

consisting of the original dynamics :math:`\mb{\dot x}` and the adjoint dynamics :math:`\mb{\dot \lambda}`. 
The canonical equations can be iteratively solved in forward and backward time for given initial values of the optimization variables. 
In each iteration and depending on the optimization variables of the actual problem to be solved, the gradients

.. math::

    \mb{d}_{\mb{u}} &= H_{\mb{u}}(\mb{x}, \mb{u}, \mb{p}, \mb{\lambda}, t, \mb{\mu}, \mb{c}) 
    \\
    \mb{d}_{\mb{p}} &= \bar V_{\mb{p}}(\mb{x}(T), \mb{p}, T, \mb{\mu}_T, \mb{c}_T) +     \int_0^T H_{\mb{p}}(\mb{x}, \mb{u}, \mb{p}, \mb{\lambda}, t, \mb{\mu}, \mb{c}) \, {\rm d}t 
    \\
    d_T &= \bar V_T(\mb{x}(T), \mb{p}, T, \mb{\mu}_T, \mb{c}_T) + H(\mb{x}(T), \mb{u}(T), \mb{p}, \mb{\lambda}(T), T, \mb{\mu}(T), \mb{c}(T))

with respect to the controls :math:`\mb{u}`, parameters :math:`\mb{p}`, 
and end time :math:`T` are used to formulate a line search problem

.. math::

   \min_{\alpha} \bar{J} \left( 
     \mb{\psi}_{\mb{u}} \left( \mb{u} - \alpha \mb{d}_{\mb{u}} \right),
     \mb{\psi}_{\mb{p}} \left( \mb{p} - \gamma_{\mb{p}} \alpha \mb{d}_{\mb{p}} \right),
     \psi_{T} \left( T - \gamma_{T} \alpha d_{T} \right); 
     \mb{\bar \mu}, \mb{\bar c}, \mb{x}_0
   \right)

with projection functions :math:`\mb{\psi}_{\mb{u}}`,
:math:`\mb{\psi}_{\mb{p}}` and :math:`\psi_{T}` and,
finally, to update the optimization variables according to

.. math::

   \mb{u} \leftarrow \mb{\psi}_{\mb{u}} \left( \mb{u} - \alpha \mb{d}_{\mb{u}} \right)
   \,,\quad
   \mb{p} \leftarrow \mb{\psi}_{\mb{p}} \left( \mb{p} - \gamma_{\mb{p}} \alpha \mb{d}_{\mb{p}} \right)
   \,,\quad
   T \leftarrow \psi_{T} \left( T - \gamma_{T} \alpha d_{T} \right) \,.

See :footcite:`InTech_GraichenKaepernick2012,Kaepernick2014,Englert:OE:2019`
for a detailed description of the projected gradient algorithm. 
GRAMPC provides two methods for the approximate solution of the line search
problem, which are explained in :ref:`sec:AlgOpt:LineSearch`.

.. _sec:AlgOpt:Structure:

Algorithmic structure
~~~~~~~~~~~~~~~~~~~~~

.. prf:algorithm:: Basic algorithmic structure of GRAMPC.
    :label: alg:AlgOpt:GrampcAlgorithm

    Optional: Shift trajectories by sampling time :math:`\Delta t`

    **For** :math:`i = 1` to :math:`i_{max}` **do**
        **For** :math:`j = 1` to :math:`j_{max}` **do**
            **If** :math:`i > 1` and :math:`j = 1` **then**
                Set :math:`\mb x^{i|j} =\mb x^{i-1}`
            **else**
                Compute :math:`\mb x^{i|j}` by forward time integration of system dynamics

                Evaluate all constraints

            **End If**
            
            Compute :math:`\mb \lambda^{i|j}` by backward time integration of adjoint system

            Evaluate gradients :math:`\mb{d_u}^{i|j}`, :math:`\mb{d_p}^{i|j}` and :math:`d_T^{i|j}`

            Solve line search problem to determine step size :math:`\alpha^{i|j}`

            Update controls :math:`\mb u^{i|j+1}`, parameters :math:`\mb p^{i|j+1}` and end time :math:`T^{i|j+1}`

            **If** minimization is converged **then**
                Break inner loop

            **End If**

        **End For**

        Set :math:`\mb u^i = \mb u^{i|j+1}`, :math:`\mb p^i = \mb p^{i|j+1}` and :math:`\mb T^i = \mb T^{i|j+1}`

        Compute :math:`\mb x^i` by forward time integration of system dynamics

        Evaluate all constraints

        Update multipliers :math:`\mb{\bar{\mu}}^{i+1}` and penalties :math:`\mb{\bar{c}}^{i+1}`

        **If** minimization is converged & constraint thresholds are satisfied **then**
            Break outer loop

        **End If**

    **End For**

    Compute cost :math:`J` and norm of constraints


The basic structure of the algorithm that is implemented in the main calling function ``grampc_run`` is outlined in
:prf:ref:`alg:AlgOpt:GrampcAlgorithm`.
The projected gradient method is realized in the inner loop and consists
of the forward and backward integration of the canonical equations as
well as the update of the optimization variables based on the gradient
and the approximate solution of the line search problem. The outer loop
corresponds to the augmented Lagrangian method consisting of the
solution of the inner minimization problem and the update of the
multipliers and penalty parameters.

As an alternative to the augmented Lagrangian framework, the user can
choose external penalty functions in GRAMPC that handle the equality
and inequality constraints as “soft” constraints. In this case, the
multipliers :math:`\mb{\bar \mu}` are fixed at zero and only the
penalty parameters are updated in the outer loop. Note that the user can
set the options ``PenaltyIncreaseFactor`` and ``PenaltyDecreaseFactor`` to :math:`1.0` in order to keep the penalty
parameters at the initial value ``PenaltyMin``. The single steps of the algorithm and
the related options are described in more detail in the following sections.

The following options can be used to adjust the basic algorithm. The
corresponding default values are listed in
:numref:`tab:ListOfOptions` in the appendix.

-  ``MaxMultIter``: Sets the maximum number of augmented Lagrangian iterations
   :math:`i_\text{max} \geq 1`. If the option ``ConvergenceCheck`` is activated, the
   algorithm evaluates the convergence criterion and terminates if the
   inner minimization converged and all constraints are satisfied within
   the tolerance defined by ``ConstraintsAbsTol``.

-  ``MaxGradIter``: Sets the maximum number of gradient iterations :math:`j_\text{max} \geq 1`. If the option ``ConvergenceCheck`` is activated, the algorithm terminates the inner loop
   as soon as the convergence criterion is fulfilled.

-  ``EqualityConstraints``: Equality constraints
   :math:`\mb{g}(\mb{x}(t), \mb{u}(t), \mb{p}, t) = \mb{0}`
   can be disabled by the option value ``off``.

-  ``InequalityConstraints``: To disable inequality constraints
   :math:`\mb{h}(\mb{x}(t), \mb{u}(t), \mb{p}, t) \le \mb{0}`,
   set this option to ``off``.

-  ``TerminalEqualityConstraints``: To disable terminal equality constraints
   :math:`\mb{g}_T(\mb{x}(T), \mb{p}, T) = \mb{0}`,
   set this option to ``off``.

-  ``TerminalInequalityConstraints``: To disable terminal inequality constraints
   :math:`\mb{h}_T(\mb{x}(T), \mb{p}, T) \le \mb{0}`,
   set this option to ``off``.

-  ``ConstraintsHandling``: State constraints are handled either by means of the augmented
   Lagrangian approach (option value ``auglag``) or as soft constraints by outer
   penalty functions (option value ``extpen``).

-  ``OptimControl``: Specifies whether the cost functional should be minimized with
   respect to the control variable :math:`\mb{u}`.

-  ``OptimParam``: Specifies whether the cost functional should be minimized with
   respect to the optimization parameters :math:`\mb{p}`.

-  ``OptimTime``: Specifies whether the cost functional should be minimized with
   respect the horizon length :math:`T` (free end time problem) or if
   :math:`T` is kept constant.

.. footbibliography::