.. _sec:OptimizationProblem:

Optimization problem and parameters
-----------------------------------

GRAMPC allows one to cope with optimal control problems of the
following type

.. math::
    :label: eq:OCP
    
   	\!\!\!\min_{\mb{u}, \mb{p}, T} \quad & J(\mb{u}, \mb{p}, T;\mb{x}_0) = V(\mb{x}(T), \mb{p}, T) + \int_0^T l(\mb{x}(t), \mb{u}(t), \mb{p}, t) \, {\rm d}t

   	\textrm{s.t.} \quad & \mb{M} \mb{\dot x}(t) = \mb{f}(\mb{x}(t), \mb{u}(t), \mb{p}, t) \,, \quad \mb{x}(t_0) = \mb{x}_0
   	
   	& \mb{g}(\mb{x}(t), \mb{u}(t), \mb{p}, t) = \mb{0} \,, \quad \mb{g}_T(\mb{x}(T), \mb{p}, T) = \mb{0}
   	
   	& \mb{h}(\mb{x}(t), \mb{u}(t), \mb{p}, t) \le \mb{0} \,, \quad \mb{h}_T(\mb{x}(T), \mb{p}, T) \le \mb{0}
   	
   	& \mb{u}(t) \in \left[\mb{u}_{\min}, \mb{u}_{\max}\right]
   	
   	& \mb{p} \in \left[\mb{p}_{\min}, \mb{p}_{\max}\right] \,,\quad  T \in \left[T_{\min}, T_{\max}\right]

in the context of model predictive control, moving horizon estimation, and/or parameter estimation. 
The cost functional :math:`J(\mb{u}, \mb{p}, T;\mb{x}_0)` to be minimized consists of the continuously differentiable terminal cost (Mayer term)
:math:`V:\mathbb{R}^{N_{\mb{x}}}\times\mathbb{R}^{N_{\mb{p}}}\times\mathbb{R}\rightarrow \mathbb{R}` and integral cost (Lagrange term)
:math:`l:\mathbb{R}^{N_{\mb{x}}}\times\mathbb{R}^{N_{\mb{u}}}\times\mathbb{R}^{N_{\mb{p}}}\times\mathbb{R} \rightarrow \mathbb{R}`
with the state variables :math:`\mb{x}\in \mathbb{R}^{N_{\mb{x}}}`, the control variables :math:`\mb{u} \in \mathbb{R}^{N_{\mb{u}}}`,
the parameters :math:`\mb{p}\in\mathbb{R}^{N_{\mb{p}}}`, and the end time :math:`T\in \mathbb{R}`.

The cost functional :math:`J(\mb{u}, \mb{p}, T;\mb{x}_0, t_0)` is
minimized with respect to the optimization variables :math:`(\mb{u}, \mb{p}, T)` subject to the system
dynamics :math:`\mb{M} \mb{\dot x}(t) = \mb{f}(\mb{x}(t), \mb{u}(t), \mb{p}, t)` with the mass matrix
:math:`\mb{M}\in\mathbb{R}^{N_{\mb{x}}\times N_{\mb{x}}}`,
the continuously differentiable right hand side
:math:`\mb{f}:\mathbb{R}^{N_{\mb{x}}}\times\mathbb{R}^{N_{\mb{u}}}\times\mathbb{R}^{N_{\mb{p}}}\times\mathbb{R} \rightarrow \mathbb{R}^{N_{\mb{x}}}`, and the initial state :math:`\mb{x}_0`. 
GRAMPC allows one to formulate terminal equality and inequality constraints
:math:`\mb{g}_T:\mathbb{R}^{N_{\mb{x}}}\times\mathbb{R}^{N_{\mb{p}}}\times\mathbb{R} \rightarrow \mathbb{R}^{N_{\mb{g}_T}}`
and
:math:`\mb{h}_T:\mathbb{R}^{N_{\mb{x}}}\times\mathbb{R}^{N_{\mb{p}}}\times\mathbb{R} \rightarrow \mathbb{R}^{N_{\mb{h}_T}}`
as well as general equality and inequality constraints
:math:`\mb{g}:\mathbb{R}^{N_{\mb{x}}}\times\mathbb{R}^{N_{\mb{u}}}\times\mathbb{R}^{N_{\mb{p}}}\times\mathbb{R} \rightarrow \mathbb{R}^{N_{\mb{g}}}`
and :math:`\mb{h}:\mathbb{R}^{N_{\mb{x}}}\times\mathbb{R}^{N_{\mb{u}}}\times\mathbb{R}^{N_{\mb{p}}}\times\mathbb{R} \rightarrow \mathbb{R}^{N_{\mb{h}}}`.
In addition, the optimization variables are limited by the box
constraints :math:`\mb{u}(t) \in \left[\mb{u}_{\min}, \mb{u}_{\max}\right]`, :math:`\mb{p} \in \left[\mb{p}_{\min}, \mb{p}_{\max}\right]` and :math:`T \in \left[T_{\min}, T_{\max}\right]`

The terminal cost :math:`V`, the integral cost :math:`l`, as well as the system dynamics :math:`\mb{f}` and all constraints :math:`(\mb{g}, \mb{g}_\text{T}, \mb{h}, \mb{h}_\text{T})` contain an explicit time dependency with regard to the internal time :math:`t \in [0, T]`.
In the context of MPC, the internal time is distinguished from the global time :math:`t_0+t \in [t_0, t_0+T]` where the initial time :math:`t_0` and initial state :math:`\mb{x}_0` correspond to the sampling instant :math:`t_k` that is incremented by the sampling time :math:`\Delta t>0` in each MPC step :math:`k`.

A detailed description about the implementation of the optimization problem :math:numref:`eq:OCP` in C code is given in :ref:`sec:ProblemImplementation`. 
In addition, some parts of the problem can be configured by parameters (cf. the GRAMPC data structure ``param``) and therefore do not require repeated compiling. 
A list of all parameters with types and allowed values is provided in the appendix (:numref:`tab:ListOfParameters`). 
Except for the horizon length ``Thor`` and the sampling time ``dt``, all parameters are optional and initialized to default values. 
A description of all parameters is as follows:

-  ``x0``: Initial state vector :math:`\mb{x}(t_0)=\mb{x}_0` at the corresponding sampling time :math:`t_0`.

-  ``xdes``: Desired (constant) setpoint vector for the state variables :math:`\mb{x}`.

-  ``u0``: Initial value of the control vector :math:`\mb{u}(t) = \mb{u}_0 = \text{const.}`, :math:`t \in [0,T]` that is used in the first iteration of GRAMPC.

-  ``udes``: Desired (constant) setpoint vector for the control variables :math:`\mb{u}`.

-  ``umin``, ``umax`` : Lower and upper bounds for the control variables :math:`\mb{u}`.

-  ``p0``: Initial value of the parameter vector :math:`\mb{p} = \mb{p}_0` that is used in the first iteration of GRAMPC.

-  ``pmin``, ``pmax``: Lower and upper bounds for the parameters :math:`\mb{p}`.

-  ``Thor``: Prediction horizon :math:`T` or initial value if the end time is optimized.

-  ``Tmin``, ``Tmax``: Lower and upper bound for the prediction horizon :math:`T`.

-  ``dt``: Sampling time :math:`\Delta t` of the considered system for model predictive control or moving horizon estimation. Required for prediction of next state ``grampc.sol.xnext`` and for the control shift, see :ref:`sec:AlgOpt:ControlShift`.

-  ``t0``: Current sampling instance :math:`t_0` that is provided in the ``grampc.param`` structure.

-  ``userparam``: Further problem-specific parameters, e.g. system parameters or weights in the cost functions that are passed to the problem functions via a ``void``-pointer in C or ``typeRNum`` array in MATLAB.

Although GRAMPC uses a continuous-time formulation of the optimization problem :math:numref:`eq:OCP`, all trajectories are internally stored in discretized form with ``Nhor`` steps (cf. :ref:`sec:AlgOpt:Integration`). 
This raises the question of whether all constraints are evaluated for the last trajectory point or only the terminal ones. In general, the constraints should be formulated in such a way that there are no conflicts. 
However, numerical difficulties can arise in some problems if constraints are formulated twice for the last point. Therefore, GRAMPC does not evaluate the constraints :math:`\mb{g}` and :math:`\mb{h}` for the last trajectory point if terminal constraints are defined, i.e. :math:`N_{\mb{g}_T} + N_{\mb{h}_T} > 0`. 
In contrast, if no terminal constraints are defined, the functions :math:`\mb{g}` and :math:`\mb{h}` are evaluated for all points. 
Note that the opposite behavior is easy to implement by including :math:`\mb{g}` and :math:`\mb{h}` in the terminal constraints :math:`\mb{g}_T` and :math:`\mb{h}_T`.
