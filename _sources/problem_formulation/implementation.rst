.. _sec:ProblemImplementation:

Problem implementation
----------------------

Regardless of what kind of problem statement is considered (i.e. model predictive control, moving horizon estimation or parameter estimation), the optimization problem :math:numref:`eq:OCP` must be implemented in C, cf. :numref:`fig:grampcGeneralStructure` 
For this purpose, the C file template ``probfct_TEMPLATE.c`` is provided within the folder ``<grampc_root>/examples/TEMPLATE``, which allows one to describe the structure of the optimization problem :math:numref:`eq:OCP`, the cost functional :math:`J` to be minimized, the system dynamics :math:`f`, and the constraints :math:`\mb{g}, \mb{g}_\text{T}, \mb{h}, \mb{h}_\text{T}`.
The number of C functions to be provided depends on the type of dynamics :math:`f` of the specific problem at hand.

.. _sec:ProblemImplementation:Explicit:

Problems involving explicit ODEs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A special case of the system dynamics :math:`f` and certainly the most important one concerns ordinary differential equations (ODEs) with explicit appearance of the first-order derivatives

.. math::

   \mb{\dot x}(t) = \mb{f}(\mb{x}(t), \mb{u}(t), \mb{p}, t)\,,

corresponding to the identity mass matrix :math:`\mb{M}=\mb{I}` in OCP :math:numref:`eq:OCP`. 
In this case, the following functions of the C template file ``probfct_TEMPLATE.c`` have to be provided:

-  ``ocp_dim``: Definition of the dimensions of the considered problem, i.e. the
   number of state variables :math:`{N_{\mb{x}}}`, control
   variables :math:`{N_{\mb{u}}}`, parameters
   :math:`{N_{\mb{p}}}`, equality constraints
   :math:`{N_{\mb{g}}}`, inequality constraints
   :math:`{N_{\mb{h}}}`, terminal equality constraints
   :math:`{N_{{\mb{g}}_T}}`, and terminal inequality constraints
   :math:`N_{{\mb{h}}_T}`.

-  ``ffct``: Formulation of the system dynamics function
   :math:`\mb{f}`.

-  ``dfdx_vec``, ``dfdu_vec``, ``dfdp_vec``: Matrix vector products
   :math:`(\frac{\partial^{} \mb{f}}{\partial \mb{x}^{}})^\mathsf{T}\mb{v}`,
   :math:`(\frac{\partial^{} \mb{f}}{\partial \mb{u}^{}})^\mathsf{T}\mb{v}`
   and
   :math:`(\frac{\partial^{} \mb{f}}{\partial \mb{p}^{}})^\mathsf{T}\mb{v}`
   for the system dynamics with an arbitrary vector
   :math:`\mb{v}` of dimension :math:`N_x`.

-  ``Vfct``, ``lfct``: Terminal and integral cost functions :math:`V` and :math:`l` of
   the cost functional :math:`J`.

-  ``dVdx``, ``dVdp``, ``dVdT``, ``dldx``, ``dldu``, ``dldp``: Gradients
   :math:`\frac{\partial^{} V}{\partial \mb{x}^{}}`,
   :math:`\frac{\partial^{} V}{\partial \mb{p}^{}}`,
   :math:`\frac{\partial^{} V}{\partial T^{}}`,
   :math:`\frac{\partial^{} l}{\partial \mb{x}^{}}`,
   :math:`\frac{\partial^{} l}{\partial \mb{u}^{}}`, and
   :math:`\frac{\partial^{} l}{\partial \mb{p}^{}}` of the cost
   functions ``Vfct`` and ``lfct``.

-  ``hfct``, ``gfct``: Inequality and equality constraint functions
   :math:`\mb{h}` and :math:`\mb{g}` as defined in
   :math:numref:`eq:OCP`.

-  ``hTfct``, ``gTfct``: Terminal inequality and equality constraint functions
   :math:`\mb{h}_T` and :math:`\mb{g}_T` as defined in
   :math:numref:`eq:OCP`.

-  ``dhdx_vec``, ``dhdu_vec``, ``dhdp_vec``: Matrix products
   :math:`(\frac{\partial^{} \mb{h}}{\partial \mb{x}^{}})^\mathsf{T}\mb{v}`,
   :math:`(\frac{\partial^{} \mb{h}}{\partial \mb{u}^{}})^\mathsf{T}\mb{v}`,
   and
   :math:`(\frac{\partial^{} \mb{h}}{\partial \mb{p}^{}})^\mathsf{T}\mb{v}`
   for the inequality constraints with an arbitrary vector
   :math:`\mb{v}` of dimension :math:`N_h`.

-  ``dgdx_vec``, ``dgdu_vec``, ``dgdp_vec``: Matrix product functions
   :math:`(\frac{\partial^{} \mb{g}}{\partial \mb{x}^{}})^\mathsf{T}\mb{v}`,
   :math:`(\frac{\partial^{} \mb{g}}{\partial \mb{u}^{}})^\mathsf{T}\mb{v}`,
   and
   :math:`(\frac{\partial^{} \mb{g}}{\partial \mb{p}^{}})^\mathsf{T}\mb{v}`
   for the equality constraints with an arbitrary vector
   :math:`\mb{v}` of dimension :math:`N_g`.

-  ``dhTdx_vec``, ``dhTdu_vec``, ``dhTdp_vec``: Matrix product functions
   :math:`(\frac{\partial^{} \mb{h}_T}{\partial \mb{x}^{}})^\mathsf{T}\mb{v}`,
   :math:`(\frac{\partial^{} \mb{h}_T}{\partial \mb{p}^{}})^\mathsf{T}\mb{v}`,
   and
   :math:`(\frac{\partial^{} \mb{h}_T}{\partial T^{}})^\mathsf{T}\mb{v}`
   for the terminal inequality constraints with an arbitrary vector
   :math:`\mb{v}` of dimension :math:`N_{h_T}`.

-  ``dgTdx_vec``, ``dgTdu_vec``, ``dgTdp_vec``: Matrix product functions
   :math:`(\frac{\partial^{} \mb{g}_T}{\partial \mb{x}^{}})^\mathsf{T}\mb{v}`,
   :math:`(\frac{\partial^{} \mb{g}_T}{\partial \mb{p}^{}})^\mathsf{T}\mb{v}`,
   and
   :math:`(\frac{\partial^{} \mb{g}_T}{\partial T^{}})^\mathsf{T}\mb{v}`
   for the terminal equality constraints with an arbitrary vector
   :math:`\mb{v}` of dimension :math:`N_{g_T}`.

The respective problem function templates only have to be filled in if the corresponding constraints and cost functions are defined for the problem at hand and depending on the actual choice of optimization variables (:math:`{\mb{u}}`, :math:`{\mb{p}}`, and/or :math:`{T}`). 
For example, if only the control :math:`\mb{u}` is optimized, the partial derivatives with respect to :math:`\mb{p}` and :math:`T` are not required.

The gradients :math:`\frac{\partial^{} V}{\partial \mb{x}^{}}`,
:math:`\frac{\partial^{} V}{\partial \mb{p}^{}}`,
:math:`\frac{\partial^{} V}{\partial T^{}}`,
:math:`\frac{\partial^{} l}{\partial \mb{x}^{}}`,
:math:`\frac{\partial^{} l}{\partial \mb{u}^{}}`, and
:math:`\frac{\partial^{} l}{\partial \mb{p}^{}}` as well as the
matrix product functions listed above appear in the partial derivatives
:math:`H_{\mb{x}}`, :math:`H_{\mb{u}}` and
:math:`H_{\mb{p}}` of the Hamiltonian :math:`H` within the
gradient method (see :ref:`sec:AlgOpt:ProjGrad`). The
matrix product formulation is chosen over the definition of Jacobians
(e.g. :math:`(\frac{\partial^{} \mb{f}}{\partial \mb{x}^{}})^\mathsf{T}\mb{v}` instead of :math:`\frac{\partial^{} \mb{f}}{\partial \mb{x}^{}}`)
in order to avoid unnecessary zero multiplications for sparse matrices or alternatively the usage of sparse numerics.

.. _`sec:ProblemImplementation:Discrete`:

Problems with discrete-time systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. versionadded:: v2.3

GRAMPC also allows to consider discrete-time system dynamics of the form

.. math::

   	\mb{x}_{k+1} = \mb{f}(\mb{x}_k, \mb{u}_k, \mb{p}, t_k)

with :math:`\mb{x}_k = \mb x(t_k)` and :math:`\mb{u}_k = \mb u(t_k)`.
For this case the same functions as in :ref:`sec:ProblemImplementation:Explicit` are used to implement the problem description but the option ``Integrator`` is set to ``discrete``.
Since the discrete ``ffct`` computes the system dynamics for a fixed sampling time :math:`\Delta t` (Parameter ``dt``), the end time :math:`T` (Parameter ``Thor``) of the optimization problem must satisfy

.. math::
   
      T = (N_{\text{hor}} - 1) \Delta t

with the number of discretization points along the horizon :math:`N_{\text{hor}}` (Option ``Nhor``).
As a consequence it is no longer possible to consider a free end time, where :math:`T` is an optimization variable, since this would change the sampling time :math:`\Delta t`.
Note that for continuous-time system dynamics the sampling time :math:`\Delta t` can be chosen independently of the end time :math:`T` and the number of discretizaton points :math:`N_{\text{hor}}`, because it is only used for predicting the next state ``grampc.sol.xnext`` and for the control shift.

.. note::

      For discrete-time systems the option ``Integrator`` is set to ``discrete``, the settings ``Nhor``, ``Thor`` and ``dt`` need to satisfy the relation ``Thor = (Nhor-1)*dt`` and the option ``OptimTime`` must be set to ``off``.

An MPC example that compares continuous-time and discrete-time versions of a double integrator is included in ``<grampc_root>/examples/Continuous_vs_Discrete``.
Furthermore, a discrete-time formulation of the helicopter example is provided in ``<grampc_root>/examples/Continuous_vs_Discrete``.

.. _`sec:ProblemImplementation:Implicit`:

Problems involving semi-implicit ODEs and DAEs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Beside explicit ODEs, GRAMPC supports semi-implicit ODEs with mass
matrix :math:`\mb{M} \neq \mb{I}` and DAEs with
:math:`\mb{M}` being singular. The underlying numerical
integrations of the dynamics :math:`f` is
carried out using the integrator
RODAS :footcite:`Hairer:Book:1996:Stiff, Rodas:Webpage:2018`. In
this case, additional C functions must be provided and several
RODAS-specific options must be set, cf. :ref:`sec:AlgOpt:IntegrationRodas` and :ref:`sec:setting_opt_par`. 
Especially, the option
``IMAS = 1`` must be set to indicate that a mass matrix is given. The
numerical integrations performed with RODAS can be accelerated by
providing partial derivatives. In summary, the following additional C
functions are used by GRAMPC for semi-implicit ODEs and DAE systems:

-  ``Mfct``, ``Mtrans``: Definition of the mass matrix :math:`\mb{M}` and its
   transpose :math:`\mb{M}^\mathsf{T}`, which is required for
   the adjoint dynamics, cf. :ref:`sec:AlgOpt:ProjGrad` in
   the projected gradient algorithm. The matrices must be specified
   column-wise. If the mass matrix has a band structure, only the
   respective elements above and below the main diagonal are specified.
   This only applies if the options ``IMAS = 1`` and ``MLJAC < N_x`` are
   selected. Non-existent elements above or below the main diagonal must
   be filled with zeros so that the same number of elements is specified
   for each column.

-  ``dfdx``, ``dfdxtrans``: The Jacobians
   :math:`\frac{\partial^{} \mb{f}}{\partial \mb{x}^{}}` and
   :math:`(\frac{\partial^{} \mb{f}}{\partial \mb{x}^{}})^\mathsf{T}=\frac{\partial^{2} H}{\partial \mb{x}\partial \mb{\lambda}}`
   are provided by these functions if the option ``IJAC = 1`` is set.
   This allows one to evaluate the right hand sides of the canonical
   equations time efficiently. The Jacobians must be implemented in
   vector form by arranging the successive columns for
   :math:`\frac{\partial^{} \mb{f}}{\partial \mb{x}^{}}` and
   :math:`(\frac{\partial^{} \mb{f}}{\partial \mb{x}^{}})^\mathsf{T}`.
   If the option ``MLJAC < N_x`` is set to exploit the band structure of
   the Jacobians, only the corresponding elements above and below the
   main diagonal must be specified.

-  ``dfdt``, ``dHdxdt``: The partial derivatives
   :math:`\frac{\partial^{} \mb{f}}{\partial t^{}}` and
   :math:`\frac{\partial^{2} H}{\partial \mb{x}\partial t}`
   allow for evaluating the right hand sides of the canonical equations
   time efficiently, if the problem explicitly depends on time
   :math:`t`. These functions must only be provided if the options
   ``IFCN = 1`` and ``IDFX = 1`` are used.

An MPC example with a semi-implicit system dynamics is included in ``<grampc_root>/examples/Reactor_PDE``.
The problem formulation is derived from a quasi-linear
diffusion-convection-reaction system that is spatially discretized using
finite elements.

.. _`sec:finite_diff`:

Finite differences and gradient check
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. versionadded:: v2.3

For the purpose of rapid prototyping, the partial derivatives can also be approximated by forward finite differences.
To this end, the header ``finite_diff.h`` provides a helper function

.. code-block:: c

    void finite_diff_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam,
                         typeRNum *memory, ctypeRNum step_size, const typeFiniteTarget target, ctypeInt func_out_size, const typeFiniteDiffFctPtr func);

where the first eight arguments are the same as for typical functions in ``probfct.h`` and the other arguments are:

- ``memory``: Sufficient user-provided memory of size at most :math:`3 \cdot \max\{ N_x, N_u, N_p, N_g, N_h, N_{g_T}, N_{h_T} \}` that is needed for storing intermediate values.
  If used within a probfct, this memory should be provided via the ``userparam`` pointer to avoid dynamic memory allocation during the optimization.

- ``step_size``: Small value used as step size for the finite differences.

- ``target``: Argument with respect to which the finite differences are computed. Options are ``DX`` for states, ``DU`` for inputs, ``DP`` for parameters and ``DT`` for time.

- ``func_out_size``: Size of the function's ``out`` argument, e.g. ``Nh`` for ``hfct``.

- ``func``: Pointer to the function that is differentiated. For ``Vfct``, ``gTfct`` and ``hTfct`` wrappers are provided with an unused dummy argument in place of the control ``u``.

The intended usage is shown in the template file ``probfct_TEMPLATE_finite_diff.c``.
In addition, the helper function

.. code-block:: C

    void finite_diff_dfdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam,
                          typeRNum *memory, ctypeRNum step_size, ctypeInt ml, ctypeInt mu, typeBoolean transpose);

approximates the partial derivatives of the ``ffct`` for semi-implicit ODEs and DAEs with the arguments:

- ``memory``: Sufficient allocated memory for storing intermediate values of size :math:`3 \cdot N_x`. If used within a probfct, this memory should be provided via the ``userparam`` pointer to avoid dynamic memory allocation during the optimization.

- ``step_size``: Small value used as step size for the finite differences.

- ``ml``: Number of lower non-zero diagonals if banded. Set to ``Nx`` for full matrix.

- ``mu``: Number of upper non-zero diagonals if banded. Set to ``Nx`` for full matrix.

- ``transpose``: Set to ``0`` for approximating ``dfdx`` and to ``1`` for approximating ``dfdxtrans``.

Finally, it is recommended to validate the user-supplied analytic gradients by comparison to finite differences.
For this task, GRAMPC provides a helper function

.. code-block:: C

    void grampc_check_gradients(typeGRAMPC *grampc, ctypeRNum tolerance, ctypeRNum step_size);

that checks the gradients at the initial point defined by ``t0``, ``x0``, ``u0`` and ``p0``.
The tolerance is applied element-wise to the relative error

.. math::

    \left| \frac{\nabla \mb{f}_{fd,i} - \nabla \mb{f}_i}{\max(1, \nabla \mb{f}_i)} \right|


where :math:`\nabla \mb{f}_{fd,i}` is the finite difference approximation and :math:`\nabla \mb{f}_i` the user-supplied derivative.
If the tolerance is violated, the function prints messages like

::

    dfdx_vec: element (out_index,vec_index) exceeds supplied tolerance with 4.234546e-5

so that one is able to correct the wrong derivative.
Recommended settings for the gradient check include a step size of :math:`\sqrt{\text{eps}}` where *eps* is the floating point machine precision, and a tolerance of 1e-6.
Nevertheless, the gradient checker can produce false positives, so every flagged gradient should be carefully checked.


.. _`sec:ball-on-plate`:

Example: Ball-on-plate system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The appropriate definition of the C functions of the template
``probfct_TEMPLATE`` is described for the example of a ball-on-plate
system :footcite:`Richter2012` in the context of MPC. 
The problem is also included in ``<grampc_root>/examples(BallOnPlate)``. 
The underlying optimization problem reads as

.. math::
    :label: BallOnPlate
    
    % cost
    \min_{u(\cdot)} \quad &J(u;\mb{x}_0) = \frac{1}{2}\Delta\mb{x}^\mathsf{T}(T) \mb{P} \Delta \mb{x}(T)+\frac{1}{2}\int_{0}^{T} \Delta \mb{x}^\mathsf{T}(\tau) \mb{Q} \Delta\mb{x}(\tau) + R \Delta u^2 \,{\rm d}\tau \\
    % dynamics and intial condition
    \text{s.t.} \quad & \begin{bmatrix} \dot x_1 \\ \dot x_2\end{bmatrix} 
    = \begin{bmatrix} 0 & 1 \\ 0 & 0\end{bmatrix}
    \begin{bmatrix} x_1\\ x_2\end{bmatrix} +
    \begin{bmatrix} -0.04 \\ -7.01\end{bmatrix}u \,, \quad \begin{bmatrix} x_1(0)\\ x_2(0)\end{bmatrix} = \begin{bmatrix} x_{k,1}\\ x_{k,2}\end{bmatrix}\\
    % state constraints
    & \begin{bmatrix} -0.2\\-0.1\end{bmatrix}\le\begin{bmatrix} x_1\\ x_2 \end{bmatrix}\le\begin{bmatrix} 0.01\\0.1\end{bmatrix} \,, \quad |u| \le 0.0524

The cost functional in :math:numref:`BallOnPlate` penalizes the state and input error :math:`\Delta \mb{x}=\mb{x}-\mb{x}_\text{des}` and :math:`\Delta u=u-u_\text{des}` in a quadratic manner using the weights

.. math::
   \mb{P}=\mb{Q} = \begin{bmatrix} 100 & 0 \\ 0 & 10\end{bmatrix}, \quad
   R = 1.

The system dynamics in :math:numref:`BallOnPlate`
describes a simplified linear model of a single axis of a ball-on-plate
system :footcite:`Richter2012`. An optimal solution of the
optimization problem has to satisfy the input and state
constraints present in :math:numref:`BallOnPlate`.

The user must provide the C functions and to describe the general
structure and the system dynamics of the optimization problem, also see :ref:`sec:ProblemImplementation`. 
As shown in :numref:`lis:GenStruct`, the C function ``ocp_dim`` is used to
define the number of states, control inputs, and number of (terminal)
inequality and equality constraints. Note that GRAMPC uses the generic
type ``typeInt`` for integer values. The word size of this integer type can be
changed in the header file ``grampc_macro.h`` within the folder ``<grampc_root>/include``. 
This is particularly advantageous with regard to implementing GRAMPC on embedded hardware.

.. code-block:: c
    :caption: Settings of the general structure of optimization problem :math:numref:`BallOnPlate`.
    :name: lis:GenStruct

    /** OCP dimensions **/
    void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, 
                 typeInt *NgT, typeInt *NhT, typeUSERPARAM *userparam)
    { 
        *Nx = 2;
        *Nu = 1;
        *Np = 0;
        *Nh = 4;
        *Ng = 0;
        *NhT = 0;	
        *NgT = 0;
    }

The system dynamics are described by the C function ``ffct`` shown in
:numref:`lis:SysDyn`. The example is given in explicit
ODE form, for which the functions ``Mfct`` and ``Mtrans`` for the mass matrix :math:`M` are
not required. Similar to the generic integer type ``typeInt``, the data type ``typeRNum`` is
used to adress floating point numbers of different word sizes, e.g.
float or double (cf. the header file ``grampc_macro.h`` included in the folder ``<grampc_root>/include``).

.. code-block:: c
    :caption: Formulation of the system dynamics of optimization problem :math:numref:`BallOnPlate`.
    :name: lis:SysDyn

    /** System function f(t,x,u,p,grampc.param,userparam) **/
    void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, 
              ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
    {
        out[0] = x[1]-0.04*u[0];
        out[1] =     -7.01*u[0];
    }

The cost functions are defined via the functions ``lfct`` and ``Vfct``, cf. :numref:`lis:Cost`.
Note that the input argument ``userparam`` is used to parametrize the cost functional in a generic way.

.. code-block:: c
    :caption: Formulation of the cost functional of optimization problem :math:numref:`BallOnPlate`.
    :name: lis:Cost

    /** Integral cost l(t,x(t),u(t),p,grampc.param,userparam) **/
    void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, 
              ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
    {
        ctypeRNum* pCost = (ctypeRNum*)userparam;
        ctypeRNum* xdes = param->xdes;
        ctypeRNum* udes = param->udes;
        out[0] = (pCost[0] * (x[0] - xdes[0]) * (x[0] - xdes[0])
                + pCost[1] * (x[1] - xdes[1]) * (x[1] - xdes[1])
                + pCost[2] * (u[0] - udes[0]) * (u[0] - udes[0])) / 2;
    }

    /** Terminal cost V(T,x(T),p,grampc.param,userparam) **/
    void Vfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, 
              const typeGRAMPCparam *param, typeUSERPARAM *userparam)
    {
        ctypeRNum* pCost = (ctypeRNum*)userparam;
        ctypeRNum* xdes = param->xdes;
        out[0] = (pCost[3] * (x[0] - xdes[0]) * (x[0] - xdes[0])
                + pCost[4] * (x[1] - xdes[1]) * (x[1] - xdes[1])) / 2;
    }

:numref:`lis:Constraints` shows the formulation
of the inequality constraints
:math:`\mb{h}(\mb{x}(t), \mb{u}(t), \mb{p}, t) \le \mb{0}`.
For the sake of completeness,
:numref:`lis:Constraints` also contains the
corresponding functions for equality constraints
:math:`\mb{g}(\mb{x}(t), \mb{u}(t), \mb{p}, t) = \mb{0}`,
terminal inequality constraints
:math:`\mb{h}_T(\mb{x}(T), \mb{p}, T) \le \mb{0}`
as well as terminal equality constraints
:math:`\mb{g}_T(\mb{x}(T), \mb{p}, T) = \mb{0}`,
which are not defined for the ball-on-plate example. Similar to the
formulation of the cost functional, the input argument ``userparam`` is used to
parametrize the inequality constraints.

.. code-block:: c
    :caption: Formulation of the state constraints of optimization problem :math:numref:`BallOnPlate`.
    :name: lis:Constraints

    /** Inequality constraints h(t,x(t),u(t),p,grampc.param,userparam) <= 0 **/
    void hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, 
              ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
    {
        ctypeRNum* pSys = (ctypeRNum*)userparam;

        out[0] =  pSys[5] - x[0];
        out[1] = -pSys[6] + x[0];
        out[2] =  pSys[7] - x[1];
        out[3] = -pSys[8] + x[1];
    }

    /** Equality constraints g(t,x(t),u(t),p,grampc.param,userparam) = 0 **/
    void gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, 
              ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
    {
    }

    /** Terminal inequality constraints hT(T,x(T),p,grampc.param,userparam) <= 0 **/
    void hTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, 
               const typeGRAMPCparam *param, typeUSERPARAM *userparam)
    {
    }

    /** Terminal equality constraints gT(T,x(T),p,grampc.param,userparam) = 0 **/
    void gTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, 
               const typeGRAMPCparam *param, typeUSERPARAM *userparam)
    {
    }

The Jacobians of the single functions defined in
:numref:`lis:SysDyn` to :numref:`lis:Constraints`
with respect to state :math:`\mb{x}` and control
:math:`\mb{u}` are required for evaluating the optimality
conditions of optimization problem
:math:numref:`BallOnPlate` within the gradient
algorithm, see :ref:`sec:AlgOpt:ProjGrad`. If applicable,
the Jacobians of the above-mentioned functions are also required with
respect to the optimization variables :math:`\mb{p}` and
:math:`T`. :numref:`lis:SysDynJacx` shows the
corresponding Jacobians for the ball-on-plate example. For the matrix
product functions ``dfdx_vec``, ``dfdu_vec``, and ``dhdx_vec``, the pointer to a generic vector ``vec`` is
passed as input argument that corresponds to the adjoint state,
respectively a vector that accounts for the Lagrange multiplier and
penalty term of state constraints,
cf. :ref:`sec:AlgOpt:BasicAlgorithm`. Note that ``vec`` is of
appropriate dimension for the respective matrix product function,
i.e. of dimension :math:`N_{\mb{x}}` or
:math:`N_{\mb{h}}`. A complete C function template and further
examples concerning the problem formulation are included in the GRAMPC
software package.

.. code-block:: c
    :caption: Jacobians of the system dynamics and inequality constraint.
    :name: lis:SysDynJacx

    /** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
    void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, 
                  ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
    {
        out[0] = 0;
        out[1] = vec[0];
    }

    /** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
    void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, 
                  ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
    {
        out[0] = (typeRNum)(-0.04)*vec[0] - (typeRNum)(7.01)*vec[1];
    }


    /** Gradient dl/dx **/
    void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, 
              const typeGRAMPCparam *param, typeUSERPARAM *userparam)
    {
        ctypeRNum* pCost = (ctypeRNum*)userparam;
        ctypeRNum* xdes = param->xdes;
        
        out[0] = pCost[0] * (x[0] - xdes[0]);
        out[1] = pCost[1] * (x[1] - xdes[1]);
    }
    /** Gradient dl/du **/
    void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, 
              const typeGRAMPCparam *param, typeUSERPARAM *userparam)
    {
        ctypeRNum* pCost = (ctypeRNum*)userparam;
        ctypeRNum* udes = param->udes;
        
        out[0] = pCost[2] * (u[0] - udes[0]);
    }

    /** Gradient dV/dx **/
    void dVdx(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, 
              const typeGRAMPCparam *param, typeUSERPARAM *userparam)
    {
        ctypeRNum* pCost = (ctypeRNum*)userparam;
        ctypeRNum* xdes = param->xdes;
        
        out[0] = pCost[3] * (x[0] - xdes[0]);
        out[1] = pCost[4] * (x[1] - xdes[1]);
    }

    /** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
    void dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, 
                  ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
    {
        out[0] = -vec[0] + vec[1];
        out[1] = -vec[2] + vec[3];
    }
    ...


.. footbibliography::