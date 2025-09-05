.. _sec:AlgOpt:Integration:

Numerical Integration
---------------------

GRAMPC employs a continuous-time formulation of the optimization
problem. However, internally, all time-dependent functions are stored in
discretized form with :math:`N_\text{hor}` elements and numerical
integration is performed to compute the cost functional and to solve the
differential equations.

.. _sec:AlgOpt:IntegrationCostODE:

Integration of cost functional and explicit ODEs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The approximate line search method requires the evaluation of the cost
functional. To this end, the integral cost is integrated numerically
with the trapezoidal rule, Simpson rule or discrete summation and the
terminal cost is added. Note that the discrete summation does not
compute the Riemann sum but the plain sum with
:math:`\sum_k l(\mb{x}_k, \mb{u}_k, \mb{p}, t_k)`.
The cost values are additionally evaluated at the end of the
optimization as an add-on information for the user.

The gradient algorithm involves the sequential forward integration of
the system dynamics and backward integration of the adjoint dynamics.
To this end, GRAMPC provides the following four explicit Runge-Kutta (ERK) methods with fixed step size and the Butcher tableaus

.. math::
   \text{erk1 (Euler's first-order method)}: \quad
   \begin{array}{c|c}
      0 & \\
      \hline
      & 1 \\
   \end{array}
   
.. math::
   \text{erk2 (Heun's second-order method)}: \quad
   \begin{array}{c|cc}
      0 & & \\
      1 & 1 & \\
      \hline
      & 1/2 & 1/2 \\
   \end{array}
   
.. math::
   \text{erk3 (Kutta's third-order method)}: \quad
   \begin{array}{c|ccc}
      0 & & & \\
      1/2 & 1/2 & & \\
      1 & -1 & 2 & \\
      \hline
      & 1/6 & 2/3 & 1/6 \\
   \end{array}
   
.. math::
   \text{erk4 (Kutta's fourth-order method)}: \quad
   \begin{array}{c|cccc}
      0 & & & & \\
      1/2 & 1/2 & & & \\
      1/2 & 0 & 1/2 & & \\
      1 & 0 & 0 & 1 & \\
      \hline
      & 1/6 & 1/3 & 1/3 & 1/6 \\
   \end{array}

In addition, a 4th-order explicit Runge-Kutta method with variable step size is available (option value ``ruku45``).
For discrete-time system dynamics, the integrator is set to ``discrete``.
For semi-implicit ODEs and DAEs, the Rosenbrock solver RODAS :footcite:`Hairer:Book:1996:Stiff,Rodas:Webpage:2018` with variable step size has to be used (option value ``rodas``).

.. note::

    For discrete-time systems the option ``Integrator`` is set to ``discrete``, the settings ``Nhor``, ``Thor`` and ``dt`` need to satisfy the relation ``Thor = (Nhor-1)*dt`` and the option ``OptimTime`` must be set to ``off``.
    See also :ref:`sec:ProblemImplementation:Discrete`.

The following options can be used to adjust the numerical integrations:

-  ``Nhor``: Number of discretization points within the time interval :math:`[0,T]`.

-  ``IntegralCost``, ``TerminalCost``: Indicate if the integral and/or terminal cost functions are defined.

-  ``IntegratorCost``: This option specifies the integration scheme for the cost
   functionals. Possible values are ``trapezoidal``, ``simpson`` and ``discrete``.

.. versionchanged:: v2.3

    Fixed typo in ``trapezoidal`` option.
    Renamed ``euler`` to ``erk1`` and ``heun`` to ``erk2``.
    Removed ``modeuler``.
    Added ``erk3`` and ``erk4``.

-  ``Integrator``: This option specifies the integration scheme for the system and adjoint dynamics.
   Possible values are ``erk1``, ``erk2``, ``erk3``, ``erk4``, ``discrete`` with fixed step size and ``ruku45``, ``rodas`` with variable step size.

-  ``IntegratorMinStepSize``: Minimum step size for RODAS and the Runge-Kutta integrator.

-  ``IntegratorMaxSteps``: Maximum number of steps for RODAS and the Runge-Kutta integrator.

-  ``IntegratorRelTol``: Relative tolerance for RODAS and the Runge-Kutta integrator with
   variable step size. Note that this option may be insignificant if the
   minimum step size is chosen too high or the maximum number of steps
   is set too low.

-  ``IntegratorAbsTol``: Absolute tolerance for RODAS and the Runge-Kutta integrator with
   variable step size. Note that this option may be insignificant if the
   minimum step size is chosen too high or the maximum number of steps
   is set too low.

.. _sec:AlgOpt:IntegrationRodas:

Integration of semi-implicit ODEs and DAEs (RODAS)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GRAMPC supports problem descriptions with ordinary differential
equations in semi-implicit form as well as differential algebraic
equations using the solver RODAS for numerical integration (see
:ref:`sec:ProblemImplementation:Implicit`). If a
semi-implicit problem or differential algebraic equations are
considered, the mass matrix :math:`\mb{M}` and its transposed
version :math:`\mb{M}^\mathsf{T}` must be defined by the C
functions ``Mfct`` and ``Mtrans``. The numerical integration can be accelerated by
additionally providing the C functions ``dfdx``, ``dfdxtrans``, ``dfdt`` and ``dHdxdt``. See
:ref:`sec:ProblemImplementation` for a detailed description of the problem implementation.

The integration with RODAS is configured by a number of flags that are
passed to the solver using the vector ``FlagsRodas`` with the elements
``[IFCN, IDFX, IJAC, IMAS, MLJAC, MUJAC, MLMAS, MUMAS]``.
See :footcite:`Rodas:Webpage:2018,Hairer:Book:1996:Stiff` for a
detailed description of these flags. The default values
:math:`[0,0,0,0,N_x,N_x,N_x,N_x]` correspond to an autonomous system
with an identity matrix as mass matrix. The following options can be
adjusted via ``FlagsRodas``:

-  ``IFCN``: Specifies if the right hand side of the system dynamics
   :math:`\mb{f} (\mb{x},\mb{u},\mb{p},t)`
   explicitly depends on time :math:`t` (``IFCN`` = :math:`1`) or if
   the problem is autonomous (``IFCN`` = :math:`0`).

-  ``IDFX``: Specifies how the computation of the partial derivatives
   :math:`\frac{\partial^{} \mb{f}}{\partial t^{}}` and
   :math:`\frac{\partial^{2} H}{\partial x\partial t}` is carried out.
   The partial derivatives are computed internally by finite differences
   (``IDFX`` = :math:`0`) or are provided by the functions ``dfdt`` and ``dHdxdt``
   (``IDFX`` = :math:`1`) as described in :ref:`sec:ProblemImplementation`.

-  ``IJAC``: Specifies how the computation of the Jacobians
   :math:`\frac{\partial^{} \mb{f}}{\partial \mb{x}^{}}`
   and
   :math:`(\frac{\partial^{} \mb{f}}{\partial \mb{x}^{}})^\mathsf{T}=\frac{\partial^{2} H}{\partial x\partial \lambda}`
   is carried out for numerically solving the canonical equations. The
   Jacobians are computed internally by finite differences
   (``IJAC`` = :math:`0`) or are provided by the functions ``dfdx`` and ``dfdxtrans``
   (``IJAC`` = :math:`1`), also see :ref:`sec:ProblemImplementation`.

-  ``IMAS``: Gives information on the mass matrix :math:`\mb{M}`, which
   is either an identity matrix (``IMAS`` = :math:`0`) or is specified
   by the function ``Mfct`` (``IMAS`` = :math:`1`). Note that the
   adjoint dynamics requires the transposed mass matrix that has to be
   provided by the function ``Mtrans``.

-  ``MLJAC``: Gives information on the banded structure of the Jacobian
   :math:`\frac{\partial^{} \mb{f}}{\partial \mb{x}^{}}`
   and
   :math:`(\frac{\partial^{} \mb{f}}{\partial \mb{x}^{}})^\mathsf{T}`,
   respectively. The Jacobian is either a full matrix
   (``MLJAC`` = :math:`N_x`) or is of banded structure. In the latter
   case, the number of non-zero diagonals below the main diagonal are
   specified by :math:`0\leq\,`\ ``MLJAC``\ :math:`\,<N_x`.

-  ``MUJAC``: Specifies the number of non-zero diagonals above the main diagonal
   of the Jacobian
   :math:`\frac{\partial^{} \mb{f}}{\partial \mb{x}^{}}`.
   This flag needs not to to be defined if ``MLJAC`` = :math:`N_x`.
   Since the partial derivative of the right hand side of the adjoint
   dynamics with respect to the adjoint state
   :math:`\mb{\lambda}` is given by
   :math:`(\frac{\partial^{} \mb{f}}{\partial \mb{x}^{}})^\mathsf{T}`,
   the meaning of the flags ``MLJAC`` and ``MUJAC`` switches in this
   case.

-  ``MLMAS`` and ``MUMAS`` : Both options have the same meaning as ``MLJAC`` and
   ``MUJAC``, but refer to the mass matrix :math:`\mb{M}`.

If a semi-implicit problem (option ``IMAS`` = :math:`1`) with Mayer term (option
``TerminalCost`` = ``on``) is considered, the terminal conditions of the adjoint system
must be specified in a specific form. To provide
RODAS :footcite:`Hairer:Book:1996:Stiff,Rodas:Webpage:2018` the
proper terminal condition :math:`\mb{\lambda}(T)`, the function ``dVdx``
must be specified as follows

.. math::

   \mb{\lambda}(T) = \underbrace{\left(\mb{M}^\mathsf{T}\right)^{-1}\bar{ V}_{\mb{x}}(\mb{x}(T), \mb{p}, T, \mb{\mu}_T, \mb{c}_T)}_{\texttt{dVdx}}

cf. Equation :math:numref:`eq:AlgOpt:OptCondLambda`.
In the case of a DAE system, the mass matrix is singular and therefore
only the elements of the mass matrix for the differential equations are
inverted. For an example of using RODAS and the respective options, take
a look at the MPC problem ``Reactor_PDE`` in the folder ``<grampc_root>/examples``.

.. footbibliography::