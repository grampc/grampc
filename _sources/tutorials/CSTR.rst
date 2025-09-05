.. _sec:Tut:CSTR:

Moving horizon estimation of a CSTR
-----------------------------------

Continuous stirred tank reactors (CSTR) are a popular class of systems
when it comes to the implementation of advanced nonlinear control
methods. In this subsection, a CSTR model is used for an example
implementation of a moving horizon estimation (MHE) used in closed loop
with MPC. Corresponding m and C files can be found in the folder ``<grampc_root>/examples/Reactor_CSTR``.

.. _problem-formulation-2:

Problem formulation
~~~~~~~~~~~~~~~~~~~

The system dynamic are given by :footcite:`Rothfuss1996`

.. math::

    \dot c_\mathrm{A} &= -k_1(T) c_\mathrm{A} - k_2(T) c_\mathrm{A}^2 + (c_{\mathrm{in}} - c_\mathrm{A}) u_1 \\
    \dot c_\mathrm{B} &= k_1(T) c_\mathrm{A} - k_1(T) c_\mathrm{B} - c_\mathrm{B} u_1 \\
    \dot{T}~ &= -\delta ( k_1(T) c_\mathrm{A} \Delta H_\mathrm{AB} + k_1(T) \Delta H_\mathrm{BC} + k_2(T) c_\mathrm{A}^2 \Delta H_\mathrm{AD}) \nonumber\\
    &\hspace{1cm}+ \alpha (T_\mathrm{C} - T) + (T_\mathrm{in} - T) u_1 \\
    \dot{T}_\mathrm{C} &= \beta (T - T_\mathrm{C}) + \gamma u_2\,,

where the monomer and product concentrations :math:`c_\mathrm{A}` and
:math:`c_\mathrm{B}`, respectively, as well as the reactor and cooling
temperature :math:`T` and :math:`T_\mathrm{C}` form the state vector
:math:`\mb{x} =[c_\mathrm{A}, c_\mathrm{B}, T, T_\mathrm{C}]`. The two functions
:math:`k_1(T)` and :math:`k_2(T)` are of Arrhenius type

.. math::

    k_i(T) = k_{i0}  \exp \left(\frac{-E_i}{T/\mathrm{^\circ C} + 273.15} \right) ,\quad i =  1,2\,.

The measured quantities are the two temperatures
:math:`\mb{y}=[T,T_C]^\mathsf{T}`. The controls
:math:`\mb{u} = [u_1, u_2]`, i.e. the normalized flow rate and
the cooling power, are assumed to be known as well.
:numref:`tab:CSTRParams` gives the parameters of the system that
are passed to the GRAMPC problem functions via ``userparam``. A more detailed
description can be found in :footcite:`Rothfuss1996`.
  
.. list-table:: Parameters of the CSTR model :footcite:`Rothfuss1996`.
    :name: tab:CSTRParams
    :widths: auto
    :header-rows: 1

    * - Parameter
      - Value
      - Unit
      - Parameter
      - Value
      - Unit
    * - :math:`\alpha`
      - 30.828
      - h :sup:`-1`
      - :math:`k_{20}`
      - 9.043 x 10 :sup:`6`
      - m :sup:`3` mol :sup:`-1` h :sup:`-1`
    * - :math:`\beta`
      - 86.688
      - h :sup:`-1`
      - :math:`E_{1}`
      - 9785.3
      - 
    * - :math:`\delta`
      - 3.522 x 10 :sup:`-4`
      - K kJ :sup:`-1`
      - :math:`E_{2}`
      - 8560.0
      - 
    * - :math:`\gamma`
      - 0.1
      - kh :sup:`-1`
      - :math:`\Delta H_{AB}`
      - 4.2
      - kJ mol :sup:`-1`
    * - :math:`T_{in}`
      - 104.9
      - °C
      - :math:`\Delta H_{BC}`
      - -11.0
      - kJ mol :sup:`-1`
    * - :math:`c_{in}`
      - 5.1 x 10 :sup:`3`
      - mol m :sup:`-3`
      - :math:`\Delta H_{AD}`
      - -41.85
      - kJ mol :sup:`-1`
    * - :math:`k_{10}`
      - 1.287 x 10 :sup:`12`
      - h :sup:`-1`
      - 
      - 
      - 

The control task at hand is the setpoint change between the two
stationary points

.. math::

   \mb{x}_\mathrm{des,1} = [1370\,\mathrm{\frac{kmol}{m^3}}, 950\,\mathrm{\frac{kmol}{m^3}}, 110.0\,\mathrm{^\circ C}, 108.6\,\mathrm{^\circ C}]^\mathsf{T},  \mb{u}_\mathrm{des,1} = [5.0\,\mathrm{h^{-1}}, -1190\,\mathrm{kJ\,h^{-1}}]^\mathsf{T}

and

.. math::
   \mb{x}_\mathrm{des,2} = [2020\,\mathrm{\frac{kmol}{m^3}}, 1070\,\mathrm{\frac{kmol}{m^3}}, 100.0\,\mathrm{^\circ C}, 97.1\,\mathrm{^\circ C}]^\mathsf{T},  \mb{u}_\mathrm{des,2} = [5.0\,\mathrm{h^{-1}}, -2540\,\mathrm{kJ\,h^{-1}}]^\mathsf{T}.

The cost functional is designed quadratically

.. math:: J(\mb{u}, \mb{x}_k) := \Delta \mb{x}(T)^\mathsf{T}\mb{P} \Delta \mb{x}(T) + \int_{0}^{T} \Delta \mb{x}(t)^\mathsf{T}\mb{Q} \Delta \mb{x}(t) \Delta \mb{u}(t)^\mathsf{T}\mb{R} \Delta \mb{u}(t)

in order to penalize the deviation of the state and control from the
desired setpoint
:math:`(\mb{x}_\mathrm{des,1},\mb{u}_\mathrm{des,1})`
with
:math:`\Delta \mb{x} = \mb{x} - \mb{x}_\mathrm{des}`
and
:math:`\Delta \mb{u} = \mb{u} - \mb{u}_\mathrm{des}`.
The weight matrices are chosen as

.. math::

    \mb{P} = \text{diag}(0.2, 1.0, 0.5,0.2),\,\, \mb{R} = \text{diag}(0.5, 5.0\mathrm{e-3}),\,\, \mb{Q} = \text{diag}(0.2, 1.0, 0.5, 0.2).

The control task will be tackled by MPC using GRAMPC. In addition, an
MHE using GRAMPC is designed to estimate the current state
:math:`\mb{\hat x}_k` w.r.t. the measured temperatures
:math:`\mb{y}=[T,T_C]^\mathsf{T}`.

In analogy to the MPC formulation :math:numref:`eq:OCP`,
moving horizon estimation is typically based on the online solution of a
dynamic optimization problem

.. math::
    :label: eq:MHE_1_orig

   	\min_{\mb{\hat x}_k} \quad & J(\mb{\hat x}_k; \mb{u}, \mb{y}) = \int_{t_k-T}^{t_k} \| \mb{\hat y}(t) - \mb{y}(t) \|^2 \, {\rm d}t

   	\,\textrm{s.t.} \quad & \mb{M}\mb{\dot{\hat x}}(t) = \mb{f}(\mb{\hat x}(t), \mb{u}(t),t) \,,\quad \mb{\hat x}(t_k) = \mb{\hat x}_k

    & \mb{\hat y}(t) = \mb{\sigma}(\mb{\hat x}(t))

that depends on the history of the measured outputs
:math:`\mb{y}(t)` and controls :math:`\mb{u}(t)` in the
past time window :math:`[t_k - T, t_k]`. The solution of :math:numref:`eq:MHE_1_orig` yields the estimate of the state
:math:`\mb{\hat{x}}_k` such that the difference between the
measured output :math:`\mb{y}(t)` and the estimated output
function :math:`\mb{\hat y}(t) = \mb{\sigma}(\mb{\hat x}(t))` is minimal in the sense of
:math:numref:`eq:MHE_1_orig`. GRAMPC solves this MHE problem by
means of parameter optimization. To this end, the state at the beginning
of the optimization horizon is defined as optimization variable,
i.e. :math:`\mb{p} = \mb{\hat{x}}(t_k - T)`.

Both MHE and MPC use a sampling rate of
:math:`\Delta t = 1\,\mathrm{s}`. A prediction horizon of
:math:`T = 20\,\mathrm{min}` with 40 discretization points is used for the
MPC, while a prediction horizon of :math:`T = 10\,\mathrm{s}` with 10
discretization points is used for the MHE. The MPC implementation uses
three gradient iterations per sampling step,
i.e. :math:`(i_\text{max},j_\text{max})=(1,3)`, while the implementation
of the MHE uses only a single gradient iteration,
i.e. :math:`(i_\text{max},j_\text{max}) =(1,1)`. Note that because the
MHE and MPC problems are defined without state constraints, the outer
augmented Lagrangian loop causes no computational overhead, as GRAMPC
skips the multiplier and penalty update. As the implementation of the
MHE is not quite as straightforward as the MPC case, the next subsection
describes the implementation process in more detail.

Implementation aspects
~~~~~~~~~~~~~~~~~~~~~~

The following lines describe the implementation of the MHE problem with
GRAMPC, the corresponding simulation results are shown in the next
subsection. In a first step, the MHE
problem :math:numref:`eq:MHE_1_orig` has to be transformed in a
more suitable representation that can be tackled with the parameter
optimization functionality of GRAMPC. To this end, a coordinate
transformation

.. math::
    :label: eq:coordTrafo

    \mb{\tilde x}(\tau) = \mb{\hat x}(t_k\!-\!T\!+\!\tau) - \mb{p} \,,\quad
    \mb{\tilde u}(\tau) = \mb{u}(t_k\!-\!T\!+\!\tau) \,,\quad
    \mb{\tilde y}(\tau) = \mb{y}(t_k\!-\!T\!+\!\tau)

is used together with the corresponding time transformation from
:math:`t\in [t_k-T,t_k]^\mathsf{T}` to the new time coordinate
:math:`\tau \in [0, T]`. In combination with the optimization variable
:math:`\mb{p} = \mb{\hat{x}}(t_k - T)` and the
homogeneous initial state
:math:`\mb{\tilde{x}}(0) = \mb{0}`, the optimization
problem can be cast into the form

.. math::
    :label: eq:MHE_orig

   	\min_{\mb{p}} \quad & J(\mb{p}; \mb{\tilde u}, \mb{\tilde y}) = \int_0^T \| \mb{\hat y}(\tau) - \mb{\tilde y}(\tau) \|^2 \, {\rm d}\tau

   	\textrm{s.t.} \quad & \mb{\dot{\tilde x}}(\tau) = \mb{f}(\mb{\tilde x}(\tau) + \mb{p}, \mb{\tilde u}(\tau), t_k\!-\!T\!+\!\tau) \,,\quad \mb{\tilde x}(0) = \mb{0}

    & \mb{\hat y}(\tau) = \mb{\sigma}(\mb{\tilde x}(\tau)+\mb{p})\,.

The implementation of this optimization problem still requires to access
the measurements :math:`\mb{\tilde y}` in the integral cost
term. This is achieved by appending the measurements to ``userparam`` (see ``startMHE.m`` in ``<grampc_root>/examples/Reactor_CSTR``) 

::

   	% init array of last MHE-Nhor measurements of the two temperatures
   	xMeas_array = repmat(grampcMPC.param.x0(3:4), 1, grampcMHE.opt.Nhor);
   	grampcMHE.userparam(end-2*grampcMHE.opt.Nhor+1:end) = xMeas_array;

The measurements are updated in each iteration of the MPC/MHE loop, e.g.

::

   	% set values of last MHE-Nhor measurements of the two temperatures
   	xMeas_temp = xtemp(end,3:4) + randn(1,2)*4; % measurement noise
   	xMeas_array = [xMeas_array(3:end), xMeas_temp];
   	grampcMHE.userparam(end-2*grampcMHE.opt.Nhor+1:end) = xMeas_array;

When the number of discretization points and the horizon length is
known, the measurements can easily be accessed in the problem
description file in the following way:

::

   	typeRNum* pSys = (typeRNum*)userparam;
   	typeRNum* pCost = &pSys[14];
   	typeRNum* pMeas = &pSys[20];
   	typeInt index = (int)floor(t / 2.777777777777778e-04 + 0.00001);
   	typeRNum meas1 = pMeas[2 * index];
   	typeRNum meas2 = pMeas[1 + 2 * index];

The pointer ``pMeas`` is set to the 20th element of the ``userparam`` vector, since the first 14
are the system parameters given in :numref:`tab:CSTRParams` and
the next six elements are the weights of the cost function. Also note
that :math:`2.778e-4\,\mathrm{h^{-1}} \approx \Delta t`. Since the order of
magnitude of the individual states and controls differs a lot, the
scaling option ``ScaleProblem`` with

.. math::

   \mb{x}_{\mathrm{scale}}  &= [500\,\mathrm{\frac{kmol}{m^3}}, 500\,\mathrm{\frac{kmol}{m^3}}, 50\,\mathrm{^\circ C}, 50\,\mathrm{^\circ C}]^\mathsf{T}& \mb{x}_{\mathrm{offset}} &= [500\,\mathrm{\frac{kmol}{m^3}}, 500\,\mathrm{\frac{kmol}{m^3}}, 50\,\mathrm{^\circ C}, 50\,\mathrm{^\circ C}]^\mathsf{T}

   \mb{p}_{\mathrm{scale}}  &= [500\,\mathrm{\frac{kmol}{m^3}}, 500\,\mathrm{\frac{kmol}{m^3}}, 50\,\mathrm{^\circ C}, 50\,\mathrm{^\circ C}]^\mathsf{T}& \mb{p}_{\mathrm{offset}} &= [500\,\mathrm{\frac{kmol}{m^3}}, 500\,\mathrm{\frac{kmol}{m^3}}, 50\,\mathrm{^\circ C}, 50\,\mathrm{^\circ C}]^\mathsf{T}

   \mb{u}_{\mathrm{scale}}  &= [16\,\mathrm{h^{-1}}, 4500\,\mathrm{kJ h^{-1}}]^\mathsf{T}&	\mb{u}_{\mathrm{offset}} &= [19\,\mathrm{h^{-1}}, -4500\,\mathrm{kJ h^{-1}}]^\mathsf{T}

is activated.

Evaluation
~~~~~~~~~~

.. figure:: ../img/tikz/CSTR_closedLoop.*
    :name: fig:CSTR_closedLoop

    Simulated MHE/MPC trajectories for the CSTR reactor example.

The moving horizon estimator is evaluated in conjunction with the MPC.
The state estimates are initialized with an initial disturbance
:math:`\delta \mb{p} = [100\,\mathrm{\frac{kmol}{m^3}}`,
:math:`100\,\mathrm{\frac{kmol}{m^3}}`,
:math:`5\,\mathrm{^\circ C}`,
:math:`7\,\mathrm{^\circ C}]^\mathsf{T}`. For a more realistic
setting, white Gaussian noise with zero mean and a standard deviation of
:math:`4\,\mathrm{^\circ C}` is added to the measurements
:math:`\mb{y}=[T,T_C]^\mathsf{T}`.

:numref:`fig:CSTR_closedLoop` shows the simulation results from
the closed loop simulation of the MHE in conjunction with the MPC. The
estimates :math:`\mb{\hat x}_k` quickly converge to the actual
states (ground truth), as e.g. can be seen in the upper right plot.
Furthermore, the cost of the MPC quickly converges to zero after each
setpoint change at :math:`t=0\,\mathrm{h}` and
:math:`t=1.5\,\mathrm{h}`, respectively, which shows the good
performance of the combined MPC/MHE problem.

The corresponding computation times of GRAMPC amount to
:math:`58\,\mathrm{\mu s}` and :math:`11\,\mathrm{\mu s}` per
MPC and MHE step, respectively, on a Windows 10 machine with Intel(R)
Core(TM) i5-5300U CPU running at 2.3GHz using the Microsoft Visual C++
2013 Professional (C) compiler.

.. footbibliography::