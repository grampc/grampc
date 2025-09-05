Using GRAMPC in Matlab/Simulink
-------------------------------

As already mentioned, the main components of GRAMPC are implemented in plain C to ensure a high level of portability. 
However, GRAMPC also provides a user-friendly interface with Matlab/Simulink to allow for a convenient MPC design.

.. _sec:grampcInterfaceMatlab:

Interface to Matlab
~~~~~~~~~~~~~~~~~~~

Each main component of GRAMPC (cf. :ref:`sec:MainComponents`) has a related Mex routine that is included in the directory ``<grampc_root>/matlab/src``, also see :numref:`fig:grampcInterfaceMatlab`. 
This allows one to run GRAMPC in Matlab/Simulink as well as altering options and parameters without recompilation.

A makefile to compile GRAMPC for use in Matlab/Simulink is provided in ``<grampc_root>/matlab``. 
The makefile compiles the source files to generate object files within ``<grampc_root>/matlab/bin``. 
In order to obtain an actual Mex compilation for a given problem, the object files must be linked against the object file obtained by compiling the problem function, since at least some of these functions depend on the actual problem formulation, e.g. the state dimension :math:`N_{\mb{x}}`.

The m-file ``startMPC.m`` in each of the examples under ``<grampc_root>/examples`` contains a flag ``compile``. 
Setting this flag to 1 leads to a compilation of the problem file and to the generation of the Mex files. 
Setting ``compile`` to 2 leads to an additional recompilation of the toolbox by calling the makefile in the directory ``<grampc_root>/matlab/src``.
The Mex files are stored in the local subfolder ``+CmexFiles``. 
The folder name begins with a plus sign allowing the user to call the functions stored in this
folder using the command ``CmexFiles.<function name>``. 
The S-function files are stored in the local application directory, since Simulink requires the S-functions to be in the Matlab path.


.. figure:: ../img/InterfaceMatlab.*
    :name: fig:grampcInterfaceMatlab

    Interface of GRAMPC to Matlab/Simulink (gray - C code, white - Matlab code)


The structure of the main components of GRAMPC in C and Matlab for setting options and parameters are slightly different, as it is not allowed (or at least not very elegant) to manipulate the input argument of Mex routines. 
Consequently, each Mex routine returns the manipulated structure variable ``grampc`` as an output argument. 
If, for example, the initial condition :math:`\mb{x}_0` should be set to a specific value in C, the function ``grampc_setparam_real_vector`` must be used, as already discussed in :ref:`sec:grampcStructure:initialization`. 
In Matlab, however, this is done using the Mex routine ``grampc_setparam_Cmex`` with the structure variable ``grampc`` as an input argument. 
The manipulated structure variable ``grampc`` is returned as an output argument including the initial condition :math:`x_0`. 
For instance, the parameter setting in C

.. code-block:: c

   ctypeRNum x0[NX] = {-1,-1};
   grampc_setparam_real_vector(grampc, "x0", x0);

would read as follows in Matlab:

::

   grampc = grampc_setparam_Cmex(grampc,'x0',[-1;-1]);

Note that the Mex routine ``grampc_setparam_Cmex`` does not distinguish between vectors and scalars, but handles the different dimensions of parameters internally.
The same applies to the Mex routine ``grampc_setopt_Cmex`` for changing algorithmic options in GRAMPC. 
The data type of the parameter or option to be set can either be double or the corresponding data type in the parameter structure ``param`` or option structure ``opt``, see also Table :numref:`tab:ListOfParameters` or
Table :numref:`tab:ListOfOptions`.

In order to simplify changing parameters and options in Matlab, GRAMPC also provides the routine ``grampc_update_struct_grampc(grampc,user)`` included in ``<grampc_root>/matlab/mfiles``. 
The purpose of this function is to allow the user to define the options and parameters to be set as structure variable instead of requiring to call the functions ``grampc_setparam_Cmex`` and ``grampc_setopt_Cmex`` manually for each chosen parameter/option. 
In detail, the structure variable ``user`` must contain the substructures ``param`` and ``opt`` that define the parameters and options to be set. 
The corresponding function call under Matlab is as follows:

::

   ...
   %% Parameter definition
   % Initial values of the states
   user.param.x0 = [-1;-1];
   ...

   %% Option definition
   % Basic algorithmic options
   user.opt.Nhor = 10;
   ...

   %% User parameter definition
   % e.g. system parameters or weights for the cost function
   userparam = [100, 10, 1, 100, 10, -0.2, 0.01, -0.1, 0.1];

   %% Grampc initialization
   grampc = CmexFiles.grampc_init_Cmex(userparam);

   %% Update grampc struct while ensuring correct data types
   grampc = grampc_update_struct_grampc(grampc,user);
   ...

Similar to :ref:`sec:grampcStructure:initialization` and :ref:`sec:CompileRun`, the following lines show the main steps to run the ball-on-plate example in Matlab. 
An executable version of this example within Matlab can be found in the folder ``<grampc_root>/examples/BallOnPlate``. 
In analogy to the C implementation, the simulation loop and the evaluation are implemented in the main file ``startMPC.m``. 
The parameters and options are defined in the separate file ``initData.m`` that is called within ``startMPC.m`` for the sake of readability and to use the settings directly in the Simulink model,
see :ref:`sec:grampcInterfaceSimulink`. Please note a template file can be found in the folder ``<grampc_root>/examples/TEMPLATES``.

| **Example** **(Matlab code for setting options and parameters, see initData.m)**

::

   %% Parameter definition
   user.param.x0    = [ 0.1, 0.01];
   user.param.xdes  = [-0.2, 0.0];

   % Initial values, setpoints and limits of the inputs
   user.param.u0    = 0;
   user.param.udes  = 0;
   user.param.umax  =  0.0524;
   user.param.umin  = -0.0524;

   % Time variables
   user.param.Thor  = 0.3;         % Prediction horizon

   user.param.dt    = 0.01;        % Sampling time
   user.param.t0    = 0.0;         % time at the current sampling step

   %% Option definition
   user.opt.Nhor        = 10;      % Number of steps for the system integration
   user.opt.MaxMultIter = 3;       % Maximum number of augmented Lagrangian iterations

   % Constraints thresholds
   user.opt.ConstraintsAbsTol = 1e-3*[1 1 1 1];

   %% User parameter definition, e.g. system parameters or weights for the cost function
   userparam = [100, 10, 180, 100, 10, -0.2, 0.2, -0.1, 0.1];

   %% Grampc initialization
   grampc = CmexFiles.grampc_init_Cmex(userparam);

   %% Update grampc struct while ensuring correct data types
   grampc = grampc_update_struct_grampc(grampc,user);

   %% Estimate and set PenaltyMin (optional)
   grampc = CmexFiles.grampc_estim_penmin_Cmex(grampc,1);ot_stat(vec,grampc,phpS);
   ...


| **Example** **(Matlab code for running GRAMPC within an MPC loop, see startMPC.m)**

::

   ...
   %% Initialization
   [grampc,Tsim] = initData;
   CmexFiles.grampc_printopt_Cmex(grampc);
   CmexFiles.grampc_printparam_Cmex(grampc);

   % init solution structure
   vec = grampc_init_struct_sol(grampc, Tsim);

   % init plots and store figure handles
   phpP = grampc_init_plot_pred(grampc,figNr);     figNr = figNr+1;
   phpT = grampc_init_plot_sim(vec,figNr);         figNr = figNr+1;
   phpS = grampc_init_plot_stat(vec,grampc,figNr); figNr = figNr+1;

   %% MPC loop
   i = 1;
   while 1
     % set current time and current state
     grampc = CmexFiles.grampc_setparam_Cmex(grampc,'t0',vec.t(i));
     grampc = CmexFiles.grampc_setparam_Cmex(grampc,'x0',vec.x(:,i));

     % run MPC and save results
     [grampc,vec.CPUtime(i)] = CmexFiles.grampc_run_Cmex(grampc);
     vec = grampc_update_struct_sol(grampc, vec, i);

     % print solver status
     printed = CmexFiles.grampc_printstatus_Cmex(grampc.sol.status,'Error');

     % check for end of simulation
     if i+1 > length(vec.t)
       break;
     end

     % simulate system
     [~,xtemp] = ode45(@CmexFiles.grampc_ffct_Cmex,vec.t(i)+[0 double(grampc.param.dt)],
                 vec.x(:,i),odeopt,vec.u(:,i),vec.p(:,i),grampc.param,grampc.userparam);
     vec.x(:,i+1) = xtemp(end,:);

     % evaluate time-dependent constraints
     vec.constr(:,i) = CmexFiles.grampc_ghfct_Cmex(vec.t(i), vec.x(:,i), vec.u(:,i),
                                            vec.p(:,i), grampc.param, grampc.userparam);

     % update iteration counter
     i = i + 1;

     % plot data
     grampc_update_plot_pred(grampc,phpP);
     grampc_update_plot_sim(vec,phpT);
     grampc_update_plot_stat(vec,grampc,phpS);
   end

Similar to the C example in :ref:`sec:CompileRun`, the structure variable ``grampc`` is initialized before the options as well as optional parameters are set. 
In addition, the plot functions (see :ref:`sec:Plotfunctions`) are initialized before GRAMPC is started within an MPC loop, where the current state of the system (new initial condition) is provided to GRAMPC. 
After computing the new controls, the status of GRAMPC is printed, see :ref:`sec:AlgOpt:StatusFlags` for more details.
Subsequently, a reference integration of the system is performed, and the constraints are evaluated before the plots are updated.

.. _sec:grampcInterfaceSimulink:

Interface to Matlab/Simulink
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GRAMPC also allows a Matlab/Simulink integration via the S-function ``grampc_run_Sfct.c`` (also included in the directory ``<grampc_root>/matlab/src``). 
A corresponding Simulink template can be found in the folder ``<grampc_root>/examples/TEMPLATES`` for a number of Matlab versions.
The directory also contains the m-file ``initData_TEMPLATE.m`` which can be used for initializing GRAMPC's options and parameters, as mentioned in the previous subsection. 
The build procedure of the Mex routines additionally compiles the S-function for the Simulink block.

The Matlab/Simulink model of GRAMPC is shown in :numref:`fig:grampcSimulink`. 
The block ``MPC-Subsystem`` contains algorithmic components of GRAMPC (implemented within the S-function ``grampc_run_Sfct.c``). 
The block ``Click to init grampc`` must be executed by a double click in order to initialize the structure variable that is required by the Matlab/Simulink model. 
This generates also the Simulink-specific structure variable ``grampc_sdata``. 
For the sake of convenience, the blocks ``Click to compile toolbox`` and ``Click to compile probfct`` are included in the model to be able to compile the GRAMPC toolbox and the specific problem directly from Matlab/Simulink.

The block ``System function`` is a Matlab Function block which implements the system dynamics in order to numerically integrate the system dynamics after each MPC step :math:`k` for the sampling time :math:`\Delta t` and to return the new state value :math:`\mb{x}_{k+1}` corresponding to the next sampling instant :math:`t_{k+1}` that is fed back to the MPC block.

Furthermore, the S-function ``grampc_run_Sfct.c`` satisfies the additional restrictions of the Matlab code generator. 
Therefore, the block can be used in models implemented for running on various hardware platforms, such as dSpace real-time systems. 
Please note that especially in case of dSpace applications, the include folders ``<grampc_root>/include`` and ``<grampc_root>/matlab/include`` as well as all C source files in ``<grampc_root>/src``,
the source file of the S-function ``<grampc_root>/matlab/src/grampc_run_Sfct.c`` and the problem function must be listed as additional build information in the ``Model Configuration Parameters`` of the Simulink model under ``Code Generation / Custom Code``. 
It is recommended to use absolute paths at least for the S-function file.

.. figure:: ../img/grampc_Simulink.*
    :name: fig:grampcSimulink

    Matlab/Simulink model of GRAMPC.

.. _sec:Plotfunctions:

Plot functions
~~~~~~~~~~~~~~

GRAMPC offers various plot functions in the folder ``<grampc_root>/matlab/mfiles``. Each plot must be initialized at first using the m-files ``grampc_init_plot_*.m``. 
During the simulation the plots can be updated by the m-files ``grampc_update_plot_*.m``. 
Beside the trajectories of the simulated system dynamics (file ending \*=sim) and trajectories along the prediction horizon (file ending \*=pred), also some statistics (file ending \*=stat) can be plotted. 
When solving OCPs instead of MPC problems, the plot along the prediction horizon shows the actual results. 
The plotted quantities depend on the parameter and option settings of GRAMPC, i.e. whether constraints are considered or not. 
The available plots are listed in more detail in the following lines. 
(Also see the example problems under ``<grampc_root>/examples`` for code samples on how to use the plot routines.)

System dynamics plot (plot_sim)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  **States**: This plot illustrates the trajectories of the state :math:`\mb{x}` along the simulation time.

-  **Adjoint states**: This plot illustrates the trajectories of the adjoint state :math:`\mb{\lambda}` along the simulation time.

-  **Controls**: This plot illustrates the trajectories of the controls :math:`\mb{u}` along the simulation time.

-  **Constraints**: This plot appears only if equality and/or inequality constraints are defined (:math:`{N_{\mb{g}}}` and/or :math:`{N_{\mb{h}}}` is larger than zero as specified in ``ocp_dim``). 
   The plot shows the trajectories of the constraints :math:`\mb{g}` and :math:`\mb{h}` along the simulation time.

-  **Lagrange multipliers**: This plot appears only if equality and/or inequality constraints are defined. 
   The plot shows the trajectories of the multipliers :math:`\mb{\mu}_{\mb{g}}` and :math:`\mb{\mu}_{\mb{h}}` along the simulation time.
   If any Lagrange multiplier reaches the limit :math:`\mu_\text{max}` (specified by ``MultiplierMax``), it indicates that the penalty parameters are too high or that the problem is not well-conditioned or that the costs are badly scaled.

-  **Penalty parameters**: This plot appears only if equality and/or inequality constraints are defined. 
   The plot shows the trajectories of the penalties :math:`\mb{c}_{\mb{g}}` and :math:`\mb{c}_{\mb{h}}` along the simulation time. 
   If any penalty reaches the maximum value :math:`c_\text{max}`, set by ``PenaltyMax``, it indicates that either the limit is not high enough or the update is too aggressive.

Prediction plot (plot_pred)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  **Predicted states**: This plot illustrates the trajectories of the state :math:`\mb{x}` along the prediction horizon.

-  **Predicted adjoint states**: This plot illustrates the trajectories of the adjoint state :math:`\mb{\lambda}` along the prediction horizon.

-  **Predicted controls**: This plot illustrates the trajectories of the controls :math:`\mb{u}` along the prediction horizon.

-  **Predicted constraints**: This plot appears only if (terminal) equality and/or inequality constraints are defined. 
   The plot shows the predicted violation of the equality constraints :math:`\mb{g}` and inequality constraints :math:`\mb{\max}(\mb{h},\mb{0})` 
   along the prediction horizon and the predicted violation of the terminal equality constraints :math:`{\mb{g}_{T}}` and inequality
   constraints :math:`\mb{\max}(\mb{h}_T,\mb{0})` at the end of the prediction horizon. 
   Please note that except for OCPs, these are not the actual but predicted internal constraint violations of the current GRAMPC iteration.

-  **Lagrange multipliers**: This plot appears only if (terminal)
   equality and/or inequality constraints are defined. 
   The plot shows the trajectories of the multipliers :math:`\mb{\mu}_{\mb{g}}` and
   :math:`\mb{\mu}_{\mb{h}}` along the prediction horizon and the multipliers
   :math:`{\mb{\mu}_{\mb{g}_T}}` and :math:`{\mb{\mu}_{\mb{h}_T}}`. 
   If any Lagrange multiplier reaches the limit :math:`\mu_\text{max}`, set by ``MultiplierMax``, it
   indicates that the penalty parameters are too high or that the problem is not well-conditioned or that the costs are badly scaled.

-  **Penalty parameters**: This plot appears only if (terminal) equality and/or inequality constraints are defined. 
   The plot shows the trajectories of the penalties :math:`\mb{c}_{\mb{g}}` and :math:`\mb{c}_{\mb{h}}` 
   along the prediction horizon and the penalties :math:`{\mb{c}_{\mb{g}_T}}` and :math:`{\mb{c}_{\mb{h}_T}}`. 
   If any penalty reaches the maximum value :math:`c_\text{max}` set by ``PenaltyMax``, it indicates that either the limit is not high enough or 
   the update is too aggressive, see also :ref:`sec:AlgOpt:UpdateMultPen`.

Statistics plot (plot_stat)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  **Costs**: This plot illustrates the costs :math:`J` along the simulation time or along the augmented Lagrangian iterations. 
   If constraints are defined, the augmented costs :math:`\bar{J}` are shown as well.

-  **Computation time**: This plot illustrates the computation time of one MPC or optimization step of GRAMPC along the simulation time or along the augmented Lagrangian iterations. 
   The time measurement is done in the function ``grampc_run_Cmex.c`` using operating system specific timer functions.
   Consequently, the time excludes the overhead resulting from the Mex interface as well as the time consumed by the plot functions.

-  **Line search step size**: This plot illustrates the step size :math:`\alpha` of the last gradient iteration along the simulation time or along the augmented Lagrangian iterations. 
   If the adaptive line search is used (see :ref:`sec:AlgOpt:LineSearchAdaptive`), the plot also illustrates the three corresponding sample points :math:`\alpha_1`, :math:`\alpha_2`, and :math:`\alpha_3`. 
   A step size equal to the maximum or minimum value :math:`\alpha_\text{max}` or :math:`\alpha_\text{min}` indicates that these values may have to be adapted or the problem may have to be scaled. 
   Additionally, if the explicit line search is chosen and the fallback strategy is not activated (see :ref:`sec:AlgOpt:LineSearchExplicit` and :ref:`sec:AlgOpt:LineSearchFallback`), a frequent use of the initial value :math:`\alpha_\text{init}` indicates an ill-conditioned problem.

-  **Gradient iterations**: This plot appears only if the option ``ConvergenceCheck`` is set to ``on``. 
   It illustrates the number of executed gradient iterations along the simulation time or along the augmented Lagrangian iterations. 
   In particular, the plot depicts whether the maximum number of gradient iterations :math:`j_\text{max}` is reached or the convergence check caused a premature termination of the minimization.

-  **Prediction horizon**: This plot appears only if the option ``OptimTime`` is set to ``on``. 
   It illustrates the prediction horizon :math:`T` along the simulation time or along the augmented Lagrangian iterations. 
   In shrinking horizon applications the value should decrease linearly after a short settling phase.

-  **Norm of constraints over horizon**: This plot appears only if constraints are defined. 
   It illustrates the norm :math:`\frac{1}{T}\sqrt{\|\mb{g}\|_{L_2}^2+\|\mb{\max}( {\mb{h}},\mb{0})\|_{L_2}^2 +\| {\mb{g}_T}\|_{2}^2+\|\mb{\max}( {\mb{h}_T},\mb{0})\|_{2}^2 }`
   over all predicted constraints plotted over the simulation time or the augmented Lagrangian iterations. 
   Especially when solving OCPs, the value should decrease continuously.

-  **Norm of penalty parameters over horizon**: This plot appears only if the number of equality :math:`{N_{\mb{g}}}`, inequality :math:`{N_{\mb{h}}}`, terminal equality :math:`{N_{\mb{g}_T}}` or terminal inequality :math:`{N_{\mb{h}_T}}` constraints is not zero. 
   It illustrates the norm :math:`\frac{1}{T}\sqrt{\|\mb{\bar c}\|_{L_2}}` over all predicted penalty parameters along the simulation time or along the augmented Lagrangian iterations.

-  **Terminal constraints**: This plot appears only if terminal constraints are defined. 
   It illustrates the violation of the terminal equality constraints :math:`{\mb{g}_{T}}` and inequality constraints :math:`\mb{\max}(\mb{h}_T,\mb{0})` along the simulation time or along the augmented Lagrangian iterations in case of OCPs.

-  **Terminal Lagrangian multipliers**: This plot appears only if terminal constraints are defined. 
   It illustrates the multipliers :math:`{\mb{\mu}_{\mb{g}_T}}` and :math:`{\mb{\mu}_{\mb{h}_T}}` along the simulation time or along the augmented Lagrangian iterations in case of OCPs.

-  **Terminal penalty parameters**: This plot appears only if terminal constraints are defined. 
   It illustrates the penalties :math:`{\mb{c}_{\mb{g}_T}}` and :math:`{\mb{c}_{\mb{h}_T}}` along the simulation time or along the augmented Lagrangian iterations in case of OCPs.

.. footbibliography::
