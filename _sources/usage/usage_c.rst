.. _chap:UsageC:

Using GRAMPC in C
-----------------

A high level of portability is provided by GRAMPC in view of different operating systems and hardware devices. 
The main components are written in plain C without the use of external libraries, also see the discussion in :ref:`sec:Structure`. 
This section describes the usage of GRAMPC in C including initialization, compilation and running the MPC framework.

.. _sec:MainComponents:

Main components of GRAMPC
~~~~~~~~~~~~~~~~~~~~~~~~~

As illustrated in Figure :numref:`fig:grampcMainComponents`, GRAMPC
contains initializing as well as running files, different integrators
for the system and adjoint dynamics (:ref:`chap:ProblemFormulation` and :ref:`chap:AlgOpt`) and
functions to alter the options and parameters (also see :ref:`chap:ProblemFormulation` and :ref:`chap:AlgOpt`).

In more detail, GRAMPC comprises the following main files (see also Appendix :ref:`appendix:FunctionInterface` for the function interface):

- ``grampc_init.c``: Functions for initializing GRAMPC.

- ``grampc_alloc.c``: Functions for memory allocation and deallocation.

- ``grampc_fixedsize.c``: Replaces the functions in ``grampc_alloc.c`` in the fixed-size mode without dynamic memory allocation.

- ``grampc_run.c``: Functions for running GRAMPC including the implemented augmented Lagrangian algorithm and the underlying gradient algorithm. Further functions of this file are concerned with the line search strategies and the update steps of the primal and dual variables as described in :ref:`chap:AlgOpt`.

- ``grampc_setparam.c``: Provides several functions for setting problem-related parameters, see the description in :ref:`chap:ProblemFormulation` and the list of parameters in Table :numref:`tab:ListOfParameters`.

- ``grampc_setopt.c``: Provides several functions for setting algorithmic options, see the description in :ref:`chap:AlgOpt` and the list of options in Table :numref:`tab:ListOfOptions`.

- ``grampc_configreader.c``: Provides an alternative method for reading parameters and options from a text file.

- ``grampc_mess.c``: Function for printing information regarding the initialization as well as execution of GRAMPC (e.g. errors, convergence behavior of the augmented Lagrangian algorithm or status of integrators).

- ``grampc_util.c``: Auxiliary functions for the GRAMPC toolbox such as memory management and trajectory interpolation. It also contains the function ``grampc_estim_penmin`` to compute an estimate of the minimal penalty parameter, cf. :ref:`sec:AlgOpt:EstimPenMin`.

- ``simpson.c``: Function for integrating the integral cost by means of the Simpson rule.

- ``trapezoidal.c``: Function for integrating the integral cost by means of the trapezoidal rule.

- ``grampc_erk.c``: Explicit Runge-Kutta methods of order 1, 2, 3, 4 with fixed step size.

- ``rodas.c``: Semi-implicit Rosenbrock integration scheme with variable step size. The code follows from an f2c conversion of the original Fortran files of RODAS :footcite:`Rodas:Webpage:2018` and are directly included in the GRAMPC source files.

- ``ruku45.c``: Runge-Kutta integration scheme of order 4 with variable step size.

- ``discrete.c``: Discrete integrator.

- ``timing.c``: Optional functions for measuring the execution time.

- ``finite_diff.c``: Finite difference utility functions and gradient checker.

.. figure:: ../img/MainComponents.*
    :name: fig:grampcMainComponents

    Main components of GRAMPC

.. _`sec:grampcStructure:initialization`:

Initialization of GRAMPC
~~~~~~~~~~~~~~~~~~~~~~~~

The global initialization of GRAMPC is done via the routine

::

   void grampc_init(typeGRAMPC **grampc, typeUSERPARAM *userparam)

where the overall structure variable contains the following
substructures (see :ref:`appendix:DataTypes` for the
definition of the data type):

-  ``sol`` (data type ``typeGRAMPCsol``): Contains the optimization variables
   :math:`(\mb{u}^{i+1},\mb{p}^{i+1},T^{i+1})` and the
   interpolated state :math:`\mb{x}^{i+1}` in the next MPC step
   as well as the corresponding cost values
   :math:`\bar J(\mb{u}^{i+1}, \mb{p}^{i+1}, T^{i+1}, \mb{\mu}^{i+1}, \mb{c}^{i+1};\mb{x}_0)` and
   :math:`J(\mb{u}^{i+1}, \mb{p}^{i+1}, T^{i+1};\mb{x}_0)`,
   respectively. The control :math:`\mb{u}^{i+1}` and the state
   :math:`\mb{x}^{i+1}` refer to the corresponding trajectories
   evaluated at the sampling time :math:`\Delta t`. In addition, the
   substructure ``sol`` contains the evaluated functions of the defined
   state constraints, some status information, and an array in which the
   number of gradient iterations are stored in each augmented Lagrangian
   step.

-  ``param`` (data type ``typeGRAMPCparam``): Contains the parameter structure of GRAMPC. 
   A detailed description of all parameters is given in :ref:`chap:ProblemFormulation`.

-  ``opt`` (data type ``typeGRAMPCopt``): Contains the option structure of GRAMPC. 
   A detailed description of all options is given in :ref:`chap:AlgOpt`.

-  ``rws`` (data type ``typeGRAMPCrws``): Contains the real-time workspace of GRAMPC including calculations of the augmented Lagrangian algorithm and the gradient algorithm along the prediction horizon.

-  ``userparam`` (data type ``typeUSERPARAM``): Can be used to define parameters, e.g. to parametrize the cost functional, the system dynamics or the state constraints.

The definition of these data types is given in Appendix :ref:`appendix:DataTypes`. 
The deallocation of is done by means of the function

.. code-block:: c

   void grampc_free(typeGRAMPC **grampc)

.. _sec:setting_opt_par:

Setting options and parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The single options in Table :ref:`tab:ListOfOptions` are set via the functions

.. code-block:: c

   void grampc_setopt_real(const typeGRAMPC *grampc, const typeChar *optName, ctypeRNum optValue)
   void grampc_setopt_real_vector(const typeGRAMPC *grampc, const typeChar *optName, ctypeRNum *optValue)
                           
   void grampc_setopt_int(const typeGRAMPC *grampc, const typeChar *optName, ctypeInt optValue)
   void grampc_setopt_int_vector(const typeGRAMPC *grampc, const typeChar *optName, ctypeInt *optValue)
                          
   void grampc_setopt_string(const typeGRAMPC *grampc, const typeChar *optName, const typeChar *optValue)

for option values with floating point, integer and string type, respectively. 
An overview of the current options can be displayed by using

.. code-block:: c

   void grampc_printopt(typeGRAMPC *grampc)

| **Example** **(Setting options in C)**

The number of gradient iterations, the integration scheme, and the relative tolerance of the integrator can be set in the following way:

.. code-block:: c

   ...
   /********* declaration *********/
   typeGRAMPC *grampc; 
   ...

   /********* option definition *********/
   /* Basic algorithmic options */
   ctypeInt MaxGradIter = 2;

   /* System integration */
   const char* Integrator = "ruku45";
   ctypeRNum IntegratorRelTol = 1e-3;
   ...
   	
   /********* grampc init *********/
   grampc_init(&grampc, userparam);
   ...

   /********* setting options *********/
   grampc_setopt_int(grampc, "MaxGradIter", MaxGradIter);

   grampc_setopt_string(grampc, "Integrator", Integrator);
   grampc_setopt_real(grampc, "IntegratorRelTol", IntegratorRelTol);
   ...

Similar to setting the GRAMPC options, the parameters in Table :ref:`tab:ListOfParameters` are set according to their data type with the following functions:

.. code-block:: c

   void grampc_setparam_real(const typeGRAMPC *grampc, const typeChar *paramName, ctypeRNum paramValue);
                             
   void grampc_setparam_real_vector(const typeGRAMPC *grampc, const typeChar *paramName, ctypeRNum *paramValue);

An overview of the parameters can be displayed by using

.. code-block:: c

   void grampc_printparam(const typeGRAMPC *grampc);

| **Example** **(Setting parameters in C)**

The sampling time, the prediction horizon, and the initial conditions for a system with two states and one control input can be set in the following way:

.. code-block:: c

   ...
   /********* parameter definition *********/
   ctypeRNum dt = (typeRNum)0.001;
   ctypeRNum Thor = 6.0;

   ctypeRNum x0[NX] = {-1,-1};
   ctypeRNum u0[NU] = {0};


   /********* setting parameters *********/
   grampc_setparam_real(grampc, "dt", dt);
   grampc_setparam_real(grampc, "Thor", Thor);
   grampc_setparam_real_vector(grampc, "x0", x0);
   grampc_setparam_real_vector(grampc, "u0", u0);
   ...


Config Reader
*************

.. versionadded:: v2.3

GRAMPC also supplies a config reader which reads ``.cfg`` files for setting the parameters and options inside the ``grampc`` struct.
The config file looks like this

.. literalinclude:: ../../../examples/Crane_2D/config_CRANE_2D.cfg
    :caption: Config file for Crane-2D example

The structure is as follows:
Whitespace or tabs are trimmed.
Therefore one can indent the parameter names or include spaces in the ``.cfg`` file.
A parameter name follows a ``=`` sign, which then follows the parameter value.
String values are plain text without whitespace.
Vectors are specified with square brackets.
The delimiter for the values is either a ``,`` or whitespace.

The parameters must be listed under ``[GRAMPC parameter]``.
However the order is not important.
The same applies for the options which must be listed under ``[GRAMPC option]``.

Calling the config reader is simple with

.. code-block:: c

    /********* Get and set options and parameter from configuration file *********/
    const char *fileName = "config_CRANE_2D.cfg";
    grampc_get_config_from_file(grampc, fileName);

and it sets all defined parameters and options for GRAMPC.
For further examples please see the Crane-2D, Crane-3D, DAE-Integrator, Double-Integrator and the Template folder.

.. _sec:CompileRun:

Compiling and calling GRAMPC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GRAMPC can be integrated into an executable program after the problem formulation (cf. :ref:`chap:ProblemFormulation`) as well as options and parameters are provided. 
A Makefile for compilation purposes is provided in the folder ``<grampc_root>/examples/TEMPLATES``. The main calling routine for GRAMPC is

.. code-block:: c

   void grampc_run(const typeGRAMPC *grampc)

The following code demonstrates how to integrate GRAMPC within an MPC loop. 
For a full implementation please refer to the examples inside ``<grampc_root>/examples/``.

| **Example** **(C code for running GRAMPC within an MPC loop)**

.. code-block:: c

    ...
    /*initialize grampc struct and set parameters and options*/
    ...
    // estimate minimum penalty parameter
    grampc_estim_penmin(grampc, 1);

    /* MPC loop */
    for (iMPC = 0; iMPC <= MaxSimIter; iMPC++) {
        grampc_run(grampc);

        /* check solver status */
        if (grampc->sol->status > 0) {
            if (grampc_printstatus(grampc->sol->status, STATUS_LEVEL_ERROR)) {
                myPrint("at iteration %i:\n -----\n", iMPC);
            }
        }

        /* reference integration or interface to real world measurements */
        get_new_state(xnext);

        /* update state and time */
        t = t + dt;
        grampc_setparam_real_vector(grampc, "x0", xnext);
    }
    ...

GRAMPC is repetitively executed until a defined simulation time is reached. 
Note that the estimate of the minimal penalty value :math:`c_\text{min}` is determined via the function ``grampc_estim_penmin`` before the augmented Lagrangian algorithm is executed, also see the discussion in :ref:`sec:AlgOpt:EstimPenMin`.

For a more convenience MPC design, there are several status flags that can be printed after each MPC step. 
As shown in the above example, the variable ``grampc->sol->status`` is used to check whether new status informations are available.
Subsequently, the function ``grampc_printstatus`` is used to visualize the corresponding informations on the console. 
In the current GRAMPC version, there are status informations concerning the integration scheme, the update of Lagrange multipliers, convergence properties of the augmented Lagrangian and gradient algorithm, and the line search method, see :ref:`sec:AlgOpt:StatusFlags` for more details. 
In addition, the GRAMPC examples provide functions to print key variables such as the system state or the controls into text files. 

The examples can be compiled by running the following commands in a (Cygwin) terminal:

::

    $ cd <grampc_root>/examples/*
    $ make

The make command compiles the file ``main_*.c`` and links it against the GRAMPC toolbox, i.e., against the GRAMPC library within ``<grampc_root>/libs``. 
Note that compiling an application example in the folder ``<grampc root>/examples/`` requires the previous compilation of the GRAMPC toolbox as described in :ref:`chap:install`.
As a result, the executable ``startMPC`` is generated, which can now be used to solve and/or design the MPC problem.

With CMake, one can specify the build target with

::

    cmake --build . --target *

Note that CMake outputs the executables for the examples inside ``<grampc_root>/build/examples/``.

Using GRAMPC without dynamic memory allocation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GRAMPC restricts the usage of dynamic memory allocation, i.e. ``calloc`` and ``realloc``, to initialization and option setting, such that no allocations are performed while the MPC is running. 
However, some microcontrollers and embedded devices do not support dynamic memory allocation at all. 
Using GRAMPC on these devices requires to replace all dynamically-sized arrays by fixed-size arrays and to remove all functions that involve dynamic memory allocation.

To this end, a header file ``fixedsize_settings.h`` must be placed in the search path of the compiler. 
This header file defines several constant parameters such as the number of states ``NX``, the number of controls ``NU``, as well as several constant options, e.g. the number of discretization steps ``NHOR``, the number of gradient iterations ``MAXGRADITER``, and the number of multiplier iterations ``MAXMULTITER``. 
Note that these options cannot be changed during run-time, but are fixed at compile-time. 
A template file is provided in the folder ``<grampc_root>/examples/TEMPLATES``.

In addition, the GRAMPC structures must be created on the stack instead of the heap. 
To this end, the macro ``TYPE_GRAMPC_POINTER`` is provided that allows to use the same code both with and without dynamic memory allocation:

.. code-block:: c

   int main()
   {
     TYPE_GRAMPC_POINTER(grampc);
     ...
     grampc_init(&grampc, userparam);
     ...
     grampc_free(&grampc);
   }

The makefile can be called with the parameter ``FIXEDSIZE=1``, which automatically defines the required preprocessor macro ``FIXEDSIZE``. 
Thus compiling and running the ball-on-plate example without dynamic memory allocation is done by executing the following commands in the terminal:

::

    $ cd <grampc_root>/examples/*
    $ make FIXEDSIZE=1
    $ ./startMPC_fixedsize

Note that in this case a separate GRAMPC library is created in the problem folder that depends on the constants defined in ``fixedsize_settings.h``. 
Changing these settings requires to recompile both the library and the application.

With CMake one can directly specify a fixed-size target in the respective ``CMakeLists.txt`` file. 
For the GRAMPC examples a ``_fixedsize`` is appended to the target name like

::

    cmake --build . --target *_fixedsize

which directly compiles the main and GRAMPC library with the fixed-size settings.

.. footbibliography::
