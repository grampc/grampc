Using GRAMPC in Python
----------------------

.. versionadded:: v2.3

The Python interface for GRAMPC is similar to the Matlab interface.
It furthermore offers the possibility to write the problem formulation in Python code without a compilation step for rapid prototyping.

Problem Formulation in C++
~~~~~~~~~~~~~~~~~~~~~~~~~~

The first step to implement an MPC in Python is writing the problem definition.
Doing so in C++ is the preferred way, since the computation speed of GRAMPC can be leveraged.
Writing the problem definition itself does not differ greatly from the C++ interface.
Instead of array pointers, e.g. :code:`typeRNum *out`, matrices from Eigen are used.
They are defined by

.. code-block:: cpp

    typedef Eigen::Matrix<typeRNum, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Matrix<typeRNum, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Ref<Vector> VectorRef;
    typedef const Eigen::Ref<const Vector>& VectorConstRef;

and offer a direct conversion to numpy for the Python interface.
Accessing an element of ``Vector``, ``VectorRef`` or ``VectorConstRef`` is the same as for an array pointer.

The class itself has to inherit from ``PyProblemDescription``.
An example of the problem specific header is given in ``<grampc root>/python/examples/Crane2D`` by

.. code-block:: cpp

    #include "pygrampc_problem_description.hpp" // <pybind11/pybind11.h> is already included here
    #include <vector>
    #include <pybind11/stl.h>
    #include <cmath>

    using namespace grampc;

    class Crane2D : public PyProblemDescription // must be inherited with public
    {
        public: // writing the class fields public reduces overhead code for writing getters and setters
            std::vector<typeRNum> Q_;
            std::vector<typeRNum> R_;
            typeRNum ScaleConstraint_;
            typeRNum MaxConstraintHeight_;
            typeRNum MaxAngularDeflection_;

        public:
            Crane2D(std::vector<typeRNum> Q, std::vector<typeRNum> R, ctypeRNum ScaleConstraint, ctypeRNum MaxConstraintHeight, ctypeRNum MaxAngularDeflection);

            ~Crane2D() {}

Differently from C and the C++ interface, no ``ocp_dim()`` function is present to define the problem dimensions.
Instead, they are passed in the constructor of ``PyProblemDescription`` with e.g. for ``Crane2D``:

.. code-block:: cpp

    Crane2D::Crane2D(std::vector<typeRNum> Q, std::vector<typeRNum> R, ctypeRNum ScaleConstraint, ctypeRNum MaxConstraintHeight, ctypeRNum MaxAngularDeflection)
    : PyProblemDescription(/*Nx*/ 6, /*Nu*/ 2, /*Np*/ 0, /*Ng*/ 0, /*Nh*/ 3, /*NgT*/ 0, /*NhT*/ 0),
      Q_(Q),
      R_(R),
      ScaleConstraint_(ScaleConstraint),
      MaxConstraintHeight_(MaxConstraintHeight),
      MaxAngularDeflection_(MaxAngularDeflection)
    {}

After writing the problem specific code, the Python wrapper with pybind11 has to be written.
The binding code is wrapped in the macro ``PYBIND11_MODULE(module_name, module_reference)`` which is the entry point for the Python interpreter when importing a pybind11 module.
For the Crane2D example it is given by

.. code-block:: cpp

    PYBIND11_MODULE(crane_problem, m)
    {
        /** Imports pygrampc so this extensions module knows of the type grampc::PyProblemDescription, otherwise an import error like
        * ImportError: generic_type: type "Crane2D" referenced unknown base type "grampc::PyProblemDescription"
        * may occur, if the problem description is imported before pygrampc.
        */
        pybind11::module_::import("pygrampc");

        pybind11::class_<Crane2D, PyProblemDescription, std::shared_ptr<Crane2D>>(m, "Crane2D")
        ...

At first, pygrampc is imported so our derived problem definition knows the already bound type ``PyProblemDefinition``.
This is just a convenience.
After that, the Python class definition is written with ``pybind11::class_<>``.
Here, we start with our custom class, then the parent we are inheriting from, and also the used capsule for reference counting.
A shared pointer is advisable since otherwise ``Crane2D`` can be prematurely destructed if no reference in Python is present.
In the end of the line, we pass the module reference, here ``m`` and the class name in Python.

The ``__init__()`` function for Python is defined by

.. code-block:: cpp

    ...
        .def(pybind11::init<std::vector<typeRNum>, std::vector<typeRNum>, typeRNum, typeRNum, typeRNum>())
    ...

where a type list of the C++ constructor is supplied.
Note that for binding ``std::vector<>`` the header ``pybind11/stl.h`` has to be included.

This is the theoretic bare minimum needed to initialize the C++ class in Python and pass to GRAMPC.
It is convenient to also expose class fields or functions to Python.
This can be done with

.. code-block:: cpp

        ...
            .def_readonly("Nx", &Crane2D::Nx)
            .def_readonly("Nu", &Crane2D::Nu)
            .def_readonly("Np", &Crane2D::Np)
            .def_readonly("Ng", &Crane2D::Ng)
            .def_readonly("Nh", &Crane2D::Nh)
            .def_readonly("NgT", &Crane2D::NgT)
            .def_readonly("NhT", &Crane2D::NhT)
            
            // make your custom fields available from python
            .def_readwrite("Q", &Crane2D::Q_)
            .def_readwrite("R", &Crane2D::R_)
            .def_readwrite("MaxAngularDeflection", &Crane2D::MaxAngularDeflection_)
            .def_readwrite("ScaleConstraint", &Crane2D::ScaleConstraint_)
            .def_readwrite("MaxConstraintHeight", &Crane2D::MaxConstraintHeight_);
    }

where in this case only fields are exposed.
For a more detailed guide please look into the pybind11 documentation.

The next step involves writing the ``CMakeLists.txt`` file for compiling our problem definition.
We start with configuring the languages and finding Python, pybind11 and Eigen.

.. code-block:: cmake

    cmake_minimum_required(VERSION 3.15)
    project(crane_problem_project LANGUAGES CXX C)

    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
    find_package(pybind11 CONFIG REQUIRED)

    # Eigen 3.4 is required
    find_package(Eigen3 3.4 REQUIRED NO_MODULE)

Then the path to the GRAMPC and PyGRAMPC header files has to be defined:

.. code-block:: cmake

    # Path to <grampc_root>
    set(GRAMPC_ROOT
        ../../../
    )
    
    include_directories(
        ${GRAMPC_ROOT}include # include directory of GRAMPC
        ${GRAMPC_ROOT}python/include # include directory of PyGRAMPC
    )

Finally the pybind11 module is defined with

.. code-block:: cmake  

    pybind11_add_module(crane_problem MODULE Crane2D.cpp)
    target_link_libraries(crane_problem PRIVATE Eigen3::Eigen) # linking against Eigen is necessary
    install(TARGETS crane_problem DESTINATION .) # copy the .pyd (Windows) or .so (Linux) to the correct installation folder where Python finds the extension

Note that ``pybind11_add_module`` is similar to ``add_executeable`` from CMake. 
Here all source files and the target is defined.
The install command is necessary to move the compiled module into ``site-packages``, so Python can directly import the problem definition.

.. important:: The CMake target in ``pybind11_add_module`` and the module name in ``PYBIND11_MODULE`` has to be same! 

To actually make the problem definition importable by Python, it needs a proper Python module definition.
This is done by supplying a ``pyproject.toml`` file given by

.. code-block:: toml

    [build-system]
    requires = ["scikit-build-core>=0.10", "pybind11"]
    build-backend = "scikit_build_core.build"

    [project]
    name = "crane_problem"
    version = "1.0"
    description="Compiled Crane2D problem for GRAMPC"

    [tool.scikit-build]
    wheel.expand-macos-universal-tags = true

In ``[build-system]`` the requirements for building this Python module are defined, namely pybind11 and scikit-build-core.
scikit-build-core takes care of building the extension with CMake and also installs the module into site-packages.

With that, the Crane2D C++ problem definition can be installed with pip with
::

    $ cd <grampc_root>/python/examples/Crane2D
    $ pip install .

and then imported by

.. code-block:: python

    from crane_problem import Crane2D

For an example of the project layout, please refer to ``<grampc root>/python/examples/Crane2D``.

Problem Formulation in Python
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Python interface also allow to write the problem definition in Python. 

.. attention:: The full speed of GRAMPC via Python is only reachable when writing the problem definition in C++! Writing the problem definition in Python can result in 100x longer computation times or more.

This leverages rapid prototyping without a compilation step, but should not be used for evaluating the computation times of GRAMPC.
PyGRAMPC provides the class ``ProblemDescription`` which redirects the C function calls to Python.
An example for defining the problem definition in Python is given in the DoubleIntegrator example by

.. code-block:: python
    :caption: Example for the ``__init__`` method in Python on the basis of the DoubleIntegrator example.
    :name: lis:PythonInit

    from pygrampc import ProblemDescription

    class DoubleIntegrator(ProblemDescription);
        def __init__(self)
            ProblemDescription.__init__(self, Nx=2, Nu=1, Np=0, Ng=0, Nh=0, NgT=2, NhT=0)

            self.CostIntegral = 0.1
            self.CostTerminal = 1.0

Like in C++, the problem dimensions are set in the constructor of ``ProblemDescription``.
Note that every field has to be set.
Just like in the C++ interface, only ``__init__()``, ``ffct()`` and ``dfdx_vec()`` are mandatory to implement.

A sample function implementation looks like

.. code-block:: python

    def ffct(self, out, t, x, u, p, param):
        out[0] = x[1]
        out[1] = u[0]

where ``t`` is a float and ``param`` the ``grampc_param`` struct.
The other parameters are numpy arrays, which point to memory allocated by GRAMPC.
Thus, the results must explicitly write into the allocated memory so statements like

.. code-block:: python

    out[0] = ...
    out[:] = ...

correctly set values to ``out``.
Note that ``out = out + 1`` computes ``out + 1`` and saves the result in the new variable ``out``, which does not point to the allocated memory from GRAMPC.
For a correct usage, please refer to the examples in ``<grampc root>/python/examples``.

Usage in Python
~~~~~~~~~~~~~~~

The usage of GRAMPC in Python is described via the ``DoubleIntegrator.py`` example in ``<grampc root>/python/examples/DoubleIntegrator``.
We first start with importing the relevant packages

.. code-block:: python

    from pygrampc import ProblemDescription, Grampc, GrampcResults
    import matplotlib.pyplot as plt
    from scipy.integrate import solve_ivp

    # alternatively define your problem definition 
    # as a C++ extension module or in a Python file 
    # and import the corresponding problem definition

Here, the ``solve_ivp`` function from Scipy is used for the reference integration. 
Scipy is easily installed with ``pip install scipy`` into the current environment.
In this example, the problem definition is directly written in the same file for rapid prototyping.
Thus, it does not need to be imported.
After that, GRAMPC is initialized:

.. code-block:: python

    if __name__ == "__main__":

        ...

        options = "DoubleIntegrator.json"

        # initialize problem and GRAMPC
        Problem = DoubleIntegrator()
        grampc = Grampc(Problem, options, plot_prediction=False)

We start with a so-called import guard with :code:`if __name__ == "__main__":` which only executes if the current file is run as a main file.
With that, the problem definition in this file can be imported without executing e.g. our test code defined here.
Then, the ``DoubleIntegrator`` class is initialized.
After that, GRAMPC is initialized with ``Problem`` and also with a path to a json file with problem specific options.
This options file is optional.
We can also turn on prediction plots for debugging the current MPC formulation.

Like in Matlab, the result struct ``GrampcResults`` is supplied, which mimics the result and statistic plots available in Matlab.
It is constructed by

.. code-block:: python

    # construct solution structure
    vec = GrampcResults(grampc, Tsim, plot_results=True, plot_statistics=True)


The specific MPC loop for the Double-Integrator is given by

.. code-block:: python
    
    dt = grampc.param.dt

    for i, t in enumerate(vec.t):
        vec.CPUtime[i] = grampc.run()
        vec.update(grampc, i)

        if i + 1 > len(vec.t) or vec.t[i] > Tsim:
            break

        # simulate system
        sol = solve_ivp(grampc.ffct, [t, t + dt], grampc.param.x0,
                        args=(grampc.sol.unext, grampc.sol.pnext, grampc.param))

        # set current time and state
        grampc.set_param({"x0": sol.y[:, -1],
                         "t0": t + dt})

        if grampc.sol.Tnext <= grampc.param.Tmin + grampc.param.dt and bool(grampc.opt.OptimTime):
            grampc.set_opt({"OptimTime": "off"})
            Tsim = vec.t[i + 1]

        # plots of the grampc predictions
        if i % plotSteps == 0:
            grampc.plot()
            vec.plot()
            if plotPause:
                input("Press Enter to continue...")

First, we call ``grampc.run()`` which returns the computation time in milliseconds.
After that, the solution struct is updated with the current results.
With ``solve_ivp`` the reference integration is carried out, with the wrapped ffct like in Matlab.
Here a custom reference integration or model can be implemented.
``Grampc`` also provides ``set_param`` and ``set_opt`` to change parameters and options in GRAMPC.
You have to always pass a Python dict with the respective key-value pairs like

.. code-block:: python

    parameters = {
        "x0": [1.05, 2.0],
        "u0": [0.05,]
    }
    grampc.set_param(parameters)

In the end of the MPC loop, the prediction, results and statistic plots are plotted, if set to :code:`True`.

To run the example, either use your preferred python code editor, or run directly from the terminal with
::

    $ cd <grampc_root>/python/examples/DoubleIntegrator
    $ python DoubleIntegrator.py

For the usage in Python code, please refer to the examples in ``<grampc root>/python/examples``.
