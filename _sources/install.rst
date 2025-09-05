.. _chap:install:

Installation and structure of GRAMPC
====================================


The following subsections describe the installation procedure of
GRAMPC for use in C and Matlab.
:ref:`sec:Structure` presents the principal structure of the toolbox.

.. _sec:installC:

Installation of GRAMPC for use in C with Cygwin
-----------------------------------------------

A convenient way to use GRAMPC in C/C++ under MS Windows is the Linux
environment Cygwin. To install Cygwin, download the setup file from the
web page http://www.cygwin.com/ and follow the installation
instructions. In the installation process when packages can be selected,
you have to choose the ``gcc`` compiler and ``make``. If Cygwin is properly installed,
open a Cygwin terminal and perform the following steps:

#. Download the current version of GRAMPC from https://github.com/grampc/grampc.

#. Unpack the archive to an arbitrary location on your computer. 
   After the unpacking procedure, a new directory with the following subfolders is created:

   -  **cpp**: Interface of GRAMPC to C++.

   -  **doc**: Contains this GRAMPC documentation.

   -  **examples**: This folder contains several executable MPC, MHE and OCP problems as well as templates for implementing your own problems.

   -  **include**: The header files of the GRAMPC project are located in this folder.

   -  **libs**: This folder is only available after compiling the GRAMPC toolbox and contains the GRAMPC library.

   -  **matlab**: Interface of GRAMPC to Matlab/Simulink, also see :ref:`sec:InstallMatlab`.

   -  **src**: The source files of the GRAMPC project are located in this folder.

   The GRAMPC directory additionally contains a makefile for building the GRAMPC toolbox. 
   For the remainder of this manual, the location of the created GRAMPC folder will be denoted by ``<grampc_root>``.

#. Compile GRAMPC by running the following commands in a terminal:

   ::
    
      $ cd <grampc_root> 
      $ make clean all

   The dollar symbol indicates the line prompt of the terminal. 
   The make command compiles the source files and generates the GRAMPC library within ``<grampc_root>/libs``, which can now be used to solve a suitable problem in C. 
   The additional argument ``clean all`` removes previously installed parts of GRAMPC.


.. _sec:installCMake:

Installation of GRAMPC for use in C/C++ with CMake
--------------------------------------------------

GRAMPC also supplies ``CMakeLists.txt`` files for building the GRAMPC library and the examples.
From the command line, if cmake is installed, type

::

    mkdir build
    cd build
    cmake ../
    cmake --build .

and the toolbox and the examples are built.
If only the toolbox shall be build, one can issue

::

    cmake --build . --target grampc

which only compiles the GRAMPC library.

.. _sec:InstallMatlab:

Installation of GRAMPC for use in Matlab
----------------------------------------

GRAMPC requires a C compiler that is supported by Matlab for a direct use. 
Details on supported compilers for the current Matlab version as well as previous releases can be found via the Mathworks homepage. 
The correct linkage of the compiler to Matlab can be checked by typing

::

   >> mex -setup

in the Matlab terminal window and subsequently selecting the
corresponding C compiler. The symbol ``>>`` denotes the Matlab prompt. 
The GRAMPC installation under Matlab proceeds in two steps:

#. After downloading and unpacking GRAMPC as described in :ref:`sec:InstallC`, go to the Matlab directory

   ::
    
      >> cd <grampc_root>/matlab

   which contains the following subfolders:

   -  **bin**: This folder is only available after compiling the GRAMPC toolbox and contains the object files of GRAMPC.

   -  **include**: The header files of the GRAMPC project for the Matlab interface.

   -  **mfiles**: Various auxiliary functions for the Matlab interface.

   -  **src**: C sources files of the Mex files which provide the interface between GRAMPC and Matlab.

   In addition, the subfolder contains the m-file ``make.m`` to start the building process.

#. Build the necessary object files for GRAMPC by executing the make function. 
   The compilation of the source files can be performed with the following options:

   -  :code:`>> make clean` removes all previously built GRAMPC files,

   -  :code:`>> make` creates the necessary object files to use GRAMPC,

   -  :code:`>> make verbose` the object files are created in verbose mode, i.e. additional information regarding the building process are provided during the compilation,

   -  :code:`>> make debug` the debug option creates the object files with additional information for use in debugging,

   -  :code:`>> make debug verbose` activates the debug option as well as the verbose mode.

   Similar to the compiling procedure in C as described in :ref:`sec:InstallC`, the ``make`` command compiles the source files to generate object files within ``<grampc root>/matlab/bin``, which can now be used to solve a suitable problem in Matlab.

.. _sec:InstallPython:

Installation of GRAMPC for use in Python
----------------------------------------

.. versionadded:: v2.3

The Python interface of GRAMPC uses pybind11 https://github.com/pybind/pybind11 and Eigen 3.4 https://eigen.tuxfamily.org/index.php?title=Main_Page. 
First, download Eigen 3.4 and follow the steps in ``INSTALL``, so its header files are available through the ``find_package()`` command in CMake.
Make sure the run the installation for Eigen with administration rights.
Then, the interface is installed through

::

   pip install .

assuming you are within ``<grampc root>``, or directly from github with

::

   pip install git+https://github.com/grampc/grampc .

The Python interface can then be imported with 

.. code-block:: python
    
    import pygrampc


.. attention:: If working on Windows, make sure to use Microsoft Visual Studio Compiler (MSVC), since Python for Windows is compiled with MSVC.


.. _sec:Structure:

Structure of GRAMPC
-------------------

The aim of GRAMPC is to be portable and executable on different
operating systems and hardware devices without the use of external
libraries. After the installation procedure as described in :ref:`sec:InstallC` and :ref:`sec:InstallMatlab`, the
GRAMPC structure shown in :numref:`fig:grampcGeneralStructure`
is available to cope with problems from optimal control, model
predictive control, moving horizon estimation, and parameter
optimization. As illustrated in :numref:`fig:grampcGeneralStructure`, the GRAMPC project is
implemented in plain C with a user-friendly interface to C++, Matlab/Simulink, and dSpace.

A specific problem can be implemented in GRAMPC using the C template
``probfct_TEMPLATE.c`` included in the folder ``<grampc root>/examples/TEMPLATES``. 
A more detailed discussion about this step can be found in :ref:`chap:ProblemFormulation`. The workspace of a GRAMPC project
as well as algorithmic options and parameters are stored by the
structure variable ``grampc``. While several parameter settings are
problem specific and need to be provided, most values are set to their
respective default value, see :ref:`chap:AlgOpt`. The GRAMPC
structure can be manipulated through a generic interface, e.g. in order
to set algorithmic options or parameters for a specific problem without
the need to recompile the ``grampc`` project every time.

.. figure:: img/GeneralStructure.*
    :name: fig:grampcGeneralStructure
    :alt: GRAMPC Structure

    General structure of GRAMPC (gray - C code, white - Matlab code).

A specific example contained in the folder ``<grampc root>/examples`` can be 
compiled in C as well as in Matlab and linked against the GRAMPC toolbox. 
A more detailed discussion on this step can be found in :ref:`chap:grampcStructure`.
