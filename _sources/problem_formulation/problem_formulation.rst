.. _chap:ProblemFormulation:

Problem formulation and implementation
======================================

GRAMPC provides a solver for nonlinear input and state constrained
optimal control problems. It is in particular tailored to the
application of real-time MPC with a moving or shrinking horizon with
focus on a memory and time efficient implementation. Other applications
concern the problem of moving horizon estimation and parameter
estimation.

This chapter describes the underlying optimization problem and how it
can be implemented for a specific application in GRAMPC. In addition,
the parameter structure of GRAMPC is introduced.

.. toctree::
    :maxdepth: 2
    
    optimization_problem
    implementation
