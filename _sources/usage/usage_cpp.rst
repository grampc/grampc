Using GRAMPC in C++
-------------------

The C++ interface provides a single class for managing the ``grampc`` struct and all related functions as outlined in :ref:`chap:UsageC`.
The main difference is the probfct handling which is described below.

Interface from C++ to C
~~~~~~~~~~~~~~~~~~~~~~~

The constructor of the class ``Grampc`` takes a pointer to a specific problem description. 
This problem description is of type ``ProblemDescription`` and holds all methods which are needed to implement the OCP.
The constructor of ``Grampc`` now passes the ``ProblemDescription*`` as userparam to ``grampc_init``.
Whenever GRAMPC calls a problem specific function like ``ffct``, the ``userparam`` pointer gets dereferenced as a ``ProblemDescription*`` which calls the corresponding method. 

.. code-block:: cpp
    :caption: Example how the bridge between C -> C++ is made.

    void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
    {
        ((grampc::ProblemDescription*)userparam)->ffct(out, t, x, u, p, param);
    }

.. _sec:CppProblemDescription:

Handling of Problem Descriptions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The C++ interface provides a virtual base class ``ProblemDescription`` under the namespace ``grampc``. 
It provides virtual methods which map the functions from ``probfct.h``.

.. note::

    ``userparam`` is not needed, because problem specific quantities can be defined as class fields.

For an example please see the folder ``<grampc_root>/cpp/examples/MassSpringDamper`` with an implementation of the mass spring damper example.
Additionally, compare with the C implementation in ``<grampc_root>/examples/MassSpringDamper``.
