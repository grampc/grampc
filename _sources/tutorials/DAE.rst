.. _sec:Tut:DAE:

Differential algebraic equations
--------------------------------

This section first introduces the DAE-system, before implementation
aspects regarding the dedicated DAE-solver RODAS are considered.
Finally, the example is evaluated.

.. _problem-formulation-3:

Problem formulation
~~~~~~~~~~~~~~~~~~~

The problem at hand is a toy example to illustrate the functionality of
GRAMPC with regard to the solution of DAEs. It consists of two
differential integrator states and one algebraic state. In addition,
this algebraic state is subject to an equality constraint. The
corresponding MPC problem is given by

.. math::
    :label: eq:DAEExample

   	\min_{u} \quad & J(\mb{x}, \mb{u}) = \int_0^T \frac{1}{2}\left( {\Delta} \mb{x} \mb{Q} {\Delta} \mb{x}^\mathsf{T}+ \mb{u} \mb{R} \mb{u}^\mathsf{T}\right)  \, {\rm d}t

   	\textrm{s.t.} \quad & {\dot x_1}(t) = u_1(t) \,, \quad \quad x_1(0) = x_{1,0}

    & {\dot x_2}(t) = u_2(t) \,, \quad \quad x_2(0) = x_{2,0}

    & \hspace{6.5mm} 0 = x_1(t) + x_2(t) - x_3(t)

    & g(\mb{x}(t)) = x_3(t) - 1 = 0

    & \mb{u}(t) \in \left[\mb{u}_{\min}, \mb{u}_{\max}\right],

where
:math:`\Delta \mb{x} = \mb{x} - \mb{x}_{\mathrm{des}}`
and the weight matrices are chosen as
:math:`\mb{Q}  = \text{diag}(500,0,0)` and
:math:`\mb{R} = \text{diag}(1,1)`, respectively. The target of
the MPC formulation is to steer the first differential state to the
desired value, while remaining on the manifold defined by
:math:`x_1(t) + x_2(t) = 1`. Note that this equation results from
substituting the algebraic equation into the constraint :math:`g(\mb{x}(t))` . 
Even though it would be possible to do this substitution and solve the resulting problem, the
purpose of this example is to illustrate the solution of a DAE.

The DAE given in :math:numref:`eq:DAEExample` can be rewritten
with a mass matrix :math:`\mb{M}`, i.e.

.. math::

    \underbrace{\begin{pmatrix}
   	1 & 0 & 0\\
   	0 & 1 & 0\\
   	0 & 0 & 0\\
   	\end{pmatrix}}_{\mb{M}} \dot{\mb{x}} = \underbrace{\begin{pmatrix}
   										  u_1\\
   										  u_2\\
   										  x_1 + x_2 - x_3\\
   										  \end{pmatrix}}_{\mb{f}(\mb{x}, \mb{u})}.

Clearly, the mass matrix is different from the identity matrix and
singular, i.e. the inverse does not exist and therefore the solver RODAS
is used to integrate the system dynamics as well as the corresponding
adjoint dynamics.

.. _implementation-aspects-1:

Implementation aspects
~~~~~~~~~~~~~~~~~~~~~~

To solve the MPC problem for a DAE, the integrator RODAS has to be used
and some additional options have to be set. Furthermore, some additional
functions have to be implemented in the ``probfct``-file.

The options are described in :ref:`sec:AlgOpt:IntegrationRodas`. In the example at
hand, the right hand side of the system dynamics is not explicitly
dependent on the time :math:`t` and therefore
``IFCN`` is set to zero. The next option concerns
the calculation of :math:`\frac{\partial f}{\partial t}` and
:math:`\frac{\partial^2 H}{\partial x \partial t}`. It can either be set
to zero (i.e. :math:`\text{IDFX}=0`) and finite differences
are utilized or set to one (i.e. :math:`\text{IDFX}=1`) and
the analytical solutions implemented in the functions ``dfdt`` and ``dHdxdt`` are called.
The third option determines if the numerical (i.e. finite differences)
or the analytical solution (i.e. ``dfdx`` and ``dfdxtrans``) is used to compute the Jacobians
:math:`\frac{\partial f}{\partial x}` and
:math:`(\frac{\partial f}{\partial x})^\mathrm{T} = \frac{\partial^2 H}{\partial x \partial \lambda}`.
The next option (``IMAS``) determines if the the mass matrix is equal to the
identity matrix (i.e. :math:`\text{IMAS}=0`) or if it is
specified by the functions ``Mfct`` and ``Mtrans`` (i.e. :math:`\text{IMAS}=1`).
In the current example, the mass matrix is singular (not the identity
matrix) and therefore the option is set to one. The remaining options
regard the size of the Jacobian and the mass matrix. The number of
non-zero lower and upper diagonals of the Jacobian are given by ``MLJAC`` and ``MUJAC`` ,
respectively. In our case, we have a full matrix and therefore set both
options to the system dimension, i.e. :math:`N_{\mb{x}}`. The
only non-zero entries of the mass matrix lie on the main diagonal. Thus,
the corresponding options (i.e. ``MLMAS`` and ``MUMAS``) are set to zero. Note that one has
to be careful, if the Jacobian or mass matrix are sparse, since the
lower and upper diagonals are padded with zeros. This is shown for the
example matrix in :numref:`fig:exampleMassMatrix` with the
corresponding code in the following Example.

.. code-block:: c
    :caption: Example (C-Code for the mass function illustrated in :numref:`fig:exampleMassMatrix`)

   	void Mfct(typeRNum *out, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
   	{
   		/* row 1 */
   		out[0]  = 0;
   		out[1]  = k;
   		out[2]  = k;
   		out[3]  = k;
   		/* row 2 */
   		out[4]  = k;
   		out[5]  = k;
   		out[6]  = k;
   		out[7]  = k;
   		/* row 3 */
   		out[8]  = k;
   		out[9]  = k;
   		out[10] = k;
   		out[11] = k;
   		/* row 4 */
   		out[12] = k;
   		out[13] = k;
   		out[14] = k;
   		out[15] = 0;
   		/* row 5 */
   		out[16] = k;
   		out[17] = k;
   		out[18] = 0;
   		out[19] = 0;		
   	};
   	

.. figure:: ../img/tikz/MassMatrixExample.*
    :name: fig:exampleMassMatrix

    Example mass matrix with the index i showing at which position in the output array (i.e. out[i]) the corresponding value has to be written, cf. the example above.

.. _evaluation-1:

Evaluation
~~~~~~~~~~

.. figure:: ../img/tikz/ExampleDAE.*
   :name: fig:ExampleDAE

   Simulated MPC trajectories for the DAE-example

:numref:`fig:ExampleDAE` shows the simulated trajectories for
three set point changes. The setpoint of the first state :math:`x_1`
changes from 1 to 0 at :math:`0\,\mathrm{s}`, from 0 to 0.5 after
:math:`1\,\mathrm{s}`, and finally from 0.5 to 1 after
:math:`2\,\mathrm{s}`. Due to the algebraic state and the equality
constraint, the trajectory of the second state :math:`x_2` has to be the
mirror image of :math:`x_1` around 0.5, which can be observed in the
upper left plot. The corresponding controls in the upper right plot are
also mirrored. The constraint violation during the simulation is shown
in the lower left plot of :numref:`fig:ExampleDAE`. The allowed
constraint violation was set to 1 x 10\ :sup:`-4`, which is
approximately met. Lastly, the lower right plot shows the original
costs and the augmented costs. Both of which quickly converge to zero after each set
point change (after 0, 1 and :math:`2\,\mathrm{s}`, respectively).

.. footbibliography::