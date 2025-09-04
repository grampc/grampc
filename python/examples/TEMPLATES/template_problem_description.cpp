/* This file is part of GRAMPC - (https://github.com/grampc/grampc)
 *
 * GRAMPC -- A software framework for embedded nonlinear model predictive
 * control using a gradient-based augmented Lagrangian approach
 *
 * Copyright 2014-2025 by Knut Graichen, Andreas Voelz, Thore Wietzke,
 * Tobias Englert (<v2.3), Felix Mesmer (<v2.3), Soenke Rhein (<v2.3),
 * Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 * All rights reserved.
 *
 * GRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
 *
 */
#include "template_problem_description.hpp"

TemplateProblemDescription::TemplateProblemDescription()
 : PyProblemDescription(/*Nx*/ -1, /*Nu*/ -1, /*Np*/ -1, /*Ng*/ -1, /*Nh*/ -1, /*NgT*/ -1, /*NhT*/ -1)
{}
/** System function f(t,x,u,p)
------------------------------------ **/
void TemplateProblemDescription::ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
{

}
/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void TemplateProblemDescription::dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{

}
/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void TemplateProblemDescription::dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{

}
/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
void TemplateProblemDescription::dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{

}


/** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam)
-------------------------------------------------- **/
void TemplateProblemDescription::lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
{

}
/** Gradient dl/dx **/
void TemplateProblemDescription::dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
{

}
/** Gradient dl/du **/
void TemplateProblemDescription::dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
{

}
/** Gradient dl/dp **/
void TemplateProblemDescription::dldp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
{

}


/** Terminal cost V(T,x,p) */
void TemplateProblemDescription::Vfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param)
{

}
/** Gradient dV/dx **/
void TemplateProblemDescription::dVdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param)
{

}
/** Gradient dV/dp **/
void TemplateProblemDescription::dVdp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param)
{

}
/** Gradient dV/dT **/
void TemplateProblemDescription::dVdT(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param)
{

}


/** Equality constraints g(t,x,u,p) = 0 */
void TemplateProblemDescription::gfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
{

}
/** Jacobian dg/dx multiplied by vector vec, i.e. (dg/dx)^T*vec or vec^T*(dg/dx) **/
void TemplateProblemDescription::dgdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{

}
/** Jacobian dg/du multiplied by vector vec, i.e. (dg/du)^T*vec or vec^T*(dg/du) **/
void TemplateProblemDescription::dgdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{

}
/** Jacobian dg/dp multiplied by vector vec, i.e. (dg/dp)^T*vec or vec^T*(dg/dp) **/
void TemplateProblemDescription::dgdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{

}


/** Inequality constraints h(t,x,u,p) < 0 */
void TemplateProblemDescription::hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
{

}
/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
void TemplateProblemDescription::dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{

}
/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
void TemplateProblemDescription::dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{

}
/** Jacobian dh/dp multiplied by vector vec, i.e. (dh/dp)^T*vec or vec^T*(dg/dp) **/
void TemplateProblemDescription::dhdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{

}


/** Terminal equality constraints gT(T,x,p) = 0 */
void TemplateProblemDescription::gTfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param)
{

}
/** Jacobian dgT/dx multiplied by vector vec, i.e. (dgT/dx)^T*vec or vec^T*(dgT/dx) **/
void TemplateProblemDescription::dgTdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{

}
/** Jacobian dgT/dp multiplied by vector vec, i.e. (dgT/dp)^T*vec or vec^T*(dgT/dp) **/
void TemplateProblemDescription::dgTdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{

}
/** Jacobian dgT/dT multiplied by vector vec, i.e. (dgT/dT)^T*vec or vec^T*(dgT/dT) **/
void TemplateProblemDescription::dgTdT_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{

}


/** Terminal inequality constraints hT(T,x,p) < 0 */
void TemplateProblemDescription::hTfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param)
{

}
/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
void TemplateProblemDescription::dhTdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{

}
/** Jacobian dhT/dp multiplied by vector vec, i.e. (dhT/dp)^T*vec or vec^T*(dhT/dp) **/
void TemplateProblemDescription::dhTdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{

}
/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
void TemplateProblemDescription::dhTdT_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{

}


/** Additional functions required for semi-implicit systems
M*dx/dt(t) = f(t0+t,x(t),u(t),p) using the solver RODAS
------------------------------------------------------- **/
/** Jacobian df/dx in vector form (column-wise) **/
void TemplateProblemDescription::dfdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
{

}
/** Jacobian df/dx in vector form (column-wise) **/
void TemplateProblemDescription::dfdxtrans(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
{

}
/** Jacobian df/dt **/
void TemplateProblemDescription::dfdt(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
{

}
/** Jacobian d(dH/dx)/dt  **/
void TemplateProblemDescription::dHdxdt(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{

}
/** Mass matrix in vector form (column-wise, either banded or full matrix) **/
void TemplateProblemDescription::Mfct(VectorRef out, const GrampcParam& param)
{

}
/** Transposed mass matrix in vector form (column-wise, either banded or full matrix) **/
void TemplateProblemDescription::Mtrans(VectorRef out, const GrampcParam& param)
{

}

PYBIND11_MODULE(template, m)
{
    /** Imports pygrampc so this extensions module knows of the type grampc::PyProblemDescription, otherwise an import error like
     * ImportError: generic_type: type "TemplateProblemDescription" referenced unknown base type "grampc::PyProblemDescription"
     * may occur, if the problem description is imported before pygrampc.
     */
    pybind11::module_::import("pygrampc");

    pybind11::class_<TemplateProblemDescription, PyProblemDescription, std::shared_ptr<TemplateProblemDescription>>(m, "TemplateProblemDescription")
        .def(pybind11::init<>())
        .def_readonly("Nx", &TemplateProblemDescription::Nx)
        .def_readonly("Nu", &TemplateProblemDescription::Nu)
        .def_readonly("Np", &TemplateProblemDescription::Np)
        .def_readonly("Ng", &TemplateProblemDescription::Ng)
        .def_readonly("Nh", &TemplateProblemDescription::Nh)
        .def_readonly("NgT", &TemplateProblemDescription::NgT)
        .def_readonly("NhT", &TemplateProblemDescription::NhT);

}
