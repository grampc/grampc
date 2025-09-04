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
#include "Crane2D.hpp"

Crane2D::Crane2D(std::vector<typeRNum> Q, std::vector<typeRNum> R, ctypeRNum ScaleConstraint, ctypeRNum MaxConstraintHeight, ctypeRNum MaxAngularDeflection)
 : PyProblemDescription(/*Nx*/ 6, /*Nu*/ 2, /*Np*/ 0, /*Ng*/ 0, /*Nh*/ 3, /*NgT*/ 0, /*NhT*/ 0),
   Q_(Q),
   R_(R),
   ScaleConstraint_(ScaleConstraint),
   MaxConstraintHeight_(MaxConstraintHeight),
   MaxAngularDeflection_(MaxAngularDeflection)
{}

/** System function f(t,x,u,p)
------------------------------------ **/
void Crane2D::ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) 
{
    out[0] = x[1];
    out[1] = u[0];
    out[2] = x[3];
    out[3] = u[1];
    out[4] = x[5];
    out[5] = -((9.81 * sin(x[4]) + cos(x[4]) * u[0] + 2 * x[3] * x[5]) / x[2]);
};

/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void Crane2D::dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) 
{
    ctypeRNum sinX = sin(x[4]);
    ctypeRNum cosX = cos(x[4]);
    ctypeRNum g = 9.81;

    out[0] = 0;
    out[1] = vec[0];
    out[2] = (g * sinX + cosX * u[0] + 2 * x[3] * x[5]) * vec[5] / (x[2]*x[2]);
    out[3] = vec[2] - (2 * x[5] * vec[5]) / x[2];
    out[4] = -((g * cosX - sinX * u[0]) * vec[5] / x[2]);
    out[5] = vec[4] - (2 * x[3] * vec[5]) / x[2];
};

/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void Crane2D::dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) 
{
    out[0] = vec[1] - (cos(x[4]) * vec[5]) / x[2];
    out[1] = vec[3];
};


/** Integral cost l(t,x(t),u(t),p,xdes,udes)
-------------------------------------------------- **/
void Crane2D::lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) 
{
    out[0] = 0.0;
    for (typeInt i = 0; i < Q_.size(); i++)
    {
        out[0] += Q_[i] * pow(x[i] - param.xdes[i], 2.0);
    }
    for (typeInt i = 0; i < R_.size(); i++)
    {
        out[0] += R_[i] * pow(u[i] - param.udes[i], 2.0);
    }
};
/** Gradient dl/dx **/
void Crane2D::dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) 
{
    for (typeInt i = 0; i < Q_.size(); i++)
    {
        out[i] = 2 * Q_[i] * (x[i] - param.xdes[i]);
    }
};
/** Gradient dl/du **/
void Crane2D::dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) 
{
    for (typeInt i = 0; i < R_.size(); i++)
    {
        out[i] = 2 * R_[i] * (u[i] - param.udes[i]);
    }
};


/** Inequality constraints h(t,x,u,p) < 0 */
void Crane2D::hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) 
{
    ctypeRNum Position = x[0] + sin(x[4]) * x[2];

    out[0] = cos(x[4]) * x[2] - ScaleConstraint_ * pow(Position, 2.0) - MaxConstraintHeight_;
    out[1] = x[5] - MaxAngularDeflection_;
    out[2] = -x[5] - MaxAngularDeflection_;
};

/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
void Crane2D::dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) 
{
    ctypeRNum tmp = ScaleConstraint_ * (x[0] + sin(x[4]) * x[2]);

    out[0] = -2 * tmp * vec[0];
    out[1] = 0;
    out[2] = (sin(x[4]) * tmp + cos(x[4])) * vec[0];
    out[3] = 0;
    out[4] = (-2 * cos(x[4]) * x[2] * tmp - sin(x[4]) * x[2]) * vec[0];
    out[5] = 0 + vec[1] - vec[2];
};

/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
void Crane2D::dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{
    out.fill(0.0);
};

PYBIND11_MODULE(crane_problem, m)
{
    /** Imports pygrampc so this extensions module knows of the type grampc::PyProblemDescription, otherwise an import error like
     * ImportError: generic_type: type "Crane2D" referenced unknown base type "grampc::PyProblemDescription"
     * may occur, if the problem description is imported before pygrampc.
     */
    pybind11::module_::import("pygrampc");

    pybind11::class_<Crane2D, PyProblemDescription, std::shared_ptr<Crane2D>>(m, "Crane2D")
        .def(pybind11::init<std::vector<typeRNum>, std::vector<typeRNum>, typeRNum, typeRNum, typeRNum>())
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