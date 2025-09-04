/* This file is part of GRAMPC - (https://github.com/grampc/grampc)
 *
 * GRAMPC -- A software framework for embedded nonlinear model predictive
 * control using a gradient-based augmented Lagrangian approach
 *
 * Copyright 2014-2025 by Knut Graichen, Andreas Voelz, Thore Wietzke,
 * Tobias Englert (<v2.3), Felix Mesmer (<v2.3), Soenke Rhein (<v2.3),
 * Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 *
 * GRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
 *
 *
 *
 *
 *
 *
 *
 *
 * This probfct-file describes the mass-spring-damper problem from
 * Kapernick, B.: Gradient-Based Nonlinear Model Predictive Control With
 * Constraint Transformation for Fast Dynamical Systems. Dissertation,
 * Ulm University. Shaker, Aachen, Germany (2016)
 *
 *                                           _T
 *		                                    /
 *      min    J(u,p,T;x0) = V(T,x(T),p) + / l(t,x(t),u(t),p) dt
 *   u(.),p,T                            _/
 *                                      0
 *             .
 *      s.t.   x(t) = f(t0+t,x(t),u(t),p), x(0) = x0
 *             u_min <= u(t) <= u_max
 *
 */

#include <exception>
#include <stdexcept>
#include "problem_description_MSD.hpp"

MassSpringDamper::MassSpringDamper(typeInt NN, typeRNum m, typeRNum c, typeRNum d, 
		const std::vector<typeRNum>& Q, const std::vector<typeRNum>& R, const std::vector<typeRNum>& P)
 : grampc::ProblemDescription(),
   NN_(NN),
   Nx_(2*NN),
   Nu_(2),
   m_(m),
   c_(c),
   d_(d),
   Q_(Q),
   R_(R),
   P_(P)
{
    if (Q_.size() != Nx_)
    {
        throw std::length_error("Q has to be length Nx = NN*2.");
    }
    if (R_.size() != Nu_)
    {
        throw std::length_error("R has to be length Nu = 2.");
    }
    if (P_.size() != Nx_)
    {
        throw std::length_error("Q has to be length Nx = NN*2.");
    }
}
/** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng), 
    inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
void MassSpringDamper::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
{
    *Nx = Nx_;
    *Nu = Nu_;
    *Np = 0;
    *Ng = 0;
    *Nh = 0;
    *NgT = 0;
    *NhT = 0;
}
/** System function f(t,x,u,p,param,userparam) **/
void MassSpringDamper::ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param)
{
	for (typeInt k = 0; k < NN_; k++)
	{
		out[k] = x[NN_ + k];
	}

	out[NN_] = -2 * c_ / m_ * x[0] + c_ / m_ * x[1] - 2 * d_ / m_ * x[NN_] + d_ / m_ * x[NN_ + 1] + 1 / m_ * u[0];

	for (typeInt k = 1; k < NN_ - 1; k++)
	{
		out[NN_ + k] = c_ / m_ * x[k - 1] - 2 * c_ / m_ * x[k] + c_ / m_ * x[k + 1] + d_ / m_ * x[NN_ + k - 1] - 2 * d_ / m_ * x[NN_ + k] + d_ / m_ * x[NN_ + k + 1];
	}

	out[2 * NN_ - 1] = c_ / m_ * x[NN_ - 2] - 2 * c_ / m_ * x[NN_ - 1] + d_ / m_ * x[2 * NN_ - 2] - 2 * d_ / m_ * x[2 * NN_ - 1] - 1 / m_ * u[1];
}
/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void MassSpringDamper::dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param)
{
	out[0] = -2 * c_ / m_ * vec[NN_] + c_ / m_ * vec[NN_ + 1];

	for (typeInt k = 1; k < NN_ - 1; k++)
	{
		out[k] = c_ / m_ * vec[NN_ + k - 1] - 2 * c_ / m_ * vec[NN_ + k] + c_ / m_ * vec[NN_ + k + 1];
	}

	out[NN_ - 1] = c_ / m_ * vec[2 * NN_ - 2] - 2 * c_ / m_ * vec[2 * NN_ - 1];

	out[NN_] = vec[0] - 2 * d_ / m_ * vec[NN_] + d_ / m_ * vec[NN_ + 1];

	for (typeInt k = 1; k < NN_ - 1; k++)
	{
		out[NN_ + k] = vec[k] + d_ / m_ * vec[NN_ + k - 1] - 2 * d_ / m_ * vec[NN_ + k] + d_ / m_ * vec[NN_ + k + 1];
	}

	out[2 * NN_ - 1] = vec[NN_ - 1] + d_ / m_ * vec[2 * NN_ - 2] - 2 * d_ / m_ * vec[2 * NN_ - 1];

}
/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void MassSpringDamper::dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param)
{
	out[0] = 1 / m_ * vec[NN_];
	out[1] = -1 / m_ * vec[2 * NN_ - 1];
}
/** Integral cost l(t,x(t),u(t),p,param,userparam) **/
void MassSpringDamper::lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param)
{
	out[0] = 0;
	for (typeInt i = 0; i < Nx_; i++)
	{
		out[0] += Q_[i] * POW2(x[i] - param->xdes[i]);
	}
	for (typeInt i = 0; i < Nu_; i++)
	{
		out[0] += R_[i] * POW2(u[i] - param->udes[i]);
	}
	out[0] *= 0.5;
}
/** Gradient dl/dx **/
void MassSpringDamper::dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param)
{
	for (typeInt i = 0; i < Nx_; i++)
	{
		out[i] = Q_[i] * (x[i] - param->xdes[i]);
	}
}
/** Gradient dl/du **/
void MassSpringDamper::dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param)
{
	for (typeInt i = 0; i < Nu_; i++)
	{
		out[i] = R_[i] * (u[i] - param->udes[i]);
	}
}
/** Terminal cost V(T,x(T),p,param,userparam) **/
void MassSpringDamper::Vfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param)
{
	out[0] = 0;
	for (typeInt i = 0; i < Nx_; i++)
	{
		out[0] += P_[i] * POW2(x[i] - param->xdes[i]);
	}
	out[0] *= 0.5;
}
/** Gradient dV/dx **/
void MassSpringDamper::dVdx(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param)
{
	for (typeInt i = 0; i < Nx_; i++)
	{
		out[i] = P_[i] * (x[i] - param->xdes[i]);
	}
}