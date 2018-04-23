/* This file is part of GRAMPC - (https://sourceforge.net/projects/grampc/)
 *
 * GRAMPC -- A software framework for embedded nonlinear model predictive
 * control using a gradient-based augmented Lagrangian approach
 *
 * Copyright (C) 2014-2018 by Tobias Englert, Knut Graichen, Felix Mesmer,
 * Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 * Developed at the Institute of Measurement, Control, and Microtechnology,
 * Ulm University. All rights reserved.
 *
 * GRAMPC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * GRAMPC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GRAMPC. If not, see <http://www.gnu.org/licenses/>
 *
 *
 *
 *
 *
 *
 *
 *
 * This probfct-file describes the CSTR reactor problem from
 * Rothfuss, R., Rudolph, J., Zeitz, M.: Flatness Based Control of a
 * Nonlinear Chemical Reactor Model. Automatica 32(10), 1433-1439 (1996)
 *
 *                                           _T
 *                                          /
 *      min    J(u,p,T;x0) = V(T,x(T),p) + / l(t,x(t),u(t),p) dt
 *   u(.),p,T                            _/
 *                                      0
 *             .
 *      s.t.   x(t) = f(t0+t,x(t),u(t),p), x(0) = x0
 *             u_min <= u(t) <= u_max
 *
 */

#include "probfct.h"

#if USE_typeRNum == USE_FLOAT
#define EXP(a)		expf(a)
#else
#define EXP(a)		expf(a)
#endif 

 /* square macro */
#define POW2(a) ((a)*(a))

/** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng), 
    inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT, typeUSERPARAM *userparam)
{
	*Nx = 4;
	*Nu = 2;
	*Np = 0;
	*Ng = 0;
	*Nh = 0;
	*NgT = 0;
	*NhT = 0;
}


/** System function f(t,x,u,p,userparam) 
    ------------------------------------ **/
void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	ctypeRNum* pSys = (ctypeRNum*)userparam;

	typeRNum k1 = pSys[0] * EXP(-(pSys[2] / (((typeRNum)273.15) + x[2])));
	typeRNum k2 = pSys[1] * EXP(-(pSys[3] / (((typeRNum)273.15) + x[2])));
	typeRNum h = -pSys[7] * (k1*pSys[10] * x[0] + k1 * pSys[11] * x[1] + k2 * x[0] * x[0] * pSys[12]);

	out[0] = (-(k2*x[0] * x[0]) + u[0] * (pSys[9] - x[0]) - k1 * x[0]) / pSys[13];
	out[1] = (k1*x[0] - k1 * x[1] - u[0] * x[1]) / pSys[13];
	out[2] = (h + u[0] * (pSys[8] - x[2]) + pSys[4] * (-x[2] + x[3])) / pSys[13];
	out[3] = (pSys[6] * u[1] + pSys[5] * (x[2] - x[3])) / pSys[13];
}


/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	ctypeRNum* pSys = (ctypeRNum*)userparam;

	typeRNum k1 = pSys[0] * EXP(-(pSys[2] / (((typeRNum)273.15) + x[2])));
	typeRNum k2 = pSys[1] * EXP(-(pSys[3] / (((typeRNum)273.15) + x[2])));
	typeRNum dk1dx3 = EXP(-(pSys[2] / (((typeRNum)273.15) + x[2])))*pSys[0] * pSys[2] / ((((typeRNum)273.15) + x[2])*(((typeRNum)273.15) + x[2]));
	typeRNum dk2dx3 = EXP(-(pSys[3] / (((typeRNum)273.15) + x[2])))*pSys[1] * pSys[3] / ((((typeRNum)273.15) + x[2])*(((typeRNum)273.15) + x[2]));
	typeRNum dhdx1 = -(pSys[7] * (k1*pSys[10] + 2 * k2*pSys[12] * x[0]));
	typeRNum dhdx2 = -(k1*pSys[7] * pSys[11]);
	typeRNum dhdx3 = -(pSys[7] * (dk2dx3*x[0] * x[0] * pSys[12] + dk1dx3 * pSys[10] * x[0] + dk1dx3 * pSys[11] * x[1]));

	out[0] = (-((u[0] + 2 * k2*x[0])*vec[0]) + k1 * (-vec[0] + vec[1]) + dhdx1 * vec[2]) / pSys[13];
	out[1] = (-((k1 + u[0])*vec[1]) + dhdx2 * vec[2]) / pSys[13];
	out[2] = (-(dk2dx3*x[0] * x[0] * vec[0]) - dk1dx3 * (x[0] * (vec[0] - vec[1]) + x[1] * vec[1]) + dhdx3 * vec[2] - (pSys[4] + u[0])*vec[2] + pSys[5] * vec[3]) / pSys[13];
	out[3] = (pSys[4] * vec[2] - pSys[5] * vec[3]) / pSys[13];
}
/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	ctypeRNum* pSys = (ctypeRNum*)userparam;

	out[0] = ((pSys[9] - x[0])*vec[0]) / pSys[13] - (x[1] * vec[1]) / pSys[13] + ((pSys[8] - x[2])*vec[2]) / pSys[13];
	out[1] = (pSys[6] * vec[3]) / pSys[13];
}
/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}


/** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam) 
    -------------------------------------------------- **/
void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	ctypeRNum* pSys = (ctypeRNum*)userparam;
	ctypeRNum* pCost = pSys + 14;

	out[0] = pCost[8] * POW2(u[0] - udes[0])
		+ pCost[9] * POW2(u[1] - udes[1])
		+ pCost[4] * POW2(x[0] - xdes[0])
		+ pCost[5] * POW2(x[1] - xdes[1])
		+ pCost[6] * POW2(x[2] - xdes[2])
		+ pCost[7] * POW2(x[3] - xdes[3]);
}
/** Gradient dl/dx **/
void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	ctypeRNum* pSys = (ctypeRNum*)userparam;
	ctypeRNum* pCost = pSys + 14;

	out[0] = 2 * pCost[4] * (x[0] - xdes[0]);
	out[1] = 2 * pCost[5] * (x[1] - xdes[1]);
	out[2] = 2 * pCost[6] * (x[2] - xdes[2]);
	out[3] = 2 * pCost[7] * (x[3] - xdes[3]);
}
/** Gradient dl/du **/
void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	ctypeRNum* pSys = (ctypeRNum*)userparam;
	ctypeRNum* pCost = pSys + 14;

	out[0] = 2 * pCost[8] * (u[0] - udes[0]);
	out[1] = 2 * pCost[9] * (u[1] - udes[1]);
}
/** Gradient dl/dp **/
void dldp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
}


/** Terminal cost V(T,x(T),p,xdes,userparam) 
    ---------------------------------------- **/
void Vfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
	ctypeRNum* pSys = (ctypeRNum*)userparam;
	ctypeRNum* pCost = pSys + 14;

	out[0] = pCost[0] * POW2(x[0] - xdes[0])
		+ pCost[1] * POW2(x[1] - xdes[1])
		+ pCost[2] * POW2(x[2] - xdes[2])
		+ pCost[3] * POW2(x[3] - xdes[3]);
}
/** Gradient dV/dx **/
void dVdx(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
	ctypeRNum* pSys = (ctypeRNum*)userparam;
	ctypeRNum* pCost = pSys + 14;

	out[0] = 2 * pCost[0] * (x[0] - xdes[0]);
	out[1] = 2 * pCost[1] * (x[1] - xdes[1]);
	out[2] = 2 * pCost[2] * (x[2] - xdes[2]);
	out[3] = 2 * pCost[3] * (x[3] - xdes[3]);
}
/** Gradient dV/dp **/
void dVdp(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
}
/** Gradient dV/dT **/
void dVdT(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
}


/** Equality constraints g(t,x(t),u(t),p,uperparam) = 0 
    --------------------------------------------------- **/
void gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian dg/dx multiplied by vector vec, i.e. (dg/dx)^T*vec or vec^T*(dg/dx) **/
void dgdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dg/du multiplied by vector vec, i.e. (dg/du)^T*vec or vec^T*(dg/du) **/
void dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dg/dp multiplied by vector vec, i.e. (dg/dp)^T*vec or vec^T*(dg/dp) **/
void dgdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}


/** Inequality constraints h(t,x(t),u(t),p,uperparam) <= 0 
    ------------------------------------------------------ **/
void hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
void dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
void dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dh/dp multiplied by vector vec, i.e. (dh/dp)^T*vec or vec^T*(dg/dp) **/
void dhdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}


/** Terminal equality constraints gT(T,x(T),p,uperparam) = 0 
    -------------------------------------------------------- **/
void gTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian dgT/dx multiplied by vector vec, i.e. (dgT/dx)^T*vec or vec^T*(dgT/dx) **/
void dgTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dgT/dp multiplied by vector vec, i.e. (dgT/dp)^T*vec or vec^T*(dgT/dp) **/
void dgTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dgT/dT multiplied by vector vec, i.e. (dgT/dT)^T*vec or vec^T*(dgT/dT) **/
void dgTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}


/** Terminal inequality constraints hT(T,x(T),p,uperparam) <= 0 
    ----------------------------------------------------------- **/
void hTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
void dhTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dhT/dp multiplied by vector vec, i.e. (dhT/dp)^T*vec or vec^T*(dhT/dp) **/
void dhTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
void dhTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}


/** Additional functions required for semi-implicit systems 
    M*dx/dt(t) = f(t0+t,x(t),u(t),p) using the solver RODAS 
    ------------------------------------------------------- **/
/** Jacobian df/dx in vector form (column-wise) **/
void dfdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian df/dx in vector form (column-wise) **/
void dfdxtrans(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian df/dt **/
void dfdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian d(dH/dx)/dt  **/
void dHdxdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *vec, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Mass matrix in vector form (column-wise, either banded or full matrix) **/
void Mfct(typeRNum *out, typeUSERPARAM *userparam)
{
}
/** Transposed mass matrix in vector form (column-wise, either banded or full matrix) **/
void Mtrans(typeRNum *out, typeUSERPARAM *userparam)
{
}