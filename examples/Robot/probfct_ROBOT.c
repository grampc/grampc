/* This file is part of GRAMPC - (https://sourceforge.net/projects/grampc/)
 *
 * GRAMPC -- A software framework for embedded nonlinear model predictive
 * control using a gradient-based augmented Lagrangian approach
 *
 * Copyright 2014-2019 by Tobias Englert, Knut Graichen, Felix Mesmer,
 * Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 * All rights reserved.
 *
 * GRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
 *
 *
 *
 *
 *
 * This probfct-file describes a planar dual-arm robot with a closed
 * kinematic chain, see Englert et al.: "A software framework for embedded
 * nonlinear model predictive control using a gradient-based augmented
 * Lagrangian approach (GRAMPC)", submitted to Optimization and Engineering
 *   
 *                                           _T
 *                                          /
 *      min    J(u,p,T;x0) = V(T,x(T),p) + / l(t,x(t),u(t),p) dt
 *   u(.),p,T                            _/
 *                                      0
 *             .
 *      s.t.   x(t) = f(t0+t,x(t),u(t),p),  x(0) = x_start
 *             g(x(t),u(t),p) = 0,          x(T) = x_goal
 *             u_min <= u(t) <= u_max
 *
 */

#include "probfct.h"

#ifndef M_PI
#define M_PI 3.14159
#endif
#define NDOF 6
#define NL 3
#define NR 3
#define NPOSE 3
#define NWEIGHTS 6

 /********* robot specific functions *********/
	/* Forward kinematics of planar robot with N revolute joints */
void kinematics(typeRNum *out, ctypeRNum *x, typeInt N, ctypeRNum *len, ctypeRNum *offset)
{
	typeInt i;
	typeRNum rot = offset[2];

	out[0] = offset[0];
	out[1] = offset[1];
	for (i = 0; i < N; ++i)
	{
		rot += x[i];
		out[0] += len[i] * cos(rot);
		out[1] += len[i] * sin(rot);
	}
	out[2] = rot;
}
/* Kinematic jacobian of planar robot with N revolute joints */
void jacobian(typeRNum *out, ctypeRNum *x, typeInt N, ctypeRNum *len, ctypeRNum *offset)
{
	typeInt i, j;
	for (i = 0; i < N; ++i)
	{
		out[0 * N + i] = 0;
		out[1 * N + i] = 0;
		out[2 * N + i] = 1;
	}
	typeRNum rot = offset[2];
	for (i = 0; i < N; ++i)
	{
		rot += x[i];
		for (j = 0; j <= i; ++j)
		{
			out[0 * N + j] += len[i] * -sin(rot);
			out[1 * N + j] += len[i] * cos(rot);
		}
	}
}


/** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng), 
    inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT, typeUSERPARAM *userparam)
{
	*Nx = NDOF;
	*Nu = NDOF;
	*Np = 0;
	*Ng = NPOSE;
	*Nh = 0;
	*NgT = NDOF;
	*NhT = 0;
}


/** System function f(t,x,u,p,userparam) 
    ------------------------------------ **/
void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	typeInt i;
	/* Single-integrator for each degree of freedom */
	for (i = 0; i < NDOF; ++i)
	{
		out[i] = u[i];
	}
}
/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	typeInt i;
	for (i = 0; i < NDOF; ++i)
	{
		out[i] = 0;
	}
}
/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	typeInt i;
	for (i = 0; i < NDOF; ++i)
	{
		out[i] = vec[i];
	}
}
/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}


/** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam) 
    -------------------------------------------------- **/
void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	ctypeRNum *pCost = (ctypeRNum *)userparam;
	typeInt i;

	out[0] = 0;
	for (i = 0; i < NDOF; ++i)
	{
		out[0] += pCost[2] * (x[i] - xdes[i]) * (x[i] - xdes[i]);
		out[0] += pCost[4] * (u[i] - udes[i]) * (u[i] - udes[i]);
	}
}
/** Gradient dl/dx **/
void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	ctypeRNum *pCost = (ctypeRNum *)userparam;
	typeInt i;

	for (i = 0; i < NDOF; ++i)
	{
		out[i] = 2 * pCost[2] * (x[i] - xdes[i]);
	}
}
/** Gradient dl/du **/
void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	ctypeRNum *pCost = (ctypeRNum *)userparam;
	typeInt i;

	for (i = 0; i < NDOF; ++i)
	{
		out[i] = 2 * pCost[4] * (u[i] - udes[i]);
	}
}
/** Gradient dl/dp **/
void dldp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
}


/** Terminal cost V(T,x(T),p,xdes,userparam) 
    ---------------------------------------- **/
void Vfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
	ctypeRNum *pCost = (ctypeRNum *)userparam;
	typeInt i;

	out[0] = 0;
	for (i = 0; i < NDOF; ++i)
	{
		out[0] += pCost[0] * (x[i] - xdes[i]) * (x[i] - xdes[i]);
	}
	out[0] += pCost[5] * T;
}
/** Gradient dV/dx **/
void dVdx(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
	ctypeRNum *pCost = (ctypeRNum *)userparam;
	typeInt i;

	for (i = 0; i < NDOF; ++i)
	{
		out[i] = 2 * pCost[0] * (x[i] - xdes[i]);
	}
}
/** Gradient dV/dp **/
void dVdp(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
}
/** Gradient dV/dT **/
void dVdT(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
	ctypeRNum *pCost = (ctypeRNum *)userparam;
	out[0] = pCost[5];
}


/** Equality constraints g(t,x(t),u(t),p,uperparam) = 0 
    --------------------------------------------------- **/
void gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	ctypeRNum *pSys = (ctypeRNum *)userparam + NWEIGHTS;
	typeInt i;
	typeRNum kinl[NPOSE];
	typeRNum kinr[NPOSE];
	kinematics(kinl, x, NL, pSys, pSys + NDOF);
	kinematics(kinr, x + NL, NR, pSys + NL, pSys + NDOF + NL);

	for (i = 0; i < NPOSE; ++i)
	{
		out[i] = (kinl[i] - kinr[i]);
	}
	out[2] += M_PI; /* orientation should be inverted between left and right arm */
}
/** Jacobian dg/dx multiplied by vector vec, i.e. (dg/dx)^T*vec or vec^T*(dg/dx) **/
void dgdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
	ctypeRNum *pSys = (ctypeRNum *)userparam + NWEIGHTS;
	typeRNum jacl[NPOSE*NL];
	typeRNum jacr[NPOSE*NR];
	typeInt i, j;
	jacobian(jacl, x, NL, pSys, pSys + NDOF);
	jacobian(jacr, x + NL, NR, pSys + NL, pSys + NDOF + NL);

	for (j = 0; j < NL; ++j)
	{
		out[j] = 0;
		for (i = 0; i < NPOSE; ++i)
		{
			out[j] += jacl[i*NL + j] * vec[i];
		}
	}
	for (j = 0; j < NR; ++j)
	{
		out[NL + j] = 0;
		for (i = 0; i < NPOSE; ++i)
		{
			out[NL + j] += -jacr[i*NR + j] * vec[i];
		}
	}
}
/** Jacobian dg/du multiplied by vector vec, i.e. (dg/du)^T*vec or vec^T*(dg/du) **/
void dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
	typeInt i;
	for (i = 0; i < NDOF; ++i)
	{
		out[i] = 0;
	}
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
	ctypeRNum *xdes = (ctypeRNum *)userparam + NWEIGHTS + 2 * NDOF;
	typeInt i;
	for (i = 0; i < NDOF; ++i)
	{
		out[i] = x[i] - xdes[i];
	}
}
/** Jacobian dgT/dx multiplied by vector vec, i.e. (dgT/dx)^T*vec or vec^T*(dgT/dx) **/
void dgTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
	typeInt i;
	for (i = 0; i < NDOF; ++i)
	{
		out[i] = vec[i];
	}
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