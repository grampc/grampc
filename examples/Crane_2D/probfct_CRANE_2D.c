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
 * This probfct-file describes the crane 2D problem from
 * Kapernick, B., Graichen, K.: Model predictive control of an overhead
 * crane using constraint substitution. In: Proc. American Control 
 * Conference (ACC), pp. 3973-3978, 2013.
 *
 *                                           _T
 *                                          /
 *      min    J(u,p,T;x0) = V(T,x(T),p) + / l(t,x(t),u(t),p) dt
 *   u(.),p,T                            _/
 *                                      0
 *             .
 *      s.t.   x(t) = f(t0+t,x(t),u(t),p), x(0) = x0
 *             h(x)  <= 0
 *             u_min <= u(t) <= u_max
 *
 */

#include "probfct.h"

 /* Data type-dependent functions */
#if USE_typeRNum == USE_FLOAT
#define SIN(a)		sinf(a)
#define COS(a)		cosf(a)
#else
#define SIN(a)		sin(a)
#define COS(a)		cos(a)
#endif 

 /* square macro */
#define POW2(a) ((a)*(a))

 /* gravitation constant */
#define GRAV ((typeRNum) 9.81)

/** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng), 
    inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT, typeUSERPARAM *userparam)
{
	*Nx = 6;
	*Nu = 2;
	*Np = 0;
	*Nh = 3;
	*Ng = 0;
	*NgT = 0;
	*NhT = 0;
}


/** System function f(t,x,u,p,param,userparam) 
    ------------------------------------ **/
void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
	out[0] = x[1];
	out[1] = u[0];
	out[2] = x[3];
	out[3] = u[1];
	out[4] = x[5];
	out[5] = -((GRAV * SIN(x[4]) + COS(x[4])*u[0] + 2 * x[3] * x[5]) / x[2]);
}
/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
	ctypeRNum aux1 = SIN(x[4]);
	ctypeRNum aux2 = COS(x[4]);

	out[0] = 0;
	out[1] = vec[0];
	out[2] = (GRAV * aux1 + aux2 * u[0] + 2 * x[3] * x[5])*vec[5] / (x[2] * x[2]);
	out[3] = vec[2] - (2 * x[5] * vec[5]) / x[2];
	out[4] = -(((GRAV * aux2 - aux1 * u[0])*vec[5]) / x[2]);
	out[5] = vec[4] - (2 * x[3] * vec[5]) / x[2];
}
/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
	out[0] = vec[1] - ((COS(x[4]))*vec[5]) / x[2];
	out[1] = vec[3];
}
/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}


/** Integral cost l(t,x(t),u(t),p,param,userparam) 
    -------------------------------------------------- **/
void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
	ctypeRNum *pCost = (ctypeRNum*)userparam;
    ctypeRNum* xdes = param->xdes;
    ctypeRNum* udes = param->udes;

	out[0] = pCost[0] * POW2(x[0] - xdes[0]) +
		pCost[1] * POW2(x[1] - xdes[1]) +
		pCost[2] * POW2(x[2] - xdes[2]) +
		pCost[3] * POW2(x[3] - xdes[3]) +
		pCost[4] * POW2(x[4] - xdes[4]) +
		pCost[5] * POW2(x[5] - xdes[5]) +
		pCost[6] * POW2(u[0] - udes[0]) +
		pCost[7] * POW2(u[1] - udes[1]);
}
/** Gradient dl/dx **/
void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
	ctypeRNum *pCost = (ctypeRNum*)userparam;
    ctypeRNum* xdes = param->xdes;

	out[0] = 2 * pCost[0] * (x[0] - xdes[0]);
	out[1] = 2 * pCost[1] * (x[1] - xdes[1]);
	out[2] = 2 * pCost[2] * (x[2] - xdes[2]);
	out[3] = 2 * pCost[3] * (x[3] - xdes[3]);
	out[4] = 2 * pCost[4] * (x[4] - xdes[4]);
	out[5] = 2 * pCost[5] * (x[5] - xdes[5]);
}
/** Gradient dl/du **/
void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
	ctypeRNum *pCost = (ctypeRNum*)userparam;
    ctypeRNum* udes = param->udes;

	out[0] = 2 * pCost[6] * (u[0] - udes[0]);
	out[1] = 2 * pCost[7] * (u[1] - udes[1]);
}
/** Gradient dl/dp **/
void dldp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}


/** Terminal cost V(T,x(T),p,param,userparam) 
    ---------------------------------------- **/
void Vfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Gradient dV/dx **/
void dVdx(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Gradient dV/dp **/
void dVdp(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Gradient dV/dT **/
void dVdT(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}


/** Equality constraints g(t,x(t),u(t),p,param,userparam) = 0 
    --------------------------------------------------- **/
void gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Jacobian dg/dx multiplied by vector vec, i.e. (dg/dx)^T*vec or vec^T*(dg/dx) **/
void dgdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Jacobian dg/du multiplied by vector vec, i.e. (dg/du)^T*vec or vec^T*(dg/du) **/
void dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Jacobian dg/dp multiplied by vector vec, i.e. (dg/dp)^T*vec or vec^T*(dg/dp) **/
void dgdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}


/** Inequality constraints h(t,x(t),u(t),p,param,userparam) <= 0 
    ------------------------------------------------------ **/
void hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
	ctypeRNum *pSys = (ctypeRNum*)userparam + 8;
	ctypeRNum xPos = x[0] + SIN(x[4])*x[2];

	out[0] = COS(x[4]) * x[2] - pSys[0] * POW2(xPos) - pSys[1];
	out[1] = x[5] - pSys[2];
	out[2] = -x[5] - pSys[2];
}
/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
void dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
	ctypeRNum *pSys = (ctypeRNum*)userparam + 8;
	ctypeRNum temp = pSys[0] * (x[0] + SIN(x[4])*x[2]);

	out[0] = -2 * temp * vec[0];
	out[1] = 0;
	out[2] = (SIN(x[4])*temp + COS(x[4]))*vec[0];
	out[3] = 0;
	out[4] = (-2 * temp * COS(x[4])*x[2] - SIN(x[4])*x[2])*vec[0];
	out[5] = 0 + vec[1] - vec[2];
}
/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
void dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
	out[0] = 0;
}
/** Jacobian dh/dp multiplied by vector vec, i.e. (dh/dp)^T*vec or vec^T*(dg/dp) **/
void dhdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}


/** Terminal equality constraints gT(T,x(T),p,param,userparam) = 0 
    -------------------------------------------------------- **/
void gTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Jacobian dgT/dx multiplied by vector vec, i.e. (dgT/dx)^T*vec or vec^T*(dgT/dx) **/
void dgTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Jacobian dgT/dp multiplied by vector vec, i.e. (dgT/dp)^T*vec or vec^T*(dgT/dp) **/
void dgTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Jacobian dgT/dT multiplied by vector vec, i.e. (dgT/dT)^T*vec or vec^T*(dgT/dT) **/
void dgTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}


/** Terminal inequality constraints hT(T,x(T),p,param,userparam) <= 0 
    ----------------------------------------------------------- **/
void hTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
void dhTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Jacobian dhT/dp multiplied by vector vec, i.e. (dhT/dp)^T*vec or vec^T*(dhT/dp) **/
void dhTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
void dhTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}


/** Additional functions required for semi-implicit systems 
    M*dx/dt(t) = f(t,x(t),u(t),p,param,userparam) using the solver RODAS 
    ------------------------------------------------------- **/
/** Jacobian df/dx in vector form (column-wise) **/
void dfdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Jacobian df/dx in vector form (column-wise) **/
void dfdxtrans(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Jacobian df/dt **/
void dfdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Jacobian d(dH/dx)/dt  **/
void dHdxdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Mass matrix in vector form (column-wise, either banded or full matrix) **/
void Mfct(typeRNum *out, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Transposed mass matrix in vector form (column-wise, either banded or full matrix) **/
void Mtrans(typeRNum *out, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
