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
 *
 *
 *
 * This probfct-file describes the quadrotor example originally proposed in
 * 
 * Hehn M, D’Andrea R. A flying inverted pendulum. IEEE International 
 * Conference on Robotics and Automation, Shanghai (China), 2011; 763–770.
 *
 * Also see
 *
 * B. Kapernick, K. Graichen. Nonlinear model predictive control based 
 * on constraint transformation. Optimal Control Applications and Methods, 
 * 37(4):807–828, 2016.
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
#define SIN(a)		sinf(a)
#define COS(a)		cosf(a)
#define TAN(a)		tanf(a)
#else
#define SIN(a)		sin(a)
#define COS(a)		cos(a)
#define TAN(a)		tan(a)
#endif

 /* square macro */
#define POW2(a) ((a)*(a))

 /* gravitation constant */
#define GRAV ((typeRNum) 9.81)

 /** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng), 
    inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT, typeUSERPARAM *userparam)
{
	*Nx = 9;
	*Nu = 4;
	*Np = 0;
	*Nh = 0;
	*Ng = 0;
	*NgT = 0;
	*NhT = 0;
}


/** System function f(t,x,u,p,userparam) 
    ------------------------------------ **/
void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	typeRNum aux1 = COS(x[8]);
	typeRNum aux2 = SIN(x[6]);
	typeRNum aux3 = COS(x[6]);
	typeRNum aux4 = SIN(x[7]);
	typeRNum aux5 = SIN(x[8]);
	typeRNum aux6 = COS(x[7]);
	typeRNum aux7 = TAN(x[7]);

	out[0] = x[1];
	out[1] = u[0] * (aux3*aux1*aux4 + aux5 * aux2);
	out[2] = x[3];
	out[3] = u[0] * (aux5*aux3*aux4 - aux1 * aux2);
	out[4] = x[5];
	out[5] = ((typeRNum)-9.81) + u[0] * aux6*aux3;
	out[6] = (1 * (u[1] * aux3 + u[2] * aux2)) / aux6;
	out[7] = u[2] * aux3 - u[1] * aux2;
	out[8] = u[1] * aux3*aux7 + u[2] * aux7*aux2 + u[3];
}
/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	typeRNum aux1 = COS(x[7]);
	typeRNum aux2 = 1 / aux1;
	typeRNum aux3 = TAN(x[7]);
	typeRNum aux4 = COS(x[8]);
	typeRNum aux5 = SIN(x[8]);
	typeRNum aux6 = POW2(aux2);
	typeRNum aux7 = SIN(x[6]);
	typeRNum aux8 = COS(x[6]);
	typeRNum aux9 = SIN(x[7]);
	typeRNum aux10 = vec[6] * aux9;
	typeRNum aux11 = vec[8] + aux10;

	out[0] = 0;
	out[1] = vec[0];
	out[2] = 0;
	out[3] = vec[2];
	out[4] = 0;
	out[5] = vec[4];
	out[6] = aux8 * (1 * aux2*u[2] * vec[6] - u[1] * vec[7] + u[2] * aux3*vec[8] - aux4 * vec[3] * u[0] + aux5 * u[0] * vec[1])
		- COS(x[7])*aux7*(1 * aux2*u[2] * vec[7] + vec[6] * u[1] * aux6 + vec[5] * u[0] + aux3 * (1 * aux2*u[1] * vec[8] + vec[3] * aux5*u[0] + aux4 * u[0] * vec[1]));
	out[7] = 1 * u[2] * aux7*aux6*aux11 + aux8 * (1 * u[1] * aux6*aux11 - vec[5] * aux9*u[0]
		+ COS(x[7])*u[0] * (vec[3] * aux5 + aux4 * vec[1]));
	out[8] = u[0] * (aux4*(vec[3] * aux9*aux8 + aux7 * vec[1]) + aux5 * (vec[3] * aux7 - aux9 * aux8*vec[1]));
}
/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	typeRNum aux1 = COS(x[8]);
	typeRNum aux2 = SIN(x[8]);
	typeRNum aux3 = SIN(x[6]);
	typeRNum aux4 = COS(x[6]);
	typeRNum aux5 = COS(x[7]);
	typeRNum aux6 = SIN(x[7]);
	typeRNum aux7 = 1 / aux5;
	typeRNum aux8 = aux6 * vec[8];
	typeRNum aux9 = vec[6] + aux8;

	out[0] = aux3 * (aux2*vec[1] - aux1 * vec[3]) + aux4 * (aux6*(vec[1] * aux1 + aux2 * vec[3]) + aux5 * vec[5]);
	out[1] = -(aux3*vec[7]) + 1 * aux4*aux7*aux9;
	out[2] = aux4 * vec[7] + 1 * aux3*aux7*aux9;
	out[3] = vec[8];
}
/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}


/** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam) 
    -------------------------------------------------- **/
void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	ctypeRNum *pCost = (ctypeRNum*)userparam;

	out[0] = pCost[18] * POW2(u[0] - udes[0])
		+ pCost[19] * POW2(u[1] - udes[1])
		+ pCost[20] * POW2(u[2] - udes[2])
		+ pCost[21] * POW2(u[3] - udes[3])
		+ pCost[9] * POW2(x[0] - xdes[0])
		+ pCost[10] * POW2(x[1] - xdes[1])
		+ pCost[11] * POW2(x[2] - xdes[2])
		+ pCost[12] * POW2(x[3] - xdes[3])
		+ pCost[13] * POW2(x[4] - xdes[4])
		+ pCost[14] * POW2(x[5] - xdes[5])
		+ pCost[15] * POW2(x[6] - xdes[6])
		+ pCost[16] * POW2(x[7] - xdes[7])
		+ pCost[17] * POW2(x[8] - xdes[8]);
}
/** Gradient dl/dx **/
void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	ctypeRNum *pCost = (ctypeRNum*)userparam;

	out[0] = 2 * pCost[9] * (x[0] - xdes[0]);
	out[1] = 2 * pCost[10] * (x[1] - xdes[1]);
	out[2] = 2 * pCost[11] * (x[2] - xdes[2]);
	out[3] = 2 * pCost[12] * (x[3] - xdes[3]);
	out[4] = 2 * pCost[13] * (x[4] - xdes[4]);
	out[5] = 2 * pCost[14] * (x[5] - xdes[5]);
	out[6] = 2 * pCost[15] * (x[6] - xdes[6]);
	out[7] = 2 * pCost[16] * (x[7] - xdes[7]);
	out[8] = 2 * pCost[17] * (x[8] - xdes[8]);
}
/** Gradient dl/du **/
void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	ctypeRNum *pCost = (ctypeRNum*)userparam;

	out[0] = 2 * pCost[18] * (u[0] - udes[0]);
	out[1] = 2 * pCost[19] * (u[1] - udes[1]);
	out[2] = 2 * pCost[20] * (u[2] - udes[2]);
	out[3] = 2 * pCost[21] * (u[3] - udes[3]);
}
/** Gradient dl/dp **/
void dldp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
}


/** Terminal cost V(T,x(T),p,xdes,userparam) 
    ---------------------------------------- **/
void Vfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
	ctypeRNum *pCost = (ctypeRNum*)userparam;

	out[0] = pCost[0] * (x[0] - xdes[0])*(x[0] - xdes[0]) +
		pCost[1] * (x[1] - xdes[1])*(x[1] - xdes[1]) +
		pCost[2] * (x[2] - xdes[2])*(x[2] - xdes[2]) +
		pCost[3] * (x[3] - xdes[3])*(x[3] - xdes[3]) +
		pCost[4] * (x[4] - xdes[4])*(x[4] - xdes[4]) +
		pCost[5] * (x[5] - xdes[5])*(x[5] - xdes[5]) +
		pCost[6] * (x[6] - xdes[6])*(x[6] - xdes[6]) +
		pCost[7] * (x[7] - xdes[7])*(x[7] - xdes[7]) +
		pCost[8] * (x[8] - xdes[8])*(x[8] - xdes[8]);
}
/** Gradient dV/dx **/
void dVdx(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
	ctypeRNum *pCost = (ctypeRNum*)userparam;

	out[0] = 2 * pCost[0] * (x[0] - xdes[0]);
	out[1] = 2 * pCost[1] * (x[1] - xdes[1]);
	out[2] = 2 * pCost[2] * (x[2] - xdes[2]);
	out[3] = 2 * pCost[3] * (x[3] - xdes[3]);
	out[4] = 2 * pCost[4] * (x[4] - xdes[4]);
	out[5] = 2 * pCost[5] * (x[5] - xdes[5]);
	out[6] = 2 * pCost[6] * (x[6] - xdes[6]);
	out[7] = 2 * pCost[7] * (x[7] - xdes[7]);
	out[8] = 2 * pCost[8] * (x[8] - xdes[8]);
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
