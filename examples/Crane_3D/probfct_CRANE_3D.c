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
 * This probfct-file describes the crane 3D problem from
 * Graichen, K., Egretzberger, M., Kugi, A.: A Suboptimal Approach to Real-Time
 * Model Predictive Control of Nonlinear Systems. Automatisierungstechnik (2010)
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

 /* definitions of the wieghts' position in the userparam */
#define iQ    0           /* position of integral weights (inputs) in pCost */
#define iR    iQ+10       /* position of integral weights (states) in pCost */
#define iP    iR+3        /* position of terminal weights (states) in pCost */

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
	*Nx = 10;
	*Nu = 3;
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
	ctypeRNum t1 = 1 / x[1];
	ctypeRNum t2 = COS(x[4]);
	ctypeRNum t5 = x[6] * x[7];
	ctypeRNum t6 = SIN(x[4]);
	ctypeRNum t7 = COS(x[3]);
	ctypeRNum t8 = t6 * t7;
	ctypeRNum t15 = t2 * x[1];
	ctypeRNum t23 = POW2(x[7]);
	ctypeRNum t24 = t23 * t2;
	ctypeRNum t25 = x[1] * t7;
	ctypeRNum t26 = SIN(x[3]);
	ctypeRNum t29 = t23 * x[0];
	ctypeRNum t45 = t7 * t7;
	ctypeRNum t49 = t2 * t2;
	ctypeRNum t54 = POW2(x[8]);
	ctypeRNum t57 = t26 * t6;
	ctypeRNum t65 = -2 * x[5] * x[7] * t2 - 2 * x[6] * x[9] + 2 * t5 * t26 + t24 * x[1] * t45 * t6 + 2 * x[7] * t49 * t25 * x[8] - t15 * t54 * t6 + t29 * t57 - t8 * GRAV - t57 * u[0] - u[2] * x[0] * t2 + u[2] * x[1] * t26;

	out[0] = x[5];
	out[1] = x[6];
	out[2] = x[7];
	out[3] = x[8];
	out[4] = x[9];
	out[5] = u[0];
	out[6] = u[1];
	out[7] = u[2];
	out[8] = -t1 / t2 * (2 * t5 * t8 + 2 * x[6] * t2 * x[8] + 2 * x[9] * x[7] * t15 * t7 - 2 * x[9] * x[1] * x[8] * t6 - t24 * t25 * t26 + t29 * t7 + t26 * GRAV - t7 * u[0] + t6 * x[1] * t7 * u[2]);
	out[9] = t1 * t65;
}
/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	ctypeRNum t1 = 1 / x[1];
	ctypeRNum t2 = COS(x[4]);
	ctypeRNum t3 = 1 / t2;
	ctypeRNum t5 = POW2(x[7]);
	ctypeRNum t6 = COS(x[3]);
	ctypeRNum t10 = SIN(x[3]);
	ctypeRNum t12 = SIN(x[4]);
	ctypeRNum t19 = x[6] * x[7];
	ctypeRNum t20 = t12 * t6;
	ctypeRNum t23 = x[6] * t2;
	ctypeRNum t26 = t5 * x[0];
	ctypeRNum t31 = POW2(x[1]);
	ctypeRNum t32 = 1 / t31;
	ctypeRNum t34 = t3 * vec[8];
	ctypeRNum t36 = x[5] * x[7];
	ctypeRNum t43 = t10 * t12;
	ctypeRNum t47 = u[2] * x[0];
	ctypeRNum t60 = t5 * t2;
	ctypeRNum t61 = t60 * x[1];
	ctypeRNum t62 = t6 * t6;
	ctypeRNum t63 = x[1] * t62;
	ctypeRNum t76 = 2 * t19 * t6;
	ctypeRNum t80 = t2 * t2;
	ctypeRNum t82 = x[1] * t10;
	ctypeRNum t86 = t26 * t20;
	ctypeRNum t87 = t43 * GRAV;
	ctypeRNum t88 = t20 * u[0];
	ctypeRNum t90 = u[2] * x[1] * t6;
	ctypeRNum t95 = x[9] * x[1];
	ctypeRNum t110 = x[7] * t2;
	ctypeRNum t112 = t6 * x[8];
	ctypeRNum t116 = POW2(x[8]);
	ctypeRNum t118 = t80 * x[1];
	ctypeRNum t121 = t10 * t2;
	ctypeRNum t154 = x[7] * x[0];
	ctypeRNum t178 = t110 * t6 - x[8] * t12;

	out[0] = -t1 * t3 * t5 * t6 * vec[8] + t1 * (t5 * t10 * t12 - u[2] * t2) * vec[9];
	out[1] = (2 * t19 * t20 + 2 * t23 * x[8] + t26 * t6 + t10 * GRAV - t6 * u[0]) * t32 * t34 - (-2 * t36 * t2 - 2 * x[6] * x[9] + 2 * t19 * t10 + t26 * t43 - t20 * GRAV - t43 * u[0] - t47 * t2) * t32 * vec[9];
	out[2] = 0.0e0;
	out[3] = t1 * (2 * t19 * t43 + 2 * x[9] * x[7] * t2 * x[1] * t10 - t61 + 2 * t60 * t63 + t26 * t10 - t6 * GRAV - t10 * u[0] + t12 * x[1] * t10 * u[2]) * t34 + t1 * (t76 - 2 * t61 * t20 * t10 - 2 * x[7] * t80 * t82 * x[8] + t86 + t87 - t88 + t90) * vec[9];
	out[4] = -(t76 - 2 * t95 * x[8] + t86 + t87 - t88 + t90) * t1 / t80 * vec[8] + t1 * (2 * t36 * t12 - t5 * x[1] * t62 + 2 * t5 * t80 * t63 - 4 * t110 * x[1] * t112 * t12 + x[1] * t116 - 2 * t118 * t116 + t26 * t121 - t6 * t2 * GRAV - t121 * u[0] + t47 * t12) * vec[9];
	out[5] = vec[0] - 2 * t1 * x[7] * t2 * vec[9];
	out[6] = vec[1] - 2 * t1 * (x[7] * t12 * t6 + t2 * x[8]) * t34 - 2 * (x[9] - x[7] * t10) * t1 * vec[9];
	out[7] = vec[2] - 2 * t1 * t6 * (x[6] * t12 + x[9] * t2 * x[1] - t110 * t82 + t154) * t3 * vec[8] + 2 * (-x[5] * t2 + x[6] * t10 + t110 * t63 * t12 + t118 * t112 + t154 * t43) * t1 * vec[9];
	out[8] = vec[3] + 2 * t1 * (-t23 + t95 * t12) * t34 + 2 * t2 * t178 * vec[9];
	out[9] = vec[4] - 2 * t178 * t3 * vec[8] - 2 * t1 * x[6] * vec[9];
}
/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	ctypeRNum t1 = 1 / x[1];
	ctypeRNum t2 = COS(x[4]);
	ctypeRNum t3 = 1 / t2;
	ctypeRNum t5 = COS(x[3]);
	ctypeRNum t6 = t5 * vec[8];
	ctypeRNum t8 = SIN(x[3]);
	ctypeRNum t10 = SIN(x[4]);

	out[0] = vec[5] + t1 * t3 * t6 - t1 * t8 * t10 * vec[9];
	out[1] = vec[6];
	out[2] = vec[7] - t3 * t10 * t6 - (x[0] * t2 - x[1] * t8) * t1 * vec[9];
}
/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}


/** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam) 
    -------------------------------------------------- **/
void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	typeInt i;
	ctypeRNum *Q = (ctypeRNum*)userparam + iQ;
	ctypeRNum *R = (ctypeRNum*)userparam + iR;

	out[0] = 0;
	for (i = 0; i <= 9; i++) {
		out[0] += Q[i] * POW2(x[i] - xdes[i]);
	}
	for (i = 0; i <= 2; i++) {
		out[0] += R[i] * u[i] * u[i];
	}
	out[0] *= 0.5;
}
/** Gradient dl/dx **/
void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	typeInt i;
	ctypeRNum *Q = (ctypeRNum*)userparam + iQ;

	for (i = 0; i <= 9; i++) {
		out[i] = Q[i] * (x[i] - xdes[i]);
	}
}
/** Gradient dl/du **/
void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	typeInt i;
	ctypeRNum *R = (ctypeRNum*)userparam + iR;

	for (i = 0; i <= 2; i++) {
		out[i] = R[i] * u[i];
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
	typeInt i;
	ctypeRNum *P = (ctypeRNum*)userparam + iP;

	out[0] = 0;
	for (i = 0; i <= 9; i++) {
		out[0] += P[i] * POW2(x[i] - xdes[i]);
	}
	out[0] *= 0.5;
}
/** Gradient dV/dx **/
void dVdx(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
	typeInt i;
	ctypeRNum *P = (ctypeRNum*)userparam + iP;

	for (i = 0; i <= 9; i++) {
		out[i] = P[i] * (x[i] - xdes[i]);
	}
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
