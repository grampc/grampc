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
 * This probfct-file describes the CSTR reactor problem in the context of
 * moving horizon estimation from
 * Rothfuss, R., Rudolph, J., Zeitz, M.: Flatness Based Control of a
 * Nonlinear Chemical Reactor Model. Automatica 32(10), 1433-1439 (1996)
 *
 *                              _T
 *							   /   ^             
 *      min    J(u,p;y) =     /  ||y(t) - y(t)||^2  dt
 *       p                  _/
 *                          0
 *             .
 *      s.t.   x(t) = f(t0-T+t,x(t)+p,u(t)),  x(0) = 0
 *             y(t) = w(x(t)+p)
 *
 */

#include "probfct.h"

#if USE_typeRNum == USE_FLOAT
#define EXP(a)		expf(a)
#else
#define EXP(a)		exp(a)
#endif

 /* square macro */
#define POW2(a) ((a)*(a))


/** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng), 
    inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT, typeUSERPARAM *userparam)
{
	*Nx = 4;
	*Nu = 2;
	*Np = 4;
	*Ng = 0;
	*Nh = 0;
	*NgT = 0;
	*NhT = 0;
}


/** System function f(t,x,u,p,param,userparam) 
    ------------------------------------ **/
void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
	typeRNum* pSys = (typeRNum*)userparam;

	out[0] = -(EXP(-(pSys[3] / (((typeRNum)273.15) + p[2] + x[2])))*pSys[1] * POW2(p[0] + x[0])) + u[0] *
		(-p[0] + pSys[9] - x[0]) - EXP(-(pSys[2] / (((typeRNum)273.15) + p[2] + x[2])))*pSys[0] * (p[0] + x[0]);
	out[1] = EXP(-(pSys[2] / (((typeRNum)273.15) + p[2] + x[2])))*pSys[0] * (p[0] + x[0]) - EXP(-(pSys[2] /
		(((typeRNum)273.15) + p[2] + x[2])))*pSys[0] * (p[1] + x[1]) - u[0] * (p[1] + x[1]);
	out[2] = -(pSys[7] * (EXP(-(pSys[3] / (((typeRNum)273.15) + p[2] + x[2])))*pSys[1] * pSys[12] *
		POW2(p[0] + x[0]) + EXP(-(pSys[2] / (((typeRNum)273.15) + p[2] + x[2])))*pSys[0] * pSys[10] *
		(p[0] + x[0]) + EXP(-(pSys[2] / (((typeRNum)273.15) + p[2] + x[2])))*pSys[0] * pSys[11] * (p[1] + x[1])))
		+ u[0] * (-p[2] + pSys[8] - x[2]) + pSys[4] * (-p[2] + p[3] - x[2] + x[3]);
	out[3] = pSys[6] * u[1] + pSys[5] * (p[2] - p[3] + x[2] - x[3]);
}
/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
	typeRNum* pSys = (typeRNum*)userparam;

	out[0] = -(EXP(-(pSys[2] / (((typeRNum)273.15) + p[2] + x[2])))*pSys[0] * (vec[0] - vec[1] + vec[2] *
		pSys[7] * pSys[10])) - vec[0] * u[0] - 2 * EXP(-(pSys[3] / (((typeRNum)273.15) + p[2] + x[2])))*
		pSys[1] * (vec[0] + vec[2] * pSys[7] * pSys[12])*(p[0] + x[0]);
	out[1] = -(EXP(-(pSys[2] / (((typeRNum)273.15) + p[2] + x[2])))*pSys[0] * (vec[1] + vec[2] * pSys[7] * pSys[11])) - vec[1] * u[0];
	out[2] = vec[3] * pSys[5] + vec[0] / POW2(((typeRNum)273.15) + p[2] + x[2])*(p[0] + x[0])*(-(EXP(-(pSys[2] / (((typeRNum)273.15) +
		p[2] + x[2])))*pSys[0] * pSys[2]) - EXP(-(pSys[3] / (((typeRNum)273.15) + p[2] + x[2])))*pSys[1] * pSys[3] * (p[0] + x[0])) +
		vec[1] * EXP(-(pSys[2] / (((typeRNum)273.15) + p[2] + x[2])))* pSys[0] * pSys[2] / POW2(273.15 + p[2] + x[2])*(p[0] - p[1] + x[0] -
			x[1]) + vec[2] * (-pSys[4] - u[0] - EXP(-((pSys[2] + pSys[3]) / (((typeRNum)273.15) + p[2] + x[2])))*pSys[7] / POW2(273.15 + p[2] +
				x[2])*(EXP(pSys[2] / (((typeRNum)273.15) + p[2] + x[2]))*pSys[1] * pSys[3] * pSys[12] * (POW2(p[0]) + POW2(x[0]) +
					2 * p[0] * x[0]) + EXP(pSys[3] / (((typeRNum)273.15) + p[2] + x[2]))*pSys[0] * pSys[2] * (pSys[10] * (p[0] + x[0]) +
						pSys[11] * (p[1] + x[1]))));
	out[3] = vec[2] * pSys[4] - vec[3] * pSys[5];
}
/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
	typeRNum* pSys = (typeRNum*)userparam;

	out[0] = vec[2] * pSys[8] + vec[0] * pSys[9] - vec[0] * (p[0] + x[0]) - vec[1] * (p[1] + x[1]) - \
		vec[2] * (p[2] + x[2]);
	out[1] = vec[3] * pSys[6];
}
/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
	typeRNum* pSys = (typeRNum*)userparam;

	out[0] = -(EXP(-(pSys[2] / (((typeRNum)273.15) + p[2] + x[2])))*pSys[0] * (vec[0] - vec[1] + vec[2] *
		pSys[7] * pSys[10])) - vec[0] * u[0] - 2 * EXP(-(pSys[3] / (((typeRNum)273.15) + p[2] + x[2])))*pSys[1] *
		(vec[0] + vec[2] * pSys[7] * pSys[12])*(p[0] + x[0]);
	out[1] = -(EXP(-(pSys[2] / (((typeRNum)273.15) + p[2] + x[2])))*pSys[0] * (vec[1] + vec[2] * pSys[7] * pSys[11])) - vec[1] * u[0];
	out[2] = vec[3] * pSys[5] + vec[0] / POW2(273.15 + p[2] + x[2])*(p[0] + x[0])*(-(EXP(-(pSys[2] /
		(((typeRNum)273.15) + p[2] + x[2])))*pSys[0] * pSys[2]) - EXP(-(pSys[3] / (((typeRNum)273.15) +
			p[2] + x[2])))*pSys[1] * pSys[3] * (p[0] + x[0])) + vec[1] * EXP(-(pSys[2] / (((typeRNum)273.15) + p[2] + x[2])))
		*pSys[0] * pSys[2] / POW2(273.15 + p[2] + x[2])*(p[0] - p[1] + x[0] - x[1]) + vec[2] * (-pSys[4] - u[0] -
			EXP(-((pSys[2] + pSys[3]) / (((typeRNum)273.15) + p[2] + x[2])))*pSys[7] / POW2(273.15 + p[2] + x[2])*
			(EXP(pSys[2] / (((typeRNum)273.15) + p[2] + x[2]))*pSys[1] * pSys[3] * pSys[12] * (POW2(p[0]) + POW2(x[0]) +
				2 * p[0] * x[0]) + EXP(pSys[3] / (((typeRNum)273.15) + p[2] + x[2]))*pSys[0] * pSys[2] * (pSys[10] * (p[0] + x[0]) +
					pSys[11] * (p[1] + x[1]))));
	out[3] = vec[2] * pSys[4] - vec[3] * pSys[5];
}


/** Integral cost l(t,x(t),u(t),p,param,userparam) 
    -------------------------------------------------- **/
void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
	typeRNum* pSys = (typeRNum*)userparam;
	typeRNum* pCost = &pSys[14];
	typeRNum* pMeas = &pSys[16];
	typeInt index = (int)floor(t / 2.777777777777778e-04 + 0.00001);
	typeRNum meas1 = pMeas[2 * index];
	typeRNum meas2 = pMeas[1 + 2 * index];

	out[0] = pCost[0] * POW2((p[2] + x[2] - meas1)) +
			 pCost[1] * POW2((p[3] + x[3] - meas2));
}
/** Gradient dl/dx **/
void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
	typeRNum* pSys = (typeRNum*)userparam;
	typeRNum* pCost = &pSys[14];
	typeRNum* pMeas = &pSys[16];
	typeInt index = (int)floor(t / 2.777777777777778e-04 + 0.00001);
	typeRNum meas1 = pMeas[2 * index];
	typeRNum meas2 = pMeas[1 + 2 * index];

	out[0] = 0.0;
	out[1] = 0.0;
	out[2] = 2.0 * pCost[1] * (p[2] + x[2] - meas1);
	out[3] = 2.0 * pCost[0] * (p[3] + x[3] - meas2);
}


/** Gradient dl/du **/
void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Gradient dl/dp **/
void dldp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
	typeRNum* pSys = (typeRNum*)userparam;
	typeRNum* pCost = &pSys[14];
	typeRNum* pMeas = &pSys[16];
	typeInt index = (int)floor(t / 2.777777777777778e-04 + 0.00001);
	typeRNum meas1 = pMeas[2 * index];
	typeRNum meas2 = pMeas[1 + 2 * index];

	out[0] = 0.0;
	out[1] = 0.0;
	out[2] = 2.0 * pCost[0] * (p[2] + x[2] - meas1);
	out[3] = 2.0 * pCost[1] * (p[3] + x[3] - meas2);
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
}
/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
void dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
void dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
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
