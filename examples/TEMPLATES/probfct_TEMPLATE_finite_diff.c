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
 * This probfct file provides an interface to GRAMPC. The most general
 * formulation of the optimal control problem (OCP) that can be solved
 * by GRAMPC has the following structure
 *                                           _T
 *	                                    /
 *      min    J(u,p,T;x0) = V(T,x(T),p) + / l(t,x(t),u(t),p) dt
 *   u(.),p,T                            _/
 *                                      0
 *              .
 *      s.t.   Mx(t) = f(t0+t,x(t),u(t),p), x(0) = x0
 *             g(t,x(t),u(t),p)  = 0,   gT(T,x(T),p)  = 0
 *             h(t,x(t),u(t),p) <= 0,   hT(T,x(t),p) <= 0
 *             u_min <= u(t) <= u_max
 *             p_min <=  p   <= p_max
 *             T_min <=  T   <= T_max
 *
 *
 */

#include "probfct.h"

#include "finite_diff.h"
#define FINITE_DIFF_STEP_SIZE 1e-9

/** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng), 
    inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT, typeUSERPARAM *userparam)
{
}


/** System function f(t,x,u,p,param,userparam) 
    ------------------------------------ **/
void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    finite_diff_vec(out, t, x, u, p, vec, param, userparam,
                    memory, FINITE_DIFF_STEP_SIZE, DX, param->Nx, &ffct);
}
/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    finite_diff_vec(out, t, x, u, p, vec, param, userparam,
                    memory, FINITE_DIFF_STEP_SIZE, DU, param->Nx, &ffct);
}
/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    finite_diff_vec(out, t, x, u, p, vec, param, userparam,
                    memory, FINITE_DIFF_STEP_SIZE, DP, param->Nx, &ffct);
}


/** Integral cost l(t,x(t),u(t),p,param,userparam) 
    -------------------------------------------------- **/
void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Gradient dl/dx **/
void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    typeRNum vec = 1.0;
    finite_diff_vec(out, t, x, u, p, &vec, param, userparam,
                    memory, FINITE_DIFF_STEP_SIZE, DX, 1, &lfct);
}
/** Gradient dl/du **/
void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    typeRNum vec = 1.0;
    finite_diff_vec(out, t, x, u, p, &vec, param, userparam,
                    memory, FINITE_DIFF_STEP_SIZE, DU, 1, &lfct);
}
/** Gradient dl/dp **/
void dldp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    typeRNum vec = 1.0;
    finite_diff_vec(out, t, x, u, p, &vec, param, userparam,
                    memory, FINITE_DIFF_STEP_SIZE, DP, 1, &lfct);
}


/** Terminal cost V(T,x(T),p,param,userparam) 
    ---------------------------------------- **/
void Vfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Gradient dV/dx **/
void dVdx(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    typeRNum vec = 1.0;
    finite_diff_vec(out, T, x, NULL, p, &vec, param, userparam,
                    memory, FINITE_DIFF_STEP_SIZE, DX, 1, &Vfct_wrapper);
}
/** Gradient dV/dp **/
void dVdp(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    typeRNum vec = 1.0;
    finite_diff_vec(out, T, x, NULL, p, &vec, param, userparam,
                    memory, FINITE_DIFF_STEP_SIZE, DP, 1, &Vfct_wrapper);
}
/** Gradient dV/dT **/
void dVdT(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    typeRNum vec = 1.0;
    finite_diff_vec(out, T, x, NULL, p, &vec, param, userparam,
                    memory, FINITE_DIFF_STEP_SIZE, DT, 1, &Vfct_wrapper);
}


/** Equality constraints g(t,x(t),u(t),p,param,userparam) = 0 
    --------------------------------------------------- **/
void gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Jacobian dg/dx multiplied by vector vec, i.e. (dg/dx)^T*vec or vec^T*(dg/dx) **/
void dgdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    finite_diff_vec(out, t, x, u, p, vec, param, userparam,
                    memory, FINITE_DIFF_STEP_SIZE, DX, param->Ng, &gfct);
}
/** Jacobian dg/du multiplied by vector vec, i.e. (dg/du)^T*vec or vec^T*(dg/du) **/
void dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    finite_diff_vec(out, t, x, u, p, vec, param, userparam,
                    memory, FINITE_DIFF_STEP_SIZE, DU, param->Ng, &gfct);
}
/** Jacobian dg/dp multiplied by vector vec, i.e. (dg/dp)^T*vec or vec^T*(dg/dp) **/
void dgdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    finite_diff_vec(out, t, x, u, p, vec, param, userparam,
                    memory, FINITE_DIFF_STEP_SIZE, DP, param->Ng, &gfct);
}


/** Inequality constraints h(t,x(t),u(t),p,param,userparam) <= 0 
    ------------------------------------------------------ **/
void hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
void dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    finite_diff_vec(out, t, x, u, p, vec, param, userparam,
                    memory, FINITE_DIFF_STEP_SIZE, DX, param->Nh, &hfct);
}
/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
void dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    finite_diff_vec(out, t, x, u, p, vec, param, userparam,
                    memory, FINITE_DIFF_STEP_SIZE, DU, param->Nh, &hfct);
}
/** Jacobian dh/dp multiplied by vector vec, i.e. (dh/dp)^T*vec or vec^T*(dg/dp) **/
void dhdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    finite_diff_vec(out, t, x, u, p, vec, param, userparam,
                    memory, FINITE_DIFF_STEP_SIZE, DP, param->Nh, &hfct);
}


/** Terminal equality constraints gT(T,x(T),p,param,userparam) = 0 
    -------------------------------------------------------- **/
void gTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Jacobian dgT/dx multiplied by vector vec, i.e. (dgT/dx)^T*vec or vec^T*(dgT/dx) **/
void dgTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    finite_diff_vec(out, T, x, NULL, p, vec, param, userparam,
                    memory, FINITE_DIFF_STEP_SIZE, DX, param->NgT, &gTfct_wrapper);
}
/** Jacobian dgT/dp multiplied by vector vec, i.e. (dgT/dp)^T*vec or vec^T*(dgT/dp) **/
void dgTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    finite_diff_vec(out, T, x, NULL, p, vec, param, userparam,
                    memory, FINITE_DIFF_STEP_SIZE, DP, param->NgT, &gTfct_wrapper);
}
/** Jacobian dgT/dT multiplied by vector vec, i.e. (dgT/dT)^T*vec or vec^T*(dgT/dT) **/
void dgTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    finite_diff_vec(out, T, x, NULL, p, vec, param, userparam,
                    memory, FINITE_DIFF_STEP_SIZE, DT, param->NgT, &gTfct_wrapper);
}


/** Terminal inequality constraints hT(T,x(T),p,param,userparam) <= 0 
    ----------------------------------------------------------- **/
void hTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
}
/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
void dhTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    finite_diff_vec(out, T, x, NULL, p, vec, param, userparam,
                    memory, FINITE_DIFF_STEP_SIZE, DX, param->NhT, &hTfct_wrapper);
}
/** Jacobian dhT/dp multiplied by vector vec, i.e. (dhT/dp)^T*vec or vec^T*(dhT/dp) **/
void dhTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    finite_diff_vec(out, T, x, NULL, p, vec, param, userparam,
                    memory, FINITE_DIFF_STEP_SIZE, DP, param->NhT, &hTfct_wrapper);
}
/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
void dhTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    finite_diff_vec(out, T, x, NULL, p, vec, param, userparam,
                    memory, FINITE_DIFF_STEP_SIZE, DT, param->NhT, &hTfct_wrapper);
}


/** Additional functions required for semi-implicit systems 
    M*dx/dt(t) = f(t0+t,x(t),u(t),p) using the solver RODAS,
    see GRAMPC docu for more information
    ------------------------------------------------------- **/
/** Jacobian df/dx in vector form (column-wise) **/
void dfdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    finite_diff_dfdx(out, t, x, u, p, param, userparam,
                     memory, FINITE_DIFF_STEP_SIZE, param->Nx, param->Nx, 0);
}
/** Jacobian df/dx in vector form (column-wise) **/
void dfdxtrans(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    // CAUTION: userparam must provide valid working memory of size 3*max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    typeRNum* memory = (typeRNum*)userparam;
    finite_diff_dfdx(out, t, x, u, p, param, userparam,
                     memory, FINITE_DIFF_STEP_SIZE, param->Nx, param->Nx, 1);
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
