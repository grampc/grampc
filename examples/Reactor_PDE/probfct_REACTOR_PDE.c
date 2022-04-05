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
 * This probfct-file describes the PDE reactor problem from
 * Utz, T., Graichen, K., Kugi, A.: Trajectory Planning and Receding Horizon
 * Tracking Control of a Quasilinear Diffusion-Convection-Reaction System.
 * IFAC Proceedings Volumes 43(14), 587-592 (2010)
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

/* square macro */
#define POW2(a) ((a)*(a))

/** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng), 
    inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT, typeUSERPARAM *userparam)
{
	*Nx = 11;
	*Nu = 1;
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
	out[0] = 0.005000000000000002*POW2(x[0]) + x[0] * (-19.46666666666667 + 0.0033333333333333335*x[1]) + (19.516666666666666 + 0.0016666666666666668*x[1])*x[1];
	out[1] = 0.0016666666666666672*POW2(x[0]) + 0.010000000000000004*POW2(x[1]) + x[0] * (20.516666666666666 + 0.0033333333333333335*x[1]) + x[1] * (-39.93333333333334 + 0.0033333333333333335*x[2]) + (19.516666666666666 + 0.0016666666666666668*x[2])*x[2];
	out[2] = 0.0016666666666666672*POW2(x[1]) + 0.010000000000000002*POW2(x[2]) + x[1] * (20.516666666666666 + 0.0033333333333333335*x[2]) + x[2] * (-39.93333333333334 + 0.0033333333333333327*x[3]) + (19.51666666666667 + 0.0016666666666666663*x[3])*x[3];
	out[3] = 0.0016666666666666668*POW2(x[2]) + 0.010000000000000004*POW2(x[3]) + x[2] * (20.51666666666667 + 0.0033333333333333327*x[3]) + x[3] * (-39.93333333333334 + 0.003333333333333335*x[4]) + (19.51666666666666 + 0.001666666666666668*x[4])*x[4];
	out[4] = 0.0016666666666666676*POW2(x[3]) + 0.010000000000000004*POW2(x[4]) + x[3] * (20.51666666666666 + 0.0033333333333333344*x[4]) + x[4] * (-39.93333333333334 + 0.0033333333333333327*x[5]) + (19.51666666666667 + 0.0016666666666666666*x[5])*x[5];
	out[5] = 0.0016666666666666668*POW2(x[4]) + 0.010000000000000002*POW2(x[5]) + x[4] * (20.51666666666667 + 0.0033333333333333327*x[5]) + x[5] * (-39.933333333333344 + 0.0033333333333333327*x[6]) + (19.51666666666667 + 0.0016666666666666663*x[6])*x[6];
	out[6] = 0.0016666666666666668*POW2(x[5]) + 0.010000000000000002*POW2(x[6]) + x[5] * (20.51666666666667 + 0.0033333333333333327*x[6]) + x[6] * (-39.933333333333344 + 0.0033333333333333327*x[7]) + (19.51666666666667 + 0.0016666666666666663*x[7])*x[7];
	out[7] = 0.0016666666666666668*POW2(x[6]) + 0.010000000000000007*POW2(x[7]) + x[6] * (20.51666666666667 + 0.0033333333333333327*x[7]) + x[7] * (-39.93333333333332 + 0.0033333333333333366*x[8]) + (19.516666666666648 + 0.0016666666666666685*x[8])*x[8];
	out[8] = 0.0016666666666666683*POW2(x[7]) + 0.010000000000000007*POW2(x[8]) + x[7] * (20.516666666666648 + 0.0033333333333333366*x[8]) + x[8] * (-39.93333333333332 + 0.0033333333333333327*x[9]) + (19.51666666666667 + 0.0016666666666666663*x[9])*x[9];
	out[9] = 0.0016666666666666668*POW2(x[8]) + 0.010000000000000002*POW2(x[9]) + x[8] * (20.51666666666667 + 0.0033333333333333327*x[9]) + x[9] * (-39.933333333333344 + 0.0033333333333333327*x[10]) + (19.51666666666667 + 0.0016666666666666663*x[10])*x[10];
	out[10] = 0.0016666666666666668*POW2(x[9]) + 2 * u[0] + x[9] * (20.51666666666667 + 0.0033333333333333327*x[10]) + (-20.466666666666672 + 0.005000000000000001*x[10])*x[10];
}
/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	out[0] = vec[1] * (20.516666666666666 + 0.0033333333333333344*x[0] + 0.0033333333333333335*x[1]) + vec[0] * (-19.46666666666667 + 0.010000000000000004*x[0] + 0.0033333333333333335*x[1]);
	out[1] = vec[0] * (19.516666666666666 + 0.0033333333333333335*x[0] + 0.0033333333333333335*x[1]) + vec[2] * (20.516666666666666 + 0.0033333333333333344*x[1] + 0.0033333333333333335*x[2]) + vec[1] * (-39.93333333333334 + 0.0033333333333333335*x[0] + 0.020000000000000007*x[1] + 0.0033333333333333335*x[2]);
	out[2] = vec[1] * (19.516666666666666 + 0.0033333333333333335*x[1] + 0.0033333333333333335*x[2]) + vec[3] * (20.51666666666667 + 0.0033333333333333335*x[2] + 0.0033333333333333327*x[3]) + vec[2] * (-39.93333333333334 + 0.0033333333333333335*x[1] + 0.020000000000000004*x[2] + 0.0033333333333333327*x[3]);
	out[3] = vec[2] * (19.51666666666667 + 0.0033333333333333327*x[2] + 0.0033333333333333327*x[3]) + vec[4] * (20.51666666666666 + 0.0033333333333333353*x[3] + 0.0033333333333333344*x[4]) + vec[3] * (-39.93333333333334 + 0.0033333333333333327*x[2] + 0.020000000000000007*x[3] + 0.003333333333333335*x[4]);
	out[4] = vec[3] * (19.51666666666666 + 0.003333333333333335*x[3] + 0.003333333333333336*x[4]) + vec[5] * (20.51666666666667 + 0.0033333333333333335*x[4] + 0.0033333333333333327*x[5]) + vec[4] * (-39.93333333333334 + 0.0033333333333333344*x[3] + 0.020000000000000007*x[4] + 0.0033333333333333327*x[5]);
	out[5] = vec[4] * (19.51666666666667 + 0.0033333333333333327*x[4] + 0.003333333333333333*x[5]) + vec[6] * (20.51666666666667 + 0.0033333333333333335*x[5] + 0.0033333333333333327*x[6]) + vec[5] * (-39.933333333333344 + 0.0033333333333333327*x[4] + 0.020000000000000004*x[5] + 0.0033333333333333327*x[6]);
	out[6] = vec[5] * (19.51666666666667 + 0.0033333333333333327*x[5] + 0.0033333333333333327*x[6]) + vec[7] * (20.51666666666667 + 0.0033333333333333335*x[6] + 0.0033333333333333327*x[7]) + vec[6] * (-39.933333333333344 + 0.0033333333333333327*x[5] + 0.020000000000000004*x[6] + 0.0033333333333333327*x[7]);
	out[7] = vec[6] * (19.51666666666667 + 0.0033333333333333327*x[6] + 0.0033333333333333327*x[7]) + vec[8] * (20.516666666666648 + 0.0033333333333333366*x[7] + 0.0033333333333333366*x[8]) + vec[7] * (-39.93333333333332 + 0.0033333333333333327*x[6] + 0.020000000000000014*x[7] + 0.0033333333333333366*x[8]);
	out[8] = vec[7] * (19.516666666666648 + 0.0033333333333333366*x[7] + 0.003333333333333337*x[8]) + vec[9] * (20.51666666666667 + 0.0033333333333333335*x[8] + 0.0033333333333333327*x[9]) + vec[8] * (-39.93333333333332 + 0.0033333333333333366*x[7] + 0.020000000000000014*x[8] + 0.0033333333333333327*x[9]);
	out[9] = vec[8] * (19.51666666666667 + 0.0033333333333333327*x[8] + 0.0033333333333333327*x[9]) + vec[10] * (20.51666666666667 + 0.0033333333333333335*x[9] + 0.0033333333333333327*x[10]) + vec[9] * (-39.933333333333344 + 0.0033333333333333327*x[8] + 0.020000000000000004*x[9] + 0.0033333333333333327*x[10]);
	out[10] = vec[9] * (19.51666666666667 + 0.0033333333333333327*x[9] + 0.0033333333333333327*x[10]) + vec[10] * (-20.466666666666672 + 0.0033333333333333327*x[9] + 0.010000000000000002*x[10]);
}
/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	out[0] = 2 * vec[10];
}
/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}


/** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam) 
    -------------------------------------------------- **/
void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	typeRNum* param = (typeRNum*)userparam;
	out[0] = 0.5*param[22] * POW2(u[0] - udes[0]) + 0.5*(param[11] * POW2(x[0] - xdes[0]) + param[12] * POW2(x[1] - xdes[1]) + param[13] * POW2(x[2] - xdes[2]) + param[14] * POW2(x[3] - xdes[3]) + param[15] * POW2(x[4] - xdes[4]) + param[16] * POW2(x[5] - xdes[5]) + param[17] * POW2(x[6] - xdes[6]) + param[18] * POW2(x[7] - xdes[7]) + param[19] * POW2(x[8] - xdes[8]) + param[20] * POW2(x[9] - xdes[9]) + param[21] * POW2(x[10] - xdes[10]));
}
/** Gradient dl/dx **/
void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	typeRNum* param = (typeRNum*)userparam;
	out[0] = 1 * param[11] * (x[0] - xdes[0]);
	out[1] = 1 * param[12] * (x[1] - xdes[1]);
	out[2] = 1 * param[13] * (x[2] - xdes[2]);
	out[3] = 1 * param[14] * (x[3] - xdes[3]);
	out[4] = 1 * param[15] * (x[4] - xdes[4]);
	out[5] = 1 * param[16] * (x[5] - xdes[5]);
	out[6] = 1 * param[17] * (x[6] - xdes[6]);
	out[7] = 1 * param[18] * (x[7] - xdes[7]);
	out[8] = 1 * param[19] * (x[8] - xdes[8]);
	out[9] = 1 * param[20] * (x[9] - xdes[9]);
	out[10] = 1 * param[21] * (x[10] - xdes[10]);
}
/** Gradient dl/du **/
void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	typeRNum* param = (typeRNum*)userparam;
	out[0] = 0 + 1 * param[22] * (u[0] - udes[0]);
}
/** Gradient dl/dp **/
void dldp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
}


/** Terminal cost V(T,x(T),p,xdes,userparam) 
    ---------------------------------------- **/
void Vfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
	typeRNum* param = (typeRNum*)userparam;
	out[0] = 0.5*(param[0] * POW2(x[0] - xdes[0]) + param[1] * POW2(x[1] - xdes[1]) + param[2] * POW2(x[2] - xdes[2]) + param[3] * POW2(x[3] - xdes[3]) + param[4] * POW2(x[4] - xdes[4]) + param[5] * POW2(x[5] - xdes[5]) + param[6] * POW2(x[6] - xdes[6]) + param[7] * POW2(x[7] - xdes[7]) + param[8] * POW2(x[8] - xdes[8]) + param[9] * POW2(x[9] - xdes[9]) + param[10] * POW2(x[10] - xdes[10]));
}

/** Gradient dV/dx **/
void dVdx(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
	typeRNum* param = (typeRNum*)userparam;
	out[0] = 34.6410161516297*param[0] * (x[0] - xdes[0]) - 9.282032303259404*param[1] * (x[1] - xdes[1]) + 2.4871130614079147*param[2] * (x[2] - xdes[2]) - 0.666419942372254*param[3] * (x[3] - xdes[3]) + 0.17856670808110178*param[4] * (x[4] - xdes[4]) - 0.0478468899521531*param[5] * (x[5] - xdes[5]) + 0.012820851727510638*param[6] * (x[6] - xdes[6]) - 0.003436516957889448*param[7] * (x[7] - xdes[7]) + 0.0009252161040471597*param[8] * (x[8] - xdes[8]) - 0.00026434745829918854*param[9] * (x[9] - xdes[9]) + 0.00013217372914959427*param[10] * (x[10] - xdes[10]);
	out[1] = -9.282032303259404*param[0] * (x[0] - xdes[0]) + 18.56406460651881*param[1] * (x[1] - xdes[1]) - 4.974226122815829*param[2] * (x[2] - xdes[2]) + 1.332839884744508*param[3] * (x[3] - xdes[3]) - 0.35713341616220357*param[4] * (x[4] - xdes[4]) + 0.0956937799043062*param[5] * (x[5] - xdes[5]) - 0.025641703455021276*param[6] * (x[6] - xdes[6]) + 0.006873033915778896*param[7] * (x[7] - xdes[7]) - 0.0018504322080943195*param[8] * (x[8] - xdes[8]) + 0.0005286949165983771*param[9] * (x[9] - xdes[9]) - 0.00026434745829918854*param[10] * (x[10] - xdes[10]);
	out[2] = 2.4871130614079147*param[0] * (x[0] - xdes[0]) - 4.974226122815829*param[1] * (x[1] - xdes[1]) + 17.409791429855403*param[2] * (x[2] - xdes[2]) - 4.6649395966057785*param[3] * (x[3] - xdes[3]) + 1.2499669565677127*param[4] * (x[4] - xdes[4]) - 0.3349282296650718*param[5] * (x[5] - xdes[5]) + 0.08974596209257449*param[6] * (x[6] - xdes[6]) - 0.024055618705226146*param[7] * (x[7] - xdes[7]) + 0.00647651272833012*param[8] * (x[8] - xdes[8]) - 0.0018504322080943204*param[9] * (x[9] - xdes[9]) + 0.0009252161040471603*param[10] * (x[10] - xdes[10]);
	out[3] = -0.666419942372254*param[0] * (x[0] - xdes[0]) + 1.332839884744508*param[1] * (x[1] - xdes[1]) - 4.6649395966057785*param[2] * (x[2] - xdes[2]) + 17.326918501678605*param[3] * (x[3] - xdes[3]) - 4.6427344101086465*param[4] * (x[4] - xdes[4]) + 1.2440191387559807*param[5] * (x[5] - xdes[5]) - 0.33334214491527653*param[6] * (x[6] - xdes[6]) + 0.08934944090512564*param[7] * (x[7] - xdes[7]) - 0.02405561870522615*param[8] * (x[8] - xdes[8]) + 0.0068730339157789*param[9] * (x[9] - xdes[9]) - 0.0034365169578894506*param[10] * (x[10] - xdes[10]);
	out[4] = 0.1785667080811018*param[0] * (x[0] - xdes[0]) - 0.3571334161622036*param[1] * (x[1] - xdes[1]) + 1.2499669565677127*param[2] * (x[2] - xdes[2]) - 4.642734410108647*param[3] * (x[3] - xdes[3]) + 17.320970683866875*param[4] * (x[4] - xdes[4]) - 4.641148325358851*param[5] * (x[5] - xdes[5]) + 1.243622617568532*param[6] * (x[6] - xdes[6]) - 0.33334214491527653*param[7] * (x[7] - xdes[7]) + 0.08974596209257452*param[8] * (x[8] - xdes[8]) - 0.025641703455021293*param[9] * (x[9] - xdes[9]) + 0.012820851727510648*param[10] * (x[10] - xdes[10]);
	out[5] = -0.0478468899521531*param[0] * (x[0] - xdes[0]) + 0.0956937799043062*param[1] * (x[1] - xdes[1]) - 0.33492822966507174*param[2] * (x[2] - xdes[2]) + 1.244019138755981*param[3] * (x[3] - xdes[3]) - 4.641148325358851*param[4] * (x[4] - xdes[4]) + 17.320574162679428*param[5] * (x[5] - xdes[5]) - 4.641148325358852*param[6] * (x[6] - xdes[6]) + 1.2440191387559805*param[7] * (x[7] - xdes[7]) - 0.3349282296650719*param[8] * (x[8] - xdes[8]) + 0.09569377990430628*param[9] * (x[9] - xdes[9]) - 0.04784688995215314*param[10] * (x[10] - xdes[10]);
	out[6] = 0.012820851727510638*param[0] * (x[0] - xdes[0]) - 0.025641703455021276*param[1] * (x[1] - xdes[1]) + 0.08974596209257447*param[2] * (x[2] - xdes[2]) - 0.33334214491527664*param[3] * (x[3] - xdes[3]) + 1.2436226175685319*param[4] * (x[4] - xdes[4]) - 4.641148325358852*param[5] * (x[5] - xdes[5]) + 17.32097068386688*param[6] * (x[6] - xdes[6]) - 4.642734410108646*param[7] * (x[7] - xdes[7]) + 1.249966956567713*param[8] * (x[8] - xdes[8]) - 0.35713341616220373*param[9] * (x[9] - xdes[9]) + 0.1785667080811019*param[10] * (x[10] - xdes[10]);
	out[7] = -0.003436516957889448*param[0] * (x[0] - xdes[0]) + 0.006873033915778896*param[1] * (x[1] - xdes[1]) - 0.02405561870522614*param[2] * (x[2] - xdes[2]) + 0.08934944090512566*param[3] * (x[3] - xdes[3]) - 0.3333421449152764*param[4] * (x[4] - xdes[4]) + 1.2440191387559805*param[5] * (x[5] - xdes[5]) - 4.642734410108646*param[6] * (x[6] - xdes[6]) + 17.326918501678602*param[7] * (x[7] - xdes[7]) - 4.664939596605781*param[8] * (x[8] - xdes[8]) + 1.332839884744509*param[9] * (x[9] - xdes[9]) - 0.6664199423722547*param[10] * (x[10] - xdes[10]);
	out[8] = 0.0009252161040471595*param[0] * (x[0] - xdes[0]) - 0.001850432208094319*param[1] * (x[1] - xdes[1]) + 0.006476512728330117*param[2] * (x[2] - xdes[2]) - 0.02405561870522615*param[3] * (x[3] - xdes[3]) + 0.08974596209257446*param[4] * (x[4] - xdes[4]) - 0.3349282296650718*param[5] * (x[5] - xdes[5]) + 1.249966956567713*param[6] * (x[6] - xdes[6]) - 4.66493959660578*param[7] * (x[7] - xdes[7]) + 17.4097914298554*param[8] * (x[8] - xdes[8]) - 4.974226122815829*param[9] * (x[9] - xdes[9]) + 2.487113061407915*param[10] * (x[10] - xdes[10]);
	out[9] = -0.00026434745829918854*param[0] * (x[0] - xdes[0]) + 0.0005286949165983771*param[1] * (x[1] - xdes[1]) - 0.0018504322080943197*param[2] * (x[2] - xdes[2]) + 0.006873033915778902*param[3] * (x[3] - xdes[3]) - 0.025641703455021283*param[4] * (x[4] - xdes[4]) + 0.09569377990430626*param[5] * (x[5] - xdes[5]) - 0.3571334161622038*param[6] * (x[6] - xdes[6]) + 1.3328398847445089*param[7] * (x[7] - xdes[7]) - 4.974226122815828*param[8] * (x[8] - xdes[8]) + 18.564064606518812*param[9] * (x[9] - xdes[9]) - 9.282032303259408*param[10] * (x[10] - xdes[10]);
	out[10] = 0.0001321737291495943*param[0] * (x[0] - xdes[0]) - 0.0002643474582991886*param[1] * (x[1] - xdes[1]) + 0.0009252161040471601*param[2] * (x[2] - xdes[2]) - 0.003436516957889452*param[3] * (x[3] - xdes[3]) + 0.012820851727510645*param[4] * (x[4] - xdes[4]) - 0.047846889952153145*param[5] * (x[5] - xdes[5]) + 0.17856670808110195*param[6] * (x[6] - xdes[6]) - 0.6664199423722547*param[7] * (x[7] - xdes[7]) + 2.487113061407915*param[8] * (x[8] - xdes[8]) - 9.282032303259408*param[9] * (x[9] - xdes[9]) + 34.641016151629714*param[10] * (x[10] - xdes[10]);
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
	out[0] = 0;
	out[1] = -19.46666666666667 + 0.010000000000000004*x[0] + 0.0033333333333333335*x[1];
	out[2] = 20.516666666666666 + 0.0033333333333333344*x[0] + 0.0033333333333333335*x[1];
	out[3] = 19.516666666666666 + 0.0033333333333333335*x[0] + 0.0033333333333333335*x[1];
	out[4] = -39.93333333333334 + 0.0033333333333333335*x[0] + 0.020000000000000007*x[1] + 0.0033333333333333335*x[2];
	out[5] = 20.516666666666666 + 0.0033333333333333344*x[1] + 0.0033333333333333335*x[2];
	out[6] = 19.516666666666666 + 0.0033333333333333335*x[1] + 0.0033333333333333335*x[2];
	out[7] = -39.93333333333334 + 0.0033333333333333335*x[1] + 0.020000000000000004*x[2] + 0.0033333333333333327*x[3];
	out[8] = 20.51666666666667 + 0.0033333333333333335*x[2] + 0.0033333333333333327*x[3];
	out[9] = 19.51666666666667 + 0.0033333333333333327*x[2] + 0.0033333333333333327*x[3];
	out[10] = -39.93333333333334 + 0.0033333333333333327*x[2] + 0.020000000000000007*x[3] + 0.003333333333333335*x[4];
	out[11] = 20.51666666666666 + 0.0033333333333333353*x[3] + 0.0033333333333333344*x[4];
	out[12] = 19.51666666666666 + 0.003333333333333335*x[3] + 0.003333333333333336*x[4];
	out[13] = -39.93333333333334 + 0.0033333333333333344*x[3] + 0.020000000000000007*x[4] + 0.0033333333333333327*x[5];
	out[14] = 20.51666666666667 + 0.0033333333333333335*x[4] + 0.0033333333333333327*x[5];
	out[15] = 19.51666666666667 + 0.0033333333333333327*x[4] + 0.003333333333333333*x[5];
	out[16] = -39.933333333333344 + 0.0033333333333333327*x[4] + 0.020000000000000004*x[5] + 0.0033333333333333327*x[6];
	out[17] = 20.51666666666667 + 0.0033333333333333335*x[5] + 0.0033333333333333327*x[6];
	out[18] = 19.51666666666667 + 0.0033333333333333327*x[5] + 0.0033333333333333327*x[6];
	out[19] = -39.933333333333344 + 0.0033333333333333327*x[5] + 0.020000000000000004*x[6] + 0.0033333333333333327*x[7];
	out[20] = 20.51666666666667 + 0.0033333333333333335*x[6] + 0.0033333333333333327*x[7];
	out[21] = 19.51666666666667 + 0.0033333333333333327*x[6] + 0.0033333333333333327*x[7];
	out[22] = -39.93333333333332 + 0.0033333333333333327*x[6] + 0.020000000000000014*x[7] + 0.0033333333333333366*x[8];
	out[23] = 20.516666666666648 + 0.0033333333333333366*x[7] + 0.0033333333333333366*x[8];
	out[24] = 19.516666666666648 + 0.0033333333333333366*x[7] + 0.003333333333333337*x[8];
	out[25] = -39.93333333333332 + 0.0033333333333333366*x[7] + 0.020000000000000014*x[8] + 0.0033333333333333327*x[9];
	out[26] = 20.51666666666667 + 0.0033333333333333335*x[8] + 0.0033333333333333327*x[9];
	out[27] = 19.51666666666667 + 0.0033333333333333327*x[8] + 0.0033333333333333327*x[9];
	out[28] = -39.933333333333344 + 0.0033333333333333327*x[8] + 0.020000000000000004*x[9] + 0.0033333333333333327*x[10];
	out[29] = 20.51666666666667 + 0.0033333333333333335*x[9] + 0.0033333333333333327*x[10];
	out[30] = 19.51666666666667 + 0.0033333333333333327*x[9] + 0.0033333333333333327*x[10];
	out[31] = -20.466666666666672 + 0.0033333333333333327*x[9] + 0.010000000000000002*x[10];
	out[32] = 0;
}
/** Jacobian df/dx in vector form (column-wise) **/
void dfdxtrans(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	out[0] = 0;
	out[1] = -19.46666666666667 + 0.010000000000000004*x[0] + 0.0033333333333333335*x[1];
	out[2] = 19.516666666666666 + 0.0033333333333333335*x[0] + 0.0033333333333333335*x[1];
	out[3] = 20.516666666666666 + 0.0033333333333333344*x[0] + 0.0033333333333333335*x[1];
	out[4] = -39.93333333333334 + 0.0033333333333333335*x[0] + 0.020000000000000007*x[1] + 0.0033333333333333335*x[2];
	out[5] = 19.516666666666666 + 0.0033333333333333335*x[1] + 0.0033333333333333335*x[2];
	out[6] = 20.516666666666666 + 0.0033333333333333344*x[1] + 0.0033333333333333335*x[2];
	out[7] = -39.93333333333334 + 0.0033333333333333335*x[1] + 0.020000000000000004*x[2] + 0.0033333333333333327*x[3];
	out[8] = 19.51666666666667 + 0.0033333333333333327*x[2] + 0.0033333333333333327*x[3];
	out[9] = 20.51666666666667 + 0.0033333333333333335*x[2] + 0.0033333333333333327*x[3];
	out[10] = -39.93333333333334 + 0.0033333333333333327*x[2] + 0.020000000000000007*x[3] + 0.003333333333333335*x[4];
	out[11] = 19.51666666666666 + 0.003333333333333335*x[3] + 0.003333333333333336*x[4];
	out[12] = 20.51666666666666 + 0.0033333333333333353*x[3] + 0.0033333333333333344*x[4];
	out[13] = -39.93333333333334 + 0.0033333333333333344*x[3] + 0.020000000000000007*x[4] + 0.0033333333333333327*x[5];
	out[14] = 19.51666666666667 + 0.0033333333333333327*x[4] + 0.003333333333333333*x[5];
	out[15] = 20.51666666666667 + 0.0033333333333333335*x[4] + 0.0033333333333333327*x[5];
	out[16] = -39.933333333333344 + 0.0033333333333333327*x[4] + 0.020000000000000004*x[5] + 0.0033333333333333327*x[6];
	out[17] = 19.51666666666667 + 0.0033333333333333327*x[5] + 0.0033333333333333327*x[6];
	out[18] = 20.51666666666667 + 0.0033333333333333335*x[5] + 0.0033333333333333327*x[6];
	out[19] = -39.933333333333344 + 0.0033333333333333327*x[5] + 0.020000000000000004*x[6] + 0.0033333333333333327*x[7];
	out[20] = 19.51666666666667 + 0.0033333333333333327*x[6] + 0.0033333333333333327*x[7];
	out[21] = 20.51666666666667 + 0.0033333333333333335*x[6] + 0.0033333333333333327*x[7];
	out[22] = -39.93333333333332 + 0.0033333333333333327*x[6] + 0.020000000000000014*x[7] + 0.0033333333333333366*x[8];
	out[23] = 19.516666666666648 + 0.0033333333333333366*x[7] + 0.003333333333333337*x[8];
	out[24] = 20.516666666666648 + 0.0033333333333333366*x[7] + 0.0033333333333333366*x[8];
	out[25] = -39.93333333333332 + 0.0033333333333333366*x[7] + 0.020000000000000014*x[8] + 0.0033333333333333327*x[9];
	out[26] = 19.51666666666667 + 0.0033333333333333327*x[8] + 0.0033333333333333327*x[9];
	out[27] = 20.51666666666667 + 0.0033333333333333335*x[8] + 0.0033333333333333327*x[9];
	out[28] = -39.933333333333344 + 0.0033333333333333327*x[8] + 0.020000000000000004*x[9] + 0.0033333333333333327*x[10];
	out[29] = 19.51666666666667 + 0.0033333333333333327*x[9] + 0.0033333333333333327*x[10];
	out[30] = 20.51666666666667 + 0.0033333333333333335*x[9] + 0.0033333333333333327*x[10];
	out[31] = -20.466666666666672 + 0.0033333333333333327*x[9] + 0.010000000000000002*x[10];
	out[32] = 0;
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
	out[0] = 0;
	out[1] = 0.03333333333333333;
	out[2] = 0.016666666666666666;
	out[3] = 0.016666666666666666;
	out[4] = 0.06666666666666667;
	out[5] = 0.016666666666666666;
	out[6] = 0.016666666666666666;
	out[7] = 0.06666666666666665;
	out[8] = 0.016666666666666663;
	out[9] = 0.016666666666666663;
	out[10] = 0.06666666666666667;
	out[11] = 0.01666666666666667;
	out[12] = 0.01666666666666667;
	out[13] = 0.06666666666666667;
	out[14] = 0.016666666666666663;
	out[15] = 0.016666666666666663;
	out[16] = 0.06666666666666665;
	out[17] = 0.016666666666666663;
	out[18] = 0.016666666666666663;
	out[19] = 0.06666666666666665;
	out[20] = 0.016666666666666663;
	out[21] = 0.016666666666666663;
	out[22] = 0.06666666666666668;
	out[23] = 0.01666666666666668;
	out[24] = 0.01666666666666668;
	out[25] = 0.06666666666666668;
	out[26] = 0.016666666666666663;
	out[27] = 0.016666666666666663;
	out[28] = 0.06666666666666665;
	out[29] = 0.016666666666666663;
	out[30] = 0.016666666666666663;
	out[31] = 0.033333333333333326;
	out[32] = 0;
}
/** Transposed mass matrix in vector form (column-wise, either banded or full matrix) **/
void Mtrans(typeRNum *out, typeUSERPARAM *userparam)
{
	out[0] = 0;
	out[1] = 0.03333333333333333;
	out[2] = 0.016666666666666666;
	out[3] = 0.016666666666666666;
	out[4] = 0.06666666666666667;
	out[5] = 0.016666666666666666;
	out[6] = 0.016666666666666666;
	out[7] = 0.06666666666666665;
	out[8] = 0.016666666666666663;
	out[9] = 0.016666666666666663;
	out[10] = 0.06666666666666667;
	out[11] = 0.01666666666666667;
	out[12] = 0.01666666666666667;
	out[13] = 0.06666666666666667;
	out[14] = 0.016666666666666663;
	out[15] = 0.016666666666666663;
	out[16] = 0.06666666666666665;
	out[17] = 0.016666666666666663;
	out[18] = 0.016666666666666663;
	out[19] = 0.06666666666666665;
	out[20] = 0.016666666666666663;
	out[21] = 0.016666666666666663;
	out[22] = 0.06666666666666668;
	out[23] = 0.01666666666666668;
	out[24] = 0.01666666666666668;
	out[25] = 0.06666666666666668;
	out[26] = 0.016666666666666663;
	out[27] = 0.016666666666666663;
	out[28] = 0.06666666666666665;
	out[29] = 0.016666666666666663;
	out[30] = 0.016666666666666663;
	out[31] = 0.033333333333333326;
	out[32] = 0;
}
