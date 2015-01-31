/*
 *
 * This file is part of GRAMPC.
 *
 * GRAMPC - a gradient-based MPC software for real-time applications
 *
 * Copyright (C) 2014 by Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Developed at the Institute of Measurement, Control, and
 * Microtechnology, University of Ulm. All rights reserved.
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
 * License along with GRAMPC. If not, see <http://www.gnu.org/licenses/>.
 *
 */


/*
 * 
 * File: probfct.c
 * Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Date: February 2014
 * Version: v1.0
 *
 * Problem function for grampc toolbox.
 *
 *
 * This probfct file describes the dynamics, the cost function and the corresponding
 * derivatives of an overhead crane with 6 states and 2 controls. For a more
 * detailed model, for instance, see
 *
 * B. Käpernick and K. Graichen, "Model Predictive Control of an Overhead Crane
 * Using Constraint Substitution", Proceedings of the American Control Conference,
 * Washington D.C. (USA), 2013.
 *
 *
 * This probfct file provides an interface to GRAMPC. The underlying
 * optimal control problem (OCP) of the model predictive control (MPC) formulation
 * has the following structure
 *                                  _T
 *                                 /
 *      min    J(u,xk) = V(T,x(T)) + / L(t,x(t),u(t)) dt
 *      u(.)                    _/
 *                             0
 *             .
 *      s.t.   x(t) = f(tk+t,x(t),u(t)), x(0) = xk
 *
 *             Ul <= u(t) <= Uu, t in [0,T]
 *
 * with states x(t), constrained controls u(t) and the fixed prediction horizon T.
 * The functions V(t,x), L(t,x,u) and f(t,x,u) denote the terminal and integral
 * cost and the systems dynamics. Note that no terminal conditions for the states
 * are included in the problem formulation.
 *
 * The necessary optimality conditions can then be derived by means of the
 * Hamiltonian
 *
 * The function interfaces below have the following meaning (adj denotes the costates):
 *
 * sysfct:     f(t,x,u)
 *
 * sysjacxadj: df(t,x,u)/dx' * adj
 *
 * sysjacuadj: df(t,x,u)/du' * adj
 *
 * sysjacx:    df(t,x,u)/dx
 *
 * sysjacu:    df(t,x,u)/du
 *
 * icostfct:   L(t,x,u)
 *
 * icostjacx:  dL(t,x,u)/dx
 *
 * icostjacu:  dL(t,x,u)/du
 *
 * fcostfct:   V(t,x)
 *
 * fcostjacx:  dV(t,x)/dx
 *
 */
 

#include "probfct.h"

/* my own definitions: pointers to crane cost parameter vector */
#define iQ    0           /* position of integral weights (inputs) in pCost */
#define iR    iQ+10       /* position of integral weights (states) in pCost */
#define iP    iR+3        /* position of terminal weights (states) in pCost */

#define PI	3.141592653589793


void sysdim(typeInt *Nx, typeInt *Nu)
{
  *Nx = 10;
  *Nu = 3;
}


void sysfct(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *pSys)
{
  typeRNum t1 = 0.1e1 / x[1];
  typeRNum t2 = cos(x[4]);
  typeRNum t5 = x[6] * x[7];
  typeRNum t6 = sin(x[4]);
  typeRNum t7 = cos(x[3]);
  typeRNum t8 = t6 * t7;
  typeRNum t15 = t2 * x[1];
  typeRNum t23 = pow(x[7], 0.2e1);
  typeRNum t24 = t23 * t2;
  typeRNum t25 = x[1] * t7;
  typeRNum t26 = sin(x[3]);
  typeRNum t29 = t23 * x[0];
  typeRNum t45 = t7 * t7;
  typeRNum t49 = t2 * t2;
  typeRNum t54 = pow(x[8], 0.2e1);
  typeRNum t57 = t26 * t6;
  typeRNum t65 = -0.2e1 * x[5] * x[7] * t2 - 0.2e1 * x[6] * x[9] + 0.2e1 * t5 * t26 + t24 * x[1] * t45 * t6 + 0.2e1 * x[7] * t49 * t25 * x[8] - t15 * t54 * t6 + t29 * t57 - t8 * 9.81 - t57 * u[0] - u[2] * x[0] * t2 + u[2] * x[1] * t26;
  out[0] = x[5];
  out[1] = x[6];
  out[2] = x[7];
  out[3] = x[8];
  out[4] = x[9];
  out[5] = u[0];
  out[6] = u[1];
  out[7] = u[2];
  out[8] = -t1 / t2 * (0.2e1 * t5 * t8 + 0.2e1 * x[6] * t2 * x[8] + 0.2e1 * x[9] * x[7] * t15 * t7 - 0.2e1 * x[9] * x[1] * x[8] * t6 - t24 * t25 * t26 + t29 * t7 + t26 * 9.81 - t7 * u[0] + t6 * x[1] * t7 * u[2]);
  out[9] = t1 * t65;
}


void sysjacxadj(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *adj,
		        typeRNum *u, typeRNum *pSys)
{
  typeRNum t1 = 0.1e1 / x[1];
  typeRNum t2 = cos(x[4]);
  typeRNum t3 = 0.1e1 / t2;
  typeRNum t5 = pow(x[7], 0.2e1);
  typeRNum t6 = cos(x[3]);
  typeRNum t10 = sin(x[3]);
  typeRNum t12 = sin(x[4]);
  typeRNum t19 = x[6] * x[7];
  typeRNum t20 = t12 * t6;
  typeRNum t23 = x[6] * t2;
  typeRNum t26 = t5 * x[0];
  typeRNum t31 = pow(x[1], 0.2e1);
  typeRNum t32 = 0.1e1 / t31;
  typeRNum t34 = t3 * adj[8];
  typeRNum t36 = x[5] * x[7];
  typeRNum t43 = t10 * t12;
  typeRNum t47 = u[2] * x[0];
  typeRNum t60 = t5 * t2;
  typeRNum t61 = t60 * x[1];
  typeRNum t62 = t6 * t6;
  typeRNum t63 = x[1] * t62;
  typeRNum t76 = 0.2e1 * t19 * t6;
  typeRNum t80 = t2 * t2;
  typeRNum t82 = x[1] * t10;
  typeRNum t86 = t26 * t20;
  typeRNum t87 = t43 * 9.81;
  typeRNum t88 = t20 * u[0];
  typeRNum t90 = u[2] * x[1] * t6;
  typeRNum t95 = x[9] * x[1];
  typeRNum t110 = x[7] * t2;
  typeRNum t112 = t6 * x[8];
  typeRNum t116 = pow(x[8], 0.2e1);
  typeRNum t118 = t80 * x[1];
  typeRNum t121 = t10 * t2;
  typeRNum t154 = x[7] * x[0];
  typeRNum t178 = t110 * t6 - x[8] * t12;
  out[0] = -t1 * t3 * t5 * t6 * adj[8] + t1 * (t5 * t10 * t12 - u[2] * t2) * adj[9];
  out[1] = (0.2e1 * t19 * t20 + 0.2e1 * t23 * x[8] + t26 * t6 + t10 * 9.81 - t6 * u[0]) * t32 * t34 - (-0.2e1 * t36 * t2 - 0.2e1 * x[6] * x[9] + 0.2e1 * t19 * t10 + t26 * t43 - t20 * 9.81 - t43 * u[0] - t47 * t2) * t32 * adj[9];
  out[2] = 0.0e0;
  out[3] = t1 * (0.2e1 * t19 * t43 + 0.2e1 * x[9] * x[7] * t2 * x[1] * t10 - t61 + 0.2e1 * t60 * t63 + t26 * t10 - t6 * 9.81 - t10 * u[0] + t12 * x[1] * t10 * u[2]) * t34 + t1 * (t76 - 0.2e1 * t61 * t20 * t10 - 0.2e1 * x[7] * t80 * t82 * x[8] + t86 + t87 - t88 + t90) * adj[9];
  out[4] = -(t76 - 0.2e1 * t95 * x[8] + t86 + t87 - t88 + t90) * t1 / t80 * adj[8] + t1 * (0.2e1 * t36 * t12 - t5 * x[1] * t62 + 0.2e1 * t5 * t80 * t63 - 0.4e1 * t110 * x[1] * t112 * t12 + x[1] * t116 - 0.2e1 * t118 * t116 + t26 * t121 - t6 * t2 * 9.81 - t121 * u[0] + t47 * t12) * adj[9];
  out[5] = adj[0] - 0.2e1 * t1 * x[7] * t2 * adj[9];
  out[6] = adj[1] - 0.2e1 * t1 * (x[7] * t12 * t6 + t2 * x[8]) * t34 - 0.2e1 * (x[9] - x[7] * t10) * t1 * adj[9];
  out[7] = adj[2] - 0.2e1 * t1 * t6 * (x[6] * t12 + x[9] * t2 * x[1] - t110 * t82 + t154) * t3 * adj[8] + 0.2e1 * (-x[5] * t2 + x[6] * t10 + t110 * t63 * t12 + t118 * t112 + t154 * t43) * t1 * adj[9];
  out[8] = adj[3] + 0.2e1 * t1 * (-t23 + t95 * t12) * t34 + 0.2e1 * t2 * t178 * adj[9];
  out[9] = adj[4] - 0.2e1 * t178 * t3 * adj[8] - 0.2e1 * t1 * x[6] * adj[9];
}


void sysjacuadj(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *adj,
		        typeRNum *u, typeRNum *pSys)
{
  typeRNum t1 = 0.1e1 / x[1];
  typeRNum t2 = cos(x[4]);
  typeRNum t3 = 0.1e1 / t2;
  typeRNum t5 = cos(x[3]);
  typeRNum t6 = t5 * adj[8];
  typeRNum t8 = sin(x[3]);
  typeRNum t10 = sin(x[4]);
  out[0] = adj[5] + t1 * t3 * t6 - t1 * t8 * t10 * adj[9];
  out[1] = adj[6];
  out[2] = adj[7] - t3 * t10 * t6 - (x[0] * t2 - x[1] * t8) * t1 * adj[9];
}


/*
 * Matrix must be entered rowwise, i.e. A=[a11,a12;a21,a22] becomes
 * out[0]=a11; out[1]=a12; out[2]=a21; out[3]=a22;
 * Note that Matlab/Fortran does it columnwise !
 */
void sysjacx(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *pSys)
{
}


/*
 * Matrix must be entered rowwise, i.e. A=[a11,a12;a21,a22] becomes
 * out[0]=a11; out[1]=a12; out[2]=a21; out[3]=a22;
 * Note that Matlab/Fortran does it columnwise !
 */
void sysjacu(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *pSys)
{
}


void icostfct(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *xdes,
			  typeRNum *udes, typeRNum *pCost)
{
  typeInt i;
  typeRNum *Q = pCost+iQ;
  typeRNum *R = pCost+iR;

  out[0]=0;
  for (i = 0; i <= 9; i++)	{
    out[0] += Q[i]*(x[i]-xdes[i])*(x[i]-xdes[i]);
  }
  for (i = 0; i <= 2; i++)	{
    out[0] += R[i]*u[i]*u[i];
  }
  out[0] *= 0.5;
}


void icostjacx(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *xdes,
			   typeRNum *udes, typeRNum *pCost)
{
  typeInt i;
  typeRNum *Q = pCost+iQ;
  for (i = 0; i <= 9; i++) {
    out[i] = Q[i]*(x[i]-xdes[i]);
  }
}


void icostjacu(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *xdes,
			   typeRNum *udes, typeRNum *pCost)
{
  typeInt i;
  typeRNum *R = pCost+iR;
  for (i = 0; i <= 2; i++) {
    out[i] = R[i]*u[i];
  }
}


void fcostfct(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *xdes, typeRNum *pCost)
{
  typeInt i;
  typeRNum *P = pCost+iP;
  out[0]=0;	
  for (i = 0; i <= 9; i++) {
    out[0] += P[i]*(x[i]-xdes[i])*(x[i]-xdes[i]);
  }
  out[0] *= 0.5;
}


void fcostjacx(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *xdes, typeRNum *pCost)
{
  typeInt i;
  typeRNum *P = pCost+iP;
  for (i = 0; i <= 9; i++) {
    out[i] = P[i]*(x[i]-xdes[i]);
  }
}
