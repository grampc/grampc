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


#define NN   5
#define NX   2*NN
#define NU   2

#define iQ  0
#define iR  iQ+NX
#define iP  iR+NU

#include "probfct.h"


void sysdim(typeInt *Nx, typeInt *Nu)
{
  *Nx = NX;
  *Nu = NU;
}


void sysfct(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *pSys)
{
  const typeRNum m = pSys[0];
  const typeRNum c = pSys[1];
  const typeRNum d = pSys[2];
  typeInt k;

  for (k=0; k<=NN-1; k++) {
    out[k] = x[NN+k];
  }

  out[NN] = -2*c/m * x[0] + c/m * x[1] - 2*d/m * x[NN] + d/m * x[NN+1] + 1/m * u[0];
  
  for (k=1; k<=NN-2; k++) {
    out[NN+k] = c/m * x[k-1] - 2*c/m * x[k] + c/m * x[k+1] + d/m * x[NN+k-1] - 2*d/m * x[NN+k] + d/m * x[NN+k+1];
  }

  out[2*NN-1] = c/m * x[NN-2] - 2*c/m * x[NN-1] + d/m * x[2*NN-2] - 2*d/m * x[2*NN-1] - 1/m * u[1];
}


void sysjacxadj(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *adj,
		        typeRNum *u, typeRNum *pSys)
{
  const typeRNum m = pSys[0];
  const typeRNum c = pSys[1];
  const typeRNum d = pSys[2];
  typeInt k;

  out[0] = -2*c/m * adj[NN] + c/m * adj[NN+1];

  for (k=1; k<=NN-2; k++) {
    out[k] = c/m * adj[NN+k-1] - 2*c/m * adj[NN+k] + c/m * adj[NN+k+1];
  }

  out[NN-1] = c/m * adj[2*NN-2] - 2*c/m * adj[2*NN-1];

  out[NN] = adj[0] - 2*d/m * adj[NN] + d/m * adj[NN+1];

  for (k=1; k<=NN-2; k++) {
    out[NN+k] = adj[k] +  d/m * adj[NN+k-1] - 2*d/m * adj[NN+k] + d/m * adj[NN+k+1];
  }

  out[2*NN-1] = adj[NN-1] +  d/m * adj[2*NN-2] - 2*d/m * adj[2*NN-1];

}


void sysjacuadj(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *adj,
		        typeRNum *u, typeRNum *pSys)
{
  const typeRNum m = pSys[0]; 

  out[0] =  1/m * adj[NN];
  out[1] = -1/m * adj[2*NN-1];
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
  for (i=0;i<=NX-1;i++)   out[0] += Q[i]*(x[i]-xdes[i])*(x[i]-xdes[i]);
  for (i=0;i<=NU-1;i++)   out[0] += R[i]*(u[i]-udes[i])*(u[i]-udes[i]);
  out[0] *= 0.5;
}


void icostjacx(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *xdes,
			   typeRNum *udes, typeRNum *pCost)
{
  typeInt i;
  typeRNum *Q = pCost+iQ;    
  for (i=0;i<=NX-1;i++)  out[i] = Q[i]*(x[i]-xdes[i]); 
}


void icostjacu(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *xdes,
			   typeRNum *udes, typeRNum *pCost)
{
  typeInt i;
  typeRNum *R = pCost+iR;
  for (i=0;i<=NU-1;i++)  out[i] = R[i]*(u[i]-udes[i]);
}


void fcostfct(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *xdes, typeRNum *pCost)
{
  typeInt i;
  typeRNum *P = pCost+iP;
  out[0]=0;
  for (i=0;i<=NX-1;i++)   out[0] += P[i]*(x[i]-xdes[i])*(x[i]-xdes[i]);
  out[0] *= 0.5;
}


void fcostjacx(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *xdes, typeRNum *pCost)
{
  typeInt i;
  typeRNum *P = pCost+iP;
  for (i=0;i<=NX-1;i++)   out[i] = P[i]*(x[i]-xdes[i]);
}
