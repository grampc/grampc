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


void sysdim(typeInt *Nx, typeInt *Nu)
{
  *Nx = ;
  *Nu = ;
}


void sysfct(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *pSys)
{
}


void sysjacxadj(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *adj, typeRNum *u, typeRNum *pSys)
{
}


void sysjacuadj(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *adj, typeRNum *u, typeRNum *pSys)
{
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


void icostfct(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *xdes, typeRNum *udes, typeRNum *pCost)
{
}


void icostjacx(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *xdes, typeRNum *udes, typeRNum *pCost)
{
}


void icostjacu(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *xdes, typeRNum *udes, typeRNum *pCost)
{
}


void fcostfct(typeRNum *out, typeRNum t,  typeRNum *x, typeRNum *xdes, typeRNum *pCost)
{
}


void fcostjacx(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *xdes, typeRNum *pCost)
{
}

