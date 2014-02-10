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
 * Problem function for GRAMPC.
 * Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Date: February 2014
 * Version: v1.0
 *
 *
 * This probfct file provides an interface to GRAMPC. The underlying
 * optimal control problem (OCP) of the model predictive control (MPC) formulation
 * has the following structure
 *                                  _T
 *                                 /
 *      min    J(u) = V(x(T),t) + / L(x(t),u(t),t) dt
 *      u(.)                    _/
 *                             0
 *             .
 *      s.t.   x(t) = f(x(t),u(t),t), x(0) = x0
 *
 *             Ul <= u(t) <= Uu, t in [0,T]
 *
 * with states x(t), constrained controls u(t) and the fixed prediction horizon T.
 * The functions V(x,t), L(x,u,t) and f(x,u,t) denote the terminal and integral
 * cost and the systems dynamics. Note that no terminal conditions for the states
 * are included in the problem formulation.
 *
 * The function interfaces below have the following meaning:
 *
 * sysfct:     f(x,u,t)
 *
 * sysjacxadj: df(x,u,t)/dx' * adj
 *
 * sysjacuadj: df(x,u,t)/du' * adj
 *
 * sysjacx:    df(x,u,t)/dx
 *
 * sysjacu:    df(x,u,t)/du
 *
 * icostfct:   L(x,u,t)
 *
 * icostjacx:  dL(x,u,t)/dx
 *
 * icostjacu:  dL(x,u,t)/du
 *
 * fcostfct:   V(x,t)
 *
 * fcostjacx:  dV(x,t)/dx
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

