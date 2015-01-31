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
  *Nx = 4;
  *Nu = 2;
}


void sysfct(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *pSys)
{
  typeRNum k1, k2, h;

  k1 = pSys[0]*exp(-(pSys[2]/(273.15+x[2])));
  k2 = pSys[1]*exp(-(pSys[3]/(273.15+x[2])));
  h  = -pSys[7]*(k1*pSys[10]*x[0]+k1*pSys[11]*x[1]+k2*x[0]*x[0]*pSys[12]);

  out[0]=(-(k2*x[0]*x[0])+u[0]*(pSys[9]-x[0])-k1*x[0])/pSys[13];
  out[1]=(k1*x[0]-k1*x[1]-u[0]*x[1])/pSys[13];
  out[2]=(h+u[0]*(pSys[8]-x[2])+pSys[4]*(-x[2]+x[3]))/pSys[13];
  out[3]=(pSys[6]*u[1]+pSys[5]*(x[2]-x[3]))/pSys[13];
}


void sysjacxadj(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *adj, typeRNum *u, typeRNum *pSys)
{
  typeRNum k1, k2, dk1dx3, dk2dx3, dhdx1, dhdx2, dhdx3;

  k1     = pSys[0]*exp(-(pSys[2]/(273.15+x[2])));
  k2     = pSys[1]*exp(-(pSys[3]/(273.15+x[2])));
  dk1dx3 = exp(-(pSys[2]/(273.15+x[2])))*pSys[0]*pSys[2]/((273.15+x[2])*(273.15+x[2]));
  dk2dx3 = exp(-(pSys[3]/(273.15+x[2])))*pSys[1]*pSys[3]/((273.15+x[2])*(273.15+x[2]));
  dhdx1  = -(pSys[7]*(k1*pSys[10]+2.*k2*pSys[12]*x[0]));
  dhdx2  = -(k1*pSys[7]*pSys[11]);
  dhdx3  = -(pSys[7]*(dk2dx3*x[0]*x[0]*pSys[12]+dk1dx3*pSys[10]*x[0]+dk1dx3*pSys[11]*x[1]));

  out[0]=(-((u[0]+2.*k2*x[0])*adj[0])+k1*(-adj[0]+adj[1])+dhdx1*adj[2])/pSys[13];
  out[1]=(-((k1+u[0])*adj[1])+dhdx2*adj[2])/pSys[13];
  out[2]=(-(dk2dx3*x[0]*x[0]*adj[0])-dk1dx3*(x[0]*(adj[0]-adj[1])+x[1]*adj[1])+dhdx3*adj[2]-(pSys[4]+u[0])*adj[2]+pSys[5]*adj[3])/pSys[13];
  out[3]=(pSys[4]*adj[2]-pSys[5]*adj[3])/pSys[13];
}


void sysjacuadj(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *adj, typeRNum *u, typeRNum *pSys)
{
  out[0]=((pSys[9]-x[0])*adj[0])/pSys[13]-(x[1]*adj[1])/pSys[13]+((pSys[8]-x[2])*adj[2])/pSys[13];
  out[1]=(pSys[6]*adj[3])/pSys[13];
}


/*
 * Matrix must be entered rowwise, i.e. A=[a11,a12;a21,a22] becomes
 * out[0]=a11; out[1]=a12; out[2]=a21; out[3]=a22;
 * Note that Matlab/Fortran does it columnwise !
 */
void sysjacx(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *pSys)
{
  typeRNum k1, k2, dk1dx3, dk2dx3, dhdx1, dhdx2, dhdx3;

  k1     = pSys[0]*exp(-(pSys[2]/(273.15+x[2])));
  k2     = pSys[1]*exp(-(pSys[3]/(273.15+x[2])));
  dk1dx3 = exp(-(pSys[2]/(273.15+x[2])))*pSys[0]*pSys[2]/((273.15+x[2])*(273.15+x[2]));
  dk2dx3 = exp(-(pSys[3]/(273.15+x[2])))*pSys[1]*pSys[3]/((273.15+x[2])*(273.15+x[2]));
  dhdx1  = -(pSys[7]*(k1*pSys[10]+2.*k2*pSys[12]*x[0]));
  dhdx2  = -(k1*pSys[7]*pSys[11]);
  dhdx3  = -(pSys[7]*(dk2dx3*x[0]*x[0]*pSys[12]+dk1dx3*pSys[10]*x[0]+dk1dx3*pSys[11]*x[1]));

  out[0]=(-k1-u[0]-2.*k2*x[0])/pSys[13];
  out[1]=0.;
  out[2]=(-(dk2dx3*x[0]*x[0])-dk1dx3*x[0])/pSys[13];
  out[3]=0.;
  out[4]=k1/pSys[13];
  out[5]=(-k1-u[0])/pSys[13];
  out[6]=(dk1dx3*x[0]-dk1dx3*x[1])/pSys[13];
  out[7]=0.;
  out[8]=dhdx1/pSys[13];
  out[9]=dhdx2/pSys[13];
  out[10]=(dhdx3-pSys[4]-u[0])/pSys[13];
  out[11]=pSys[4]/pSys[13];
  out[12]=0.;
  out[13]=0.;
  out[14]=pSys[5]/pSys[13];
  out[15]=-(pSys[5]/pSys[13]);
}


/*
 * Matrix must be entered rowwise, i.e. A=[a11,a12;a21,a22] becomes
 * out[0]=a11; out[1]=a12; out[2]=a21; out[3]=a22;
 * Note that Matlab/Fortran does it columnwise !
 */
void sysjacu(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *pSys)
{
  out[0]=(pSys[9]-x[0])/pSys[13];
  out[1]=0.;
  out[2]=-(x[1]/pSys[13]);
  out[3]=0.;
  out[4]=(pSys[8]-x[2])/pSys[13];
  out[5]=0.;
  out[6]=0.;
  out[7]=pSys[6]/pSys[13];
}


void icostfct(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *xdes, typeRNum *udes, typeRNum *pCost)
{
  out[0]=pCost[8]*(u[0]-udes[0])*(u[0]-udes[0]) +
         pCost[9]*(u[1]-udes[1])*(u[1]-udes[1]) +
         pCost[4]*(x[0]-xdes[0])*(x[0]-xdes[0]) +
         pCost[5]*(x[1]-xdes[1])*(x[1]-xdes[1]) +
         pCost[6]*(x[2]-xdes[2])*(x[2]-xdes[2]) +
         pCost[7]*(x[3]-xdes[3])*(x[3]-xdes[3]);
}


void icostjacx(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *xdes, typeRNum *udes, typeRNum *pCost)
{
  out[0]=2.*pCost[4]*(x[0]-xdes[0]);
  out[1]=2.*pCost[5]*(x[1]-xdes[1]);
  out[2]=2.*pCost[6]*(x[2]-xdes[2]);
  out[3]=2.*pCost[7]*(x[3]-xdes[3]);
}


void icostjacu(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *xdes, typeRNum *udes, typeRNum *pCost)
{
  out[0]=2.*pCost[8]*(u[0]-udes[0]);
  out[1]=2.*pCost[9]*(u[1]-udes[1]);
}


void fcostfct(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *xdes, typeRNum *pCost)
{
  out[0]=pCost[0]*(x[0]-xdes[0])*(x[0]-xdes[0]) +
         pCost[1]*(x[1]-xdes[1])*(x[1]-xdes[1]) +
         pCost[2]*(x[2]-xdes[2])*(x[2]-xdes[2]) +
         pCost[3]*(x[3]-xdes[3])*(x[3]-xdes[3]);
}


void fcostjacx(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *xdes, typeRNum *pCost)
{
  out[0]=2.*pCost[0]*(x[0]-xdes[0]);
  out[1]=2.*pCost[1]*(x[1]-xdes[1]);
  out[2]=2.*pCost[2]*(x[2]-xdes[2]);
  out[3]=2.*pCost[3]*(x[3]-xdes[3]);
}

