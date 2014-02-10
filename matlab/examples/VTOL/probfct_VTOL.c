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
 * Probfunction for grampc toolbox.
 *
 *
 * This probfct file describes the dynamics, the cost function and the corresponding
 * derivatives of a vertical takeoff and landing (VTOL) aircraft with 6 states and 2
 * controls. For a more detailed model, see, for instance,
 *
 * S. Sastry, "Nonlinear Systems - Analysis, Stability and Control", Springer-Verlag,
 * New-York, pp. 469-481, 1999.
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
 * The function interfaces below have the following meaning (adj denotes the costates):
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
  *Nx = 6;
  *Nu = 2;
}


void sysfct(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *pSys)
{
  out[0]=x[1];
  out[1]=-(sin(x[4])*u[0])+cos(x[4])*pSys[0]*u[1];
  out[2]=x[3];
  out[3]=-1.+cos(x[4])*u[0]+pSys[0]*sin(x[4])*u[1];
  out[4]=x[5];
  out[5]=u[1];
}


void sysjacxadj(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *adj, typeRNum *u, typeRNum *pSys)
{
  out[0]=0.;
  out[1]=adj[0];
  out[2]=0.;
  out[3]=adj[2];
  out[4]=(-(cos(x[4])*u[0])-pSys[0]*sin(x[4])*u[1])*adj[1]+(-(sin(x[4])*u[0])+cos(x[4])*pSys[0]*u[1])*adj[3];
  out[5]=adj[4];
}


void sysjacuadj(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *adj, typeRNum *u, typeRNum *pSys)
{
  out[0]=-(sin(x[4])*adj[1])+cos(x[4])*adj[3];
  out[1]=cos(x[4])*pSys[0]*adj[1]+pSys[0]*sin(x[4])*adj[3]+adj[5];
}


/*
 * Matrix must be entered rowwise, i.e. A=[a11,a12;a21,a22] becomes
 * out[0]=a11; out[1]=a12; out[2]=a21; out[3]=a22;
 * Note that Matlab/Fortran does it columnwise !
 */
void sysjacx(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *pSys)
{
  out[0] =0.; 
  out[1] =1.; 
  out[2] =0.; 
  out[3] =0.; 
  out[4] =0.;
  out[5] =0.;
	
  out[6] =0.; 
  out[7] =0.; 
  out[8] =0.; 
  out[9] =0.; 
  out[10]=-(cos(x[4])*u[0])-pSys[0]*sin(x[4])*u[1]; 
  out[11]=0.;
	
  out[12]=0.; 
  out[13]=0.; 
  out[14]=0.; 
  out[15]=1.; 
  out[16]=0.;                                       
  out[17]=0.;
	
  out[18]=0.; 
  out[19]=0.; 
  out[20]=0.; 
  out[21]=0.; 
  out[22]=-(sin(x[4])*u[0])+cos(x[4])*pSys[0]*u[1]; 
  out[23]=0.;
	
  out[24]=0.; 
  out[25]=0.; 
  out[26]=0.; 
  out[27]=0.; 
  out[28]=0.;                                       
  out[29]=1.;
	
  out[30]=0.; 
  out[31]=0.; 
  out[32]=0.; 
  out[33]=0.; 
  out[34]=0.;                                       
  out[35]=0.;
}


/*
 * Matrix must be entered rowwise, i.e. A=[a11,a12;a21,a22] becomes
 * out[0]=a11; out[1]=a12; out[2]=a21; out[3]=a22;
 * Note that Matlab/Fortran does it columnwise !
 */
void sysjacu(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *pSys)
{
  out[0] =0.;         
  out[1] =0.;
	
  out[2] =-sin(x[4]); 
  out[3] =cos(x[4])*pSys[0];
	
  out[4] =0.;         
  out[5] =0.;
	
  out[6] =cos(x[4]);  
  out[7] =pSys[0]*sin(x[4]);
	
  out[8] =0.;         
  out[9] =0.;
	
  out[10]=0.;         
  out[11]=1.;
}


void icostfct(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *xdes, typeRNum *udes, typeRNum *pCost)
{
  out[0]=pCost[12]*(u[0]-udes[0])*(u[0]-udes[0]) +
         pCost[13]*(u[1]-udes[1])*(u[1]-udes[1]) +
         pCost[6] *(x[0]-xdes[0])*(x[0]-xdes[0]) +
         pCost[7] *(x[1]-xdes[1])*(x[1]-xdes[1]) +
         pCost[8] *(x[2]-xdes[2])*(x[2]-xdes[2]) +
         pCost[9] *(x[3]-xdes[3])*(x[3]-xdes[3]) +
         pCost[10]*(x[4]-xdes[4])*(x[4]-xdes[4]) +
         pCost[11]*(x[5]-xdes[5])*(x[5]-xdes[5]);
}


void icostjacx(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *xdes, typeRNum *udes, typeRNum *pCost)
{
  out[0]=2.*pCost[6] *(x[0]-xdes[0]);
  out[1]=2.*pCost[7] *(x[1]-xdes[1]);
  out[2]=2.*pCost[8] *(x[2]-xdes[2]);
  out[3]=2.*pCost[9] *(x[3]-xdes[3]);
  out[4]=2.*pCost[10]*(x[4]-xdes[4]);
  out[5]=2.*pCost[11]*(x[5]-xdes[5]);
}


void icostjacu(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *xdes, typeRNum *udes, typeRNum *pCost)
{
  out[0]=2.*pCost[12]*(u[0]-udes[0]);
  out[1]=2.*pCost[13]*(u[1]-udes[1]);
}


void fcostfct(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *xdes, typeRNum *pCost)
{
  out[0]=pCost[0]*(x[0]-xdes[0])*(x[0]-xdes[0]) +
         pCost[1]*(x[1]-xdes[1])*(x[1]-xdes[1]) +
         pCost[2]*(x[2]-xdes[2])*(x[2]-xdes[2]) +
         pCost[3]*(x[3]-xdes[3])*(x[3]-xdes[3]) +
         pCost[4]*(x[4]-xdes[4])*(x[4]-xdes[4]) +
         pCost[5]*(x[5]-xdes[5])*(x[5]-xdes[5]);
}


void fcostjacx(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *xdes, typeRNum *pCost)
{
  out[0]=2.*pCost[0]*(x[0]-xdes[0]);
  out[1]=2.*pCost[1]*(x[1]-xdes[1]);
  out[2]=2.*pCost[2]*(x[2]-xdes[2]);
  out[3]=2.*pCost[3]*(x[3]-xdes[3]);
  out[4]=2.*pCost[4]*(x[4]-xdes[4]);
  out[5]=2.*pCost[5]*(x[5]-xdes[5]);
}

