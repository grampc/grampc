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
 * derivatives of a quadrotor with 9 states and 4 controls. For a more
 * detailed model see
 *
 * M. Hehn and R. D'Andrea, "A flying inverted pendulum", Proceedings of the IEEE
 * International Conference on Robotics and Automation, pp. 763-770, 2011.
 *
 * The functions were computed and exported by means of the computer algebra 
 * program MATHEMATICA. 
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
  *Nx = 9;
  *Nu = 4;
}


void sysfct(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *pSys)
{
  double aux1, aux2, aux3, aux4, aux5, aux6, aux7;
    
  aux1=cos(x[8]);
  aux2=sin(x[6]);
  aux3=cos(x[6]);
  aux4=sin(x[7]);
  aux5=sin(x[8]);
  aux6=cos(x[7]);
  aux7=tan(x[7]);
    
  out[0]=x[1];
  out[1]=u[0]*(aux3*aux1*aux4+aux5*aux2);
  out[2]=x[3];
  out[3]=u[0]*(aux5*aux3*aux4-aux1*aux2);
  out[4]=x[5];
  out[5]=-9.81+u[0]*aux6*aux3;
  out[6]=(1.*(u[1]*aux3+u[2]*aux2))/aux6;
  out[7]=u[2]*aux3-u[1]*aux2;
  out[8]=u[1]*aux3*aux7+u[2]*aux7*aux2+u[3];
}


void sysjacxadj(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *adj, typeRNum *u, typeRNum *pSys)
{
  double aux1, aux2, aux3, aux4, aux5, aux6, aux7, aux8, aux9, aux10, aux11;
    
  aux1=cos(x[7]);
  aux2=1/aux1;
  aux3=tan(x[7]);
  aux4=cos(x[8]);
  aux5=sin(x[8]);
  aux6=aux2*aux2;
  aux7=sin(x[6]);
  aux8=cos(x[6]);
  aux9=sin(x[7]);
  aux10=adj[6]*aux9;
  aux11=adj[8]+aux10;
    
  out[0]=0.;
  out[1]=adj[0];
  out[2]=0.;
  out[3]=adj[2];
  out[4]=0.;
  out[5]=adj[4];
  out[6]=aux8*(1.*aux2*u[2]*adj[6]-u[1]*adj[7]+u[2]*aux3*adj[8]-aux4\
	       *adj[3]*u[0]+aux5*u[0]*adj[1])-cos(x[7])*aux7*(1.*aux2*u[2]*adj[7]+\
							      adj[6]*u[1]*aux6+adj[5]*u[0]+aux3*(1.*aux2*u[1]*adj[8]+adj[3]*\
												 aux5*u[0]+aux4*u[0]*adj[1]));
  out[7]=1.*u[2]*aux7*aux6*aux11+aux8*(1.*u[1]*aux6*aux11-adj[5]*\
				       aux9*u[0]+cos(x[7])*u[0]*(adj[3]*aux5+aux4*adj[1]));
  out[8]=u[0]*(aux4*(adj[3]*aux9*aux8+aux7*adj[1])+aux5*(adj[3]*aux7-\
            aux9*aux8*adj[1]));
}


void sysjacuadj(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *adj, typeRNum *u, typeRNum *pSys)
{
    double aux1, aux2, aux3, aux4, aux5, aux6, aux7, aux8, aux9;
    
    aux1=cos(x[8]);
    aux2=sin(x[8]);
    aux3=sin(x[6]);
    aux4=cos(x[6]);
    aux5=cos(x[7]);
    aux6=sin(x[7]);
    aux7=1/aux5;
    aux8=aux6*adj[8];
    aux9=adj[6]+aux8;
    
    out[0]=aux3*(aux2*adj[1]-aux1*adj[3])+aux4*(aux6*(adj[1]*aux1+aux2*adj[3])+aux5*adj[5]);
    out[1]=-(aux3*adj[7])+1.*aux4*aux7*aux9;
    out[2]=aux4*adj[7]+1.*aux3*aux7*aux9;
    out[3]=adj[8];
}


/*
 * Matrix must be entered rowwise, i.e. A=[a11,a12;a21,a22] becomes
 * out[0]=a11; out[1]=a12; out[2]=a21; out[3]=a22;
 * Note that Matlab/Fortran does it columnwise !
 */
void sysjacx(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *pSys)
{
  out[0]=0.;
  out[1]=1.;
  out[2]=0.;
  out[3]=0.;
  out[4]=0.;
  out[5]=0.;
  out[6]=0.;
  out[7]=0.;
  out[8]=0.;
  out[9]=0.;
  out[10]=0.;
  out[11]=0.;
  out[12]=0.;
  out[13]=0.;
  out[14]=0.;
  out[15]=(-(cos(x[8])*sin(x[6])*sin(x[7]))+cos(x[6])*sin(x[8]))*u[0];
  out[16]=cos(x[6])*cos(x[7])*cos(x[8])*u[0];
  out[17]=(cos(x[8])*sin(x[6])-cos(x[6])*sin(x[7])*sin(x[8]))*u[0];
  out[18]=0.;
  out[19]=0.;
  out[20]=0.;
  out[21]=1.;
  out[22]=0.;
  out[23]=0.;
  out[24]=0.;
  out[25]=0.;
  out[26]=0.;
  out[27]=0.;
  out[28]=0.;
  out[29]=0.;
  out[30]=0.;
  out[31]=0.;
  out[32]=0.;
  out[33]=(-(cos(x[6])*cos(x[8]))-sin(x[6])*sin(x[7])*sin(x[8]))*u[0];
  out[34]=cos(x[6])*cos(x[7])*sin(x[8])*u[0];
  out[35]=(cos(x[6])*cos(x[8])*sin(x[7])+sin(x[6])*sin(x[8]))*u[0];
  out[36]=0.;
  out[37]=0.;
  out[38]=0.;
  out[39]=0.;
  out[40]=0.;
  out[41]=1.;
  out[42]=0.;
  out[43]=0.;
  out[44]=0.;
  out[45]=0.;
  out[46]=0.;
  out[47]=0.;
  out[48]=0.;
  out[49]=0.;
  out[50]=0.;
  out[51]=-(cos(x[7])*sin(x[6])*u[0]);
  out[52]=-(cos(x[6])*sin(x[7])*u[0]);
  out[53]=0.;
  out[54]=0.;
  out[55]=0.;
  out[56]=0.;
  out[57]=0.;
  out[58]=0.;
  out[59]=0.;
  out[60]=(1.*(-(sin(x[6])*u[1])+cos(x[6])*u[2]))/cos(x[7]);
  out[61]=(1.*tan(x[7])*(cos(x[6])*u[1]+sin(x[6])*u[2]))/cos(x[7]);
  out[62]=0.;
  out[63]=0.;
  out[64]=0.;
  out[65]=0.;
  out[66]=0.;
  out[67]=0.;
  out[68]=0.;
  out[69]=-(cos(x[6])*u[1])-sin(x[6])*u[2];
  out[70]=0.;
  out[71]=0.;
  out[72]=0.;
  out[73]=0.;
  out[74]=0.;
  out[75]=0.;
  out[76]=0.;
  out[77]=0.;
  out[78]=-(sin(x[6])*tan(x[7])*u[1])+cos(x[6])*tan(x[7])*u[2];
  out[79]=1.*cos(x[6])*pow(cos(x[7]),-2.)*u[1]+1.*pow(cos(x[7]),-2.)*sin(x[6])*u[2];
  out[80]=0.;
}


/*
 * Matrix must be entered rowwise, i.e. A=[a11,a12;a21,a22] becomes
 * out[0]=a11; out[1]=a12; out[2]=a21; out[3]=a22;
 * Note that Matlab/Fortran does it columnwise !
 */
void sysjacu(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *pSys)
{
  out[0]=0.;
  out[1]=0.;
  out[2]=0.;
  out[3]=0.;
  out[4]=cos(x[6])*cos(x[8])*sin(x[7])+sin(x[6])*sin(x[8]);
  out[5]=0.;
  out[6]=0.;
  out[7]=0.;
  out[8]=0.;
  out[9]=0.;
  out[10]=0.;
  out[11]=0.;
  out[12]=-(cos(x[8])*sin(x[6]))+cos(x[6])*sin(x[7])*sin(x[8]);
  out[13]=0.;
  out[14]=0.;
  out[15]=0.;
  out[16]=0.;
  out[17]=0.;
  out[18]=0.;
  out[19]=0.;
  out[20]=cos(x[6])*cos(x[7]);
  out[21]=0.;
  out[22]=0.;
  out[23]=0.;
  out[24]=0.;
  out[25]=(1.*cos(x[6]))/cos(x[7]);
  out[26]=(1.*sin(x[6]))/cos(x[7]);
  out[27]=0.;
  out[28]=0.;
  out[29]=-sin(x[6]);
  out[30]=cos(x[6]);
  out[31]=0.;
  out[32]=0.;
  out[33]=cos(x[6])*tan(x[7]);
  out[34]=sin(x[6])*tan(x[7]);
  out[35]=1.;
}


void icostfct(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *xdes, typeRNum *udes, typeRNum *pCost)
{
  out[0]=pCost[18]*(u[0]-udes[0])*(u[0]-udes[0]) +
         pCost[19]*(u[1]-udes[1])*(u[1]-udes[1]) +
         pCost[20]*(u[2]-udes[2])*(u[2]-udes[2]) +
         pCost[21]*(u[3]-udes[3])*(u[3]-udes[3]) +
         pCost[9] *(x[0]-xdes[0])*(x[0]-xdes[0]) +
         pCost[10]*(x[1]-xdes[1])*(x[1]-xdes[1]) +
         pCost[11]*(x[2]-xdes[2])*(x[2]-xdes[2]) +
         pCost[12]*(x[3]-xdes[3])*(x[3]-xdes[3]) +
         pCost[13]*(x[4]-xdes[4])*(x[4]-xdes[4]) +
         pCost[14]*(x[5]-xdes[5])*(x[5]-xdes[5]) +
         pCost[15]*(x[6]-xdes[6])*(x[6]-xdes[6]) +
         pCost[16]*(x[7]-xdes[7])*(x[7]-xdes[7]) +
         pCost[17]*(x[8]-xdes[8])*(x[8]-xdes[8]);
}


void icostjacx(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *xdes, typeRNum *udes, typeRNum *pCost)
{
  out[0]=2.*pCost[9] *(x[0]-xdes[0]);
  out[1]=2.*pCost[10]*(x[1]-xdes[1]);
  out[2]=2.*pCost[11]*(x[2]-xdes[2]);
  out[3]=2.*pCost[12]*(x[3]-xdes[3]);
  out[4]=2.*pCost[13]*(x[4]-xdes[4]);
  out[5]=2.*pCost[14]*(x[5]-xdes[5]);
  out[6]=2.*pCost[15]*(x[6]-xdes[6]);
  out[7]=2.*pCost[16]*(x[7]-xdes[7]);
  out[8]=2.*pCost[17]*(x[8]-xdes[8]);
}


void icostjacu(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *xdes, typeRNum *udes, typeRNum *pCost)
{
  out[0]=2.*pCost[18]*(u[0]-udes[0]);
  out[1]=2.*pCost[19]*(u[1]-udes[1]);
  out[2]=2.*pCost[20]*(u[2]-udes[2]);
  out[3]=2.*pCost[21]*(u[3]-udes[3]);
}


void fcostfct(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *xdes, typeRNum *pCost)
{
  out[0]=pCost[0]*(x[0]-xdes[0])*(x[0]-xdes[0]) +
         pCost[1]*(x[1]-xdes[1])*(x[1]-xdes[1]) +
         pCost[2]*(x[2]-xdes[2])*(x[2]-xdes[2]) +
         pCost[3]*(x[3]-xdes[3])*(x[3]-xdes[3]) +
         pCost[4]*(x[4]-xdes[4])*(x[4]-xdes[4]) +
         pCost[5]*(x[5]-xdes[5])*(x[5]-xdes[5]) +
         pCost[6]*(x[6]-xdes[6])*(x[6]-xdes[6]) +
         pCost[7]*(x[7]-xdes[7])*(x[7]-xdes[7]) +
         pCost[8]*(x[8]-xdes[8])*(x[8]-xdes[8]);
}


void fcostjacx(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *xdes, typeRNum *pCost)
{
  out[0]=2.*pCost[0]*(x[0]-xdes[0]);
  out[1]=2.*pCost[1]*(x[1]-xdes[1]);
  out[2]=2.*pCost[2]*(x[2]-xdes[2]);
  out[3]=2.*pCost[3]*(x[3]-xdes[3]);
  out[4]=2.*pCost[4]*(x[4]-xdes[4]);
  out[5]=2.*pCost[5]*(x[5]-xdes[5]);
  out[6]=2.*pCost[6]*(x[6]-xdes[6]);
  out[7]=2.*pCost[7]*(x[7]-xdes[7]);
  out[8]=2.*pCost[8]*(x[8]-xdes[8]);
}

