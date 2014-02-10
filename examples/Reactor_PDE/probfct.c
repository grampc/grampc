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
 * derivatives of a quasilinear diffusion-convection-reaction system.
 * More details can be found in
 *
 * U. Tilman, K. Graichen and A. Kugi, "Trajectory planning and receding horizon
 * tracking control of a quasilinear diffusion-convection-reaction system",
 * Proceedings of the 8th IFAC Symposium on Nonlinear Control Systems, Bologna (Italy),
 * pp. 587-592, 2010.
 *
 * The functions were computed and exported by means of the computer algebra 
 * program MAPLE. 
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
 

#define NX  10
#define NU  1

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
  typeRNum q0  = pSys[0];
  typeRNum q1  = pSys[1];
  typeRNum nu  = pSys[2];
  typeRNum r0  = pSys[3];
  typeRNum r1  = pSys[4];
  typeRNum rb0 = pSys[5];
  typeRNum rb1 = pSys[6];
  typeInt N = NX-1;
  typeInt k;

  typeRNum dz = 1.0/N;

  typeRNum t1,t4,t7,t13,t15,t21,t27,t31,t36,t41,t44,t54,t59,t64,t70,t74,t82,t86,t92,t99,t110,t116,t122,t128;


  t1 = x[0];
  t4 = x[1];
  t7 = t1 * t1;
  t13 = t4 * t4;
  t15 = dz * dz;
  t21 = 0.1e1 / t15;
  out[0] = (-0.2e1 * q0 * t1 + 0.2e1 * q0 * t4 - 0.3e1 * q1 * t7 + 0.4e1 * q1 * t1 * t4 - q1 * t13 + t1 * t15 * r0 + t7 * t15 * r1) * t21 / (1.0 + 0.0 * t1);

  for ( k=1 ; k < N ; k++) {
    t27 = x[k];
    t31 = t27 * t27;
    t36 = x[k + 1];
    t41 = x[k - 1];
    t44 = q1 * t27;
    t54 = nu * dz;
    t59 = 0.0 * t27;
    t64 = 0.2e1 * t27 * t15 * r0 + 0.2e1 * t31 * t15 * r1 + 0.2e1 * q0 * t36 - 0.4e1 * q0 * t27 + 0.2e1 * q0 * t41 + 0.4e1 * t44 * t36 - 0.6e1 * q1 * t31 + 0.4e1 * t44 * t41 - 0.2e1 * q1 * t36 * t41 - t54 * 1.0 * t36 + t54 * 1.0 * t41 - t54 * t59 * t36 + t54 * t59 * t41;
    out[k] = t64 * t21 / (1.0 + t59) / 0.2e1;
  }

  t70 = x[N];
  t74 = t70 * t70;
  t82 = x[N - 1];
  t86 = t70 * dz;
  t92 = q1 * t70;
  t99 = q1 * t74;
  t110 = t82 * t82;
  t116 = nu * t15;
  t122 = 0.0 * t70;
  t128 = t70 * rb0 * t15 * r0 + t74 * rb0 * t15 * r1 + 0.2e1 * q0 * u[0] * dz + 0.2e1 * q0 * rb0 * t82 - 0.2e1 * q0 * rb1 * t86 - 0.2e1 * q0 * t70 * rb0 + 0.4e1 * t92 * u[0] * dz + 0.4e1 * t92 * rb0 * t82 - 0.4e1 * t99 * rb1 * dz - 0.3e1 * t99 * rb0 - 0.2e1 * q1 * u[0] * dz * t82 - q1 * rb0 * t110 + 0.2e1 * q1 * rb1 * t86 * t82 - t116 * 1.0 * u[0] + t116 * 1.0 * rb1 * t70 - t116 * t122 * u[0] + t116 * 0.0 * t74 * rb1;
  out[N] = t128 / rb0 * t21 / (1.0 + t122);
}


void sysjacxadj(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *adj,
		        typeRNum *u, typeRNum *pSys)
{
  typeRNum q0  = pSys[0];
  typeRNum q1  = pSys[1];
  typeRNum nu  = pSys[2];
  typeRNum r0  = pSys[3];
  typeRNum r1  = pSys[4];
  typeRNum rb0 = pSys[5];
  typeRNum rb1 = pSys[6];
  typeInt N = NX-1;
  typeInt k;

  typeRNum dz = 1.0/N;

  typeRNum dfdx0,dfdx1,dfdx2;

  typeRNum t1,t2,t4,t6,t8,t12,t15,t16,t22,t27,t34;
  typeRNum t3,t7,t9,t11,t14,t20,t26,t37,t63,t65;
  typeRNum t10,t13,t17,t19,t21,t24,t28,t29,t32,t39,t42,t43,t50,t55,t74,t82,t96,t98;

  for (k=0;k<=N;k++) out[k]=0.0;

  t1 = 2 * q0;
  t2 = q1 * x[0];
  t4 = q1 * x[1];
  t6 = dz * dz;
  t8 = x[0] * t6;
  t12 = 0.1e1 / t6;
  t15 = 1.0 + 0.0 * x[0];
  t16 = 0.1e1 / t15;
  t22 = pow(x[0], 0.2e1);
  t27 = pow(x[1], 0.2e1);
  t34 = t15 * t15;
  dfdx0 = (-(typeRNum) t1 - 0.6e1 * t2 + 0.4e1 * t4 + t6 * r0 + 0.2e1 * t8 * r1) * t12 * t16 - (-0.2e1 * (typeRNum) q0 * x[0] + 0.2e1 * (typeRNum) q0 * x[1] - 0.3e1 * q1 * t22 + 0.4e1 * t2 * x[1] - q1 * t27 + t8 * r0 + t22 * t6 * r1) * t12 / t34 * 0.0;
  dfdx1 = ((typeRNum) t1 + 0.4e1 * t2 - 0.2e1 * t4) * t12 * t16;

  out[0] += dfdx0*adj[0];
  out[1] += dfdx1*adj[0];

  for (k=1;k<N;k++) {
    t1 = 2 * q0;
    t2 = q1 * x[k];
    t3 = 0.4e1 * t2;
    t4 = q1 * x[k + 1];
    t6 = nu * dz;
    t7 = t6 * 1.0;
    t8 = 0.0 * x[k];
    t9 = t6 * t8;
    t11 = dz * dz;
    t12 = 0.1e1 / t11;
    t14 = 1.0 + t8;
    t15 = 0.1e1 / t14;
    t20 = x[k] * t11;
    t26 = q1 * x[k - 1];
    t37 = pow(x[k], 0.2e1);
    t63 = 0.2e1 * t20 * r0 + 0.2e1 * t37 * t11 * r1 + 0.2e1 * (typeRNum) q0 * x[k + 1] - 0.4e1 * (typeRNum) q0 * x[k] + 0.2e1 * (typeRNum) q0 * x[k - 1] + 0.4e1 * t2 * x[k + 1] - 0.6e1 * q1 * t37 + 0.4e1 * t2 * x[k - 1] - 0.2e1 * t4 * x[k - 1] - t6 * 1.0 * x[k + 1] + t6 * 1.0 * x[k - 1] - t6 * t8 * x[k + 1] + t6 * t8 * x[k - 1];
    t65 = t14 * t14;
    dfdx0 = ((typeRNum) t1 + t3 - 0.2e1 * t4 + t7 + t9) * t12 * t15 / 0.2e1;
    dfdx1 = (0.2e1 * t11 * r0 + 0.4e1 * t20 * r1 - (typeRNum) (4 * q0) + 0.4e1 * t4 - 0.12e2 * t2 + 0.4e1 * t26 - t6 * 0.0 * x[k + 1] + t6 * 0.0 * x[k - 1]) * t12 * t15 / 0.2e1 - t63 * t12 / t65 * 0.0 / 0.2e1;
    dfdx2 = ((typeRNum) t1 + t3 - 0.2e1 * t26 - t7 - t9) * t12 * t15 / 0.2e1;

    out[k-1] += dfdx0*adj[k];
    out[k  ] += dfdx1*adj[k];
    out[k+1] += dfdx2*adj[k];
  }

  t1 = q0 * rb0;
  t2 = 0.2e1 * t1;
  t3 = q1 * x[N];
  t4 = t3 * rb0;
  t6 = q1 * u[0];
  t7 = t6 * dz;
  t9 = q1 * rb0;
  t10 = t9 * x[N - 1];
  t12 = q1 * rb1;
  t13 = x[N] * dz;
  t14 = t12 * t13;
  t17 = 0.1e1 / rb0;
  t19 = dz * dz;
  t20 = 0.1e1 / t19;
  t21 = 0.0 * x[N];
  t22 = 1.0 + t21;
  t24 = t20 / t22;
  t28 = x[N] * rb0;
  t29 = t19 * r1;
  t32 = q0 * rb1;
  t39 = dz * x[N - 1];
  t42 = nu * t19;
  t43 = 1.0 * rb1;
  t50 = rb0 * t19 * r0 + 0.2e1 * t28 * t29 - 0.2e1 * t32 * dz - t2 + 0.4e1 * t7 + 0.4e1 * t10 - 0.8e1 * t14 - 0.6e1 * t4 + 0.2e1 * t12 * t39 + t42 * t43 - t42 * 0.0 * u[0] + 0.2e1 * t42 * t21 * rb1;
  t55 = pow(x[N], 0.2e1);
  t74 = q1 * t55;
  t82 = pow(x[N - 1], 0.2e1);
  t96 = t28 * t19 * r0 + t55 * rb0 * t29 + 0.2e1 * q0 * u[0] * dz + 0.2e1 * t1 * x[N - 1] - 0.2e1 * t32 * t13 - 0.2e1 * q0 * x[N] * rb0 + 0.4e1 * t3 * u[0] * dz + 0.4e1 * t3 * rb0 * x[N - 1] - 0.4e1 * t74 * rb1 * dz - 0.3e1 * t74 * rb0 - 0.2e1 * t6 * t39 - t9 * t82 + 0.2e1 * t12 * t13 * x[N - 1] - t42 * 1.0 * u[0] + t42 * t43 * x[N] - t42 * t21 * u[0] + t42 * 0.0 * t55 * rb1;
  t98 = t22 * t22;
  dfdx0 = (t2 + 0.4e1 * t4 - 0.2e1 * t7 - 0.2e1 * t10 + 0.2e1 * t14) * t17 * t24;
  dfdx1 = t50 * t17 * t24 - t96 * t17 * t20 / t98 * 0.0;

  out[N-1] += dfdx0*adj[N];
  out[N  ] += dfdx1*adj[N];
}


void sysjacuadj(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *adj,
		        typeRNum *u, typeRNum *pSys)
{
  typeRNum q0  = pSys[0];
  typeRNum q1  = pSys[1];
  typeRNum nu  = pSys[2];
  typeRNum r0  = pSys[3];
  typeRNum r1  = pSys[4];
  typeRNum rb0 = pSys[5];
  typeRNum rb1 = pSys[6];
  typeInt N = NX-1;

  typeRNum dz = 1.0/N;

  typeRNum dfdu = (0.2e1 * q0 * dz + 0.4e1 * q1 * x[N] * dz - 0.2e1 * q1 * dz * x[N - 1] - nu * dz * dz * 1.0 - nu * dz * dz * 0.0 * x[N]) / rb0 * pow(dz, -0.2e1) / (1.0 + 0.0 * x[N]);

  out[0] = dfdu*adj[N];
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
