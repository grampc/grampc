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
 * File: heun2.c
 * Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Date: February 2014
 * Version: v1.0
 *
 * Heun integrator for GRAMPC.
 *
 */


#include "heun2.h"

void intsysHeun(typeRNum *y, typeInt pInt, typeInt Nint, typeRNum *t, typeRNum *x,
    typeRNum *u,typeGRAMPC *grampc,
    void (*pfct)(typeRNum *,typeRNum *,typeRNum *,typeRNum *,typeRNum *,typeGRAMPC *))
{
  typeInt i,j;
  typeRNum h, h2;
  typeRNum *s1  = grampc->rws->rwsIntegration;
  typeRNum *ys1 = s1  + grampc->param->Nx;
  typeRNum *s2  = ys1 + grampc->param->Nx;

  for (j = 0; j <= Nint-2; j++) {
    h  = t[pInt] - t[0];
    h2 = 0.5*h;

    /* s1 */
    (*pfct)(s1,y,t,x,u,grampc);
    /* s2 */
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      ys1[i] = y[i] + h * s1[i];
    }
    (*pfct)(s2,ys1,t+pInt,x+pInt*grampc->param->Nx,u+pInt*grampc->param->Nu,grampc);
    /* xnew */
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      y[i+pInt*grampc->param->Nx] = y[i] + h2 * (s1[i] + s2[i]);
    }

    /* next pointers */
    t += pInt;
    x += pInt*grampc->param->Nx;
    u += pInt*grampc->param->Nu;
    y += pInt*grampc->param->Nx;
  }
}
