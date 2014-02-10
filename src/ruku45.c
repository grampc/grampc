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
 * File: ruku45.c
 * Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Date: February 2014
 * Version: v1.0
 *
 * Fourth order Runge-Kutta integrator for GRAMPC.
 *
 */


#include "grampc_run.h"
#include "ruku45.h"

void intsysRuKu45(typeRNum *y, typeInt pInt, typeInt Nint, typeRNum *t, typeRNum *x,
    typeRNum *u,typeGRAMPC *grampc,
    void (*pfct)(typeRNum *,typeRNum *,typeRNum *,typeRNum *,typeRNum *,typeGRAMPC *))
{
  typeInt i;
  typeInt Nres = Nint;
  typeRNum h, hnext;
  typeRNum *s1  = grampc->rws->rwsIntegration;
  typeRNum *ys1 = s1  + grampc->param->Nx;
  typeRNum *s2  = ys1 + grampc->param->Nx;
  typeRNum *ys2 = s2  + grampc->param->Nx;
  typeRNum *s3  = ys2 + grampc->param->Nx;
  typeRNum *ys3 = s3  + grampc->param->Nx;
  typeRNum *s4  = ys3 + grampc->param->Nx;
  typeRNum *ys4 = s4  + grampc->param->Nx;
  typeRNum *s5  = ys4 + grampc->param->Nx;
  typeRNum *ys5 = s5  + grampc->param->Nx;
  typeRNum *s6  = ys5 + grampc->param->Nx;
  typeRNum *ys6 = s6  + grampc->param->Nx;
  typeRNum *s7  = ys6 + grampc->param->Nx;
  typeRNum *ys7 = s7  + grampc->param->Nx;

  typeRNum *tnow = t;
  typeRNum *xnow = x;
  typeRNum *unow = u;
  typeRNum *ynow = y;
  typeRNum tt    = t[pInt*(Nint-1)];
  typeRNum *tvec;
  typeRNum *uvec;
  typeRNum *zvec;

  /* intermediate steps */
  typeRNum t1;
  typeRNum t2;
  typeRNum *y1   = ys7  + grampc->param->Nx;
  typeRNum *y2   = y1   + grampc->param->Nx;
  typeRNum *uakt = y2   + grampc->param->Nx;
  typeRNum *xakt = uakt + grampc->param->Nu;
  typeRNum tdummy;

  /* Dormand-Prince-formula parameters */
  typeRNum a21 = 1.0/5.0,
      a31 = 3.0/40.0,
      a32 = 9.0/40.0,
      a41 = 44.0/45.0,
      a42 = -56.0/15.0,
      a43 = 32.0/9.0,
      a51 = 19372.0/6561.0,
      a52 = -25360.0/2187.0,
      a53 = 64448.0/6561.0,
      a54 = -212.0/729.0,
      a61 = 9017.0/3168.0,
      a62 = -355.0/33.0,
      a63 = 46732.0/5247.0,
      a64 = 49.0/176.0,
      a65 = -5103.0/18656.0;
  typeRNum b1 = 35.0/384.0,
      b2 = 0.0,
      b3 = 500.0/1113.0,
      b4 = 125.0/192.0,
      b5 = -2187.0/6784.0,
      b6 = 11.0/84.0,
      b7 = 0.0;
  typeRNum a71 = 35.0/384.0,
      a72 = 0.0,
      a73 = 500.0/1113.0,
      a74 = 125.0/192.0,
      a75 = -2187.0/6784.0,
      a76 = 11.0/84.0;
  typeRNum b1s = 5179.0/57600.0,
      b2s = 0.0,
      b3s = 7571.0/16695.0,
      b4s = 393.0/640.0,
      b5s = -92097.0/339200.0,
      b6s = 187.0/2100.0,
      b7s = 1.0/40.0;
  typeRNum db1 = b1-b1s,
      db2 = b2-b2s,
      db3 = b3-b3s,
      db4 = b4-b4s,
      db5 = b5-b5s,
      db6 = b6-b6s,
      db7 = b7-b7s;
  typeRNum c2 = 1.0/5.0,
      c3 = 3.0/10.0,
      c4 = 4.0/5.0,
      c5 = 8.0/9.0,
      c6 = 1.0,
      c7 = 1.0;
  typeRNum d1 = -12715105075.0/11282082432.0,
      d3 = 87487479700.0/32700410799.0,
      d4 = -10690763975.0/1880347072.0,
      d5 = 701980252875.0/199316789632.0,
      d6 = -1453857185.0/822651844.0,
      d7 = 69997945.0/29380423.0;

  typeInt reject    = 0;
  typeRNum minscale = 1.0/5.0;
  typeRNum maxscale = 10.0;
  typeRNum errold   = 0.0001;
  typeRNum bb       = 0.08;
  typeRNum k        = 5.0; 						/* Verfahrensordnung */
  typeRNum aa       = 1/k-bb*0.75;
  typeRNum safe     = 0.9;
  typeInt ende      = 0;
  typeRNum scale;
  typeRNum sk;
  typeRNum err;
  typeRNum errPow;

  typeRNum theta;
  typeRNum interp1, interp2, interp3, interp4, interp5;

  if (pInt == FWINT) {

  }
  else {

  }

  if (pInt == FWINT) {
    tvec = t;
    uvec = u;
    zvec = x;
  }
  else {
    tvec = t - (Nint-1);
    uvec = u - grampc->param->Nu*(Nint-1);
    zvec = x - grampc->param->Nx*(Nint-1);
  }

  /* initialization */
  t1 = tnow[0];
  for (i = 0; i <= grampc->param->Nx-1; i++) {
    y1[i] = ynow[i];
  }
  for (i = 0; i <= grampc->param->Nu-1; i++) {
    uakt[i] = unow[i];
  }
  if (y != x) {
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      xakt[i] = xnow[i];
    }
  }
  h = tnow[pInt] - tnow[0];

  /* Stage 1 */
  (*pfct)(s1,ynow,tnow,xnow,unow,grampc);

  while (pInt*(tt-tnow[pInt]) >= 0 && ende == 0) {
    /* START WHILE LOOP *******************************************************/

    /* Stage 2 */
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      ys1[i] = y1[i] + h*(a21*s1[i]);
    }
    tdummy = t1 + c2*h;
    interplin(uakt,tvec,uvec,tdummy,grampc->param->Nu,Nres,pInt);
    if (y != x) {
      interplin(xakt,tvec,zvec,tdummy,grampc->param->Nx,Nres,pInt);
    }
    (*pfct)(s2,ys1,&tdummy,xakt,uakt,grampc);

    /* Stage 3 */
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      ys2[i] = y1[i] + h*(a31*s1[i] + a32*s2[i]);
    }
    tdummy = t1 + c3*h;
    interplin(uakt,tvec,uvec,tdummy,grampc->param->Nu,Nres,pInt);
    if (y != x) {
      interplin(xakt,tvec,zvec,tdummy,grampc->param->Nx,Nres,pInt);
    }
    (*pfct)(s3,ys2,&tdummy,xakt,uakt,grampc);

    /* Stage 4 */
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      ys3[i] = y1[i] + h*(a41*s1[i] + a42*s2[i] + a43*s3[i]);
    }
    tdummy = t1 + c4*h;
    interplin(uakt,tvec,uvec,tdummy,grampc->param->Nu,Nres,pInt);
    if (y != x) {
      interplin(xakt,tvec,zvec,tdummy,grampc->param->Nx,Nres,pInt);
    }
    (*pfct)(s4,ys3,&tdummy,xakt,uakt,grampc);

    /* Stage 5 */
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      ys4[i] = y1[i] + h*(a51*s1[i] + a52*s2[i] + a53*s3[i] + a54*s4[i]);
    }
    tdummy = t1 + c5*h;
    interplin(uakt,tvec,uvec,tdummy,grampc->param->Nu,Nres,pInt);
    if (y != x) {
      interplin(xakt,tvec,zvec,tdummy,grampc->param->Nx,Nres,pInt);
    }
    (*pfct)(s5,ys4,&tdummy,xakt,uakt,grampc);

    /* Stage 6 */
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      ys5[i] = y1[i] + h*(a61*s1[i] + a62*s2[i] + a63*s3[i] + a64*s4[i] + a65*s5[i]);
    }
    tdummy = t1 + c6*h;
    interplin(uakt,tvec,uvec,tdummy,grampc->param->Nu,Nres,pInt);
    if (y != x) {
      interplin(xakt,tvec,zvec,tdummy,grampc->param->Nx,Nres,pInt);
    }
    (*pfct)(s6,ys5,&tdummy,xakt,uakt,grampc);

    /* Stage 7 */
    for (i = 0;i <= grampc->param->Nx-1; i++) {
      ys6[i] = y1[i] + h*(a71*s1[i] + a72*s2[i] + a73*s3[i] + a74*s4[i] + a75*s5[i] + a76*s6[i]);
    }
    tdummy = t1 + c7*h;
    interplin(uakt,tvec,uvec,tdummy,grampc->param->Nu,Nres,pInt);
    if (y != x) {
      interplin(xakt,tvec,zvec,tdummy,grampc->param->Nx,Nres,pInt);
    }
    (*pfct)(s7,ys6,&tdummy,xakt,uakt,grampc);

    /* Error */
    t2  = t1 + h;
    err = 0.0;

    /* HIER EVENTUELL NOCHMAL PR�FEN */
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      y2[i]  = y1[i] + h*(b1*s1[i] + b2*s2[i] + b3*s3[i] + b4*s4[i] + b5*s5[i] + b6*s6[i]);
      sk     = grampc->opt->IntegratorAbsTol + grampc->opt->IntegratorRelTol*MAX(fabs(y2[i]),fabs(y1[i]));
      errPow = h*(db1*s1[i] + db2*s2[i] + db3*s3[i] + db4*s4[i] + db5*s5[i] + db6*s6[i] + db7*s7[i])/sk;
      err   += errPow*errPow;
    }
    err = sqrt(err/(grampc->param->Nx));

    /* Schrittweitensteuerung */
    if (err<=1.0) {
      if (err==0.0) {
        scale=maxscale;
      }
      else {
        /* POW BEFEHLE DIE MAN RAUSNEHMEN K�NNTE ?? */
        scale = safe*pow(err,-aa)*pow(errold,bb);
        if (scale < minscale) {
          scale = minscale;
        }
        if (scale > maxscale) {
          scale = maxscale;
        }
      }
      if (reject == 1) {
        hnext = h*MIN(scale,1);
      }
      else {
        hnext = h*scale;
      }
      errold = MAX(err,0.0001);
      reject = 0;
    }
    else {
      scale  = MAX(safe*pow(err,-aa),minscale);
      h      = h*scale;
      reject = 1;
    }
    if (reject == 0) {
      while (pInt*(t2-tnow[pInt]) >= 0 && ende == 0) {
        theta = (tnow[pInt]-t1)/(t2-t1);
        for (i = 0; i <= grampc->param->Nx-1; i++)	{
          interp1 		     = y1[i];
          interp2 		     = y2[i] - y1[i];
          interp3 		     = h*s1[i] - interp2;
          interp4 		     = interp2 - h*s7[i] - interp3;
          interp5              = h*(d1*s1[i] + d3*s3[i] + d4*s4[i] + d5*s5[i] + d6*s6[i] + d7*s7[i]);
          ynow[i+pInt*grampc->param->Nx] = interp1 + theta*(interp2 + (1.0 - theta)*(interp3 + theta*(interp4 + (1 - theta)*interp5)));
        }
        if (tnow[pInt] == tt) {
          ende = 1;
        }
        else {
          Nres -= 1;
          tnow += pInt;
          xnow += pInt*grampc->param->Nx;
          unow += pInt*grampc->param->Nu;
          ynow += pInt*grampc->param->Nx;
          if (pInt == FWINT) {
            tvec += 1;
            zvec += grampc->param->Nx;
            uvec += grampc->param->Nu;
          }
        }
      }
      h = hnext;
      /* Umspeicherung der lokalen Zeit- und Zustandsvariablen */
      t1 = t2;
      for (i = 0; i <= grampc->param->Nx-1; i++) {
        y1[i] = y2[i];
      }
      for (i = 0; i <= grampc->param->Nx-1; i++) {
        s1[i] = s7[i]; /* Verwendung der vorherigen letzten Stage fuer den kommenden Schritt */
      }
    }
    /* END WHILE LOOP *********************************************************/
  }
}
