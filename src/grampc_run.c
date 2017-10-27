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
 * File: grampc_run.c
 * Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Date: February 2014
 * Version: v1.0
 *
 * MPC running file for GRAMPC.
 *
 */


#include "grampc_run.h"
#include "probfct.h"
#include "grampc_mess.h"

void grampc_run(typeGRAMPC *grampc)
{
  typeInt i, j, k, igrad,shift,shiftu;
  typeInt NLS_adapt = 2*(NLS + 1);
  typeRNum *ui;
  typeRNum *t          = grampc->rws->t;
  typeRNum *x          = grampc->rws->x;
  typeRNum *adj        = grampc->rws->adj;
  typeRNum *u          = grampc->rws->u;
  typeRNum *dHdu       = grampc->rws->dHdu;
  typeRNum *lsAdapt    = grampc->rws->lsAdapt + grampc->opt->MaxIter*NLS_adapt;
  typeRNum *uls        = grampc->rws->uls;
  typeRNum *lsExplicit = grampc->rws->lsExplicit;
  typeRNum *uprev      = grampc->rws->uprev;
  typeRNum *dHduprev   = grampc->rws->dHduprev;
  typeRNum *x_unscaled = NULL;
  typeRNum interpfact;

  /* function pointers for system and cost integration */
  void (*pIntSys)(typeRNum *, typeInt , typeInt, typeRNum *, typeRNum *, typeRNum *,
		  typeGRAMPC *, void (*)(typeRNum *,typeRNum *,typeRNum *,typeRNum *,typeRNum *,typeGRAMPC *));
  void (*pIntCost)(typeRNum *, typeRNum *, typeRNum *, typeRNum *,
		   typeGRAMPC *);

  /* check init- and setpoints */
  if (grampc->param->u0 == NULL) {
    grampc_error(U0_NOT_DEFINED);
  }
  if (grampc->param->xk == NULL) {
    grampc_error(XK_NOT_DEFINED);
  }
  if (grampc->param->udes == NULL) {
    grampc_error(UDES_NOT_DEFINED);
  }
  if (grampc->param->xdes == NULL) {
    grampc_error(XDES_NOT_DEFINED);
  }
  if (grampc->param->dt <= 0.0) {
    grampc_error(DT_NOT_DEFINED);
  }
  if (grampc->param->Thor <= 0.0) {
    grampc_error(THOR_NOT_DEFINED);
  }
  
  /* check if additional parameters are filled appropriately */
  if ((grampc->param->NpCost > 0) && (grampc->param->pCost == NULL)) {
    grampc_error(PCOST_NULL);
  }
  if ((grampc->param->NpSys > 0) && (grampc->param->pSys == NULL)) {
    grampc_error(PSYS_NULL);
  }

  /* integrator for system integration */
  if (!strncmp(grampc->opt->Integrator,"euler",VALUE_INTEGRATOR)) {
    pIntSys = &intsysEuler;
  }
  else if (!strncmp(grampc->opt->Integrator,"modeuler",VALUE_INTEGRATOR)) {
    pIntSys = &intsysModEuler;
  }
  else if (!strncmp(grampc->opt->Integrator,"heun",VALUE_INTEGRATOR)) {
    pIntSys = &intsysHeun;
  }
  else {
    pIntSys = &intsysRuKu45;
  }

  /* integrator for cost integration */
  if (!strncmp(grampc->opt->CostIntegrator,"trapezodial",VALUE_COSTINTMETHOD)) {
    pIntCost = &intCostTrapezodial;
  }
  else {
    pIntCost = &intCostSimpson;
  }

  /* initial condition */
  if (!strncmp(grampc->opt->ScaleProblem,"on",VALUE_ONOFF)) {
    x_unscaled = grampc->rws->rwsScale;
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      x[i] = (grampc->param->xk[i]-grampc->param->xOffset[i])/grampc->param->xScale[i];
    }
  }
  else {
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      x[i] = grampc->param->xk[i];
    }
  }

  /* shift of input trajectory */
  if (!strncmp(grampc->opt->ShiftControl,"on",VALUE_ONOFF)) {
    ui = u;

	/* compute how far the sampling points must be shifted */
	shift = (typeInt)(grampc->param->dt / (t[1] - t[0]));
	shiftu = shift*grampc->param->Nu;
	interpfact = (grampc->param->dt / (t[1] - t[0])) - shift;

	if (shift >= grampc->param->Nhor){
		grampc_error("Horizon too short for the current sampling time.");
	}

	/* Interpolation between the sampling points */
    for (i = 0; i <= grampc->param->Nhor-2-shift; i++) {
      for (j = 0; j <= grampc->param->Nu-1; j++) {
		  ui[j] = ui[shiftu + j] + (ui[shiftu + j + grampc->param->Nu] - ui[shiftu + j]) *interpfact;
      }
      ui += grampc->param->Nu;
    }

    /* next element: extrapolation */
	if (shift == grampc->param->Nhor - 1){
		for (j = 0; j <= grampc->param->Nu - 1; j++) {
			ui[j] = ui[shiftu + j] + (ui[shiftu + j] - ui[shiftu + j - grampc->param->Nu])* interpfact;
		}
	}
	else{
		for (j = 0; j <= grampc->param->Nu - 1; j++) {
			ui[j] = ui[shiftu + j] + (ui[shiftu + j] - ui[j - grampc->param->Nu]) / (1 / interpfact - 1);
		}
    }
	ui += grampc->param->Nu;
	i++;

	/* if there are elements left hold last value */
	for (; i <= grampc->param->Nhor - 1; i++){
		for (j = 0; j <= grampc->param->Nu - 1; j++) {
			ui[j] = ui[j - grampc->param->Nu];
		}
		ui += grampc->param->Nu;
	}
  }

  /* LOOP OVER NO. OF GRADIENT STEPS *********************************************/
  for (igrad = 0; igrad <= grampc->opt->MaxIter-1; igrad++) {
    /* forward integration of system */
    (*pIntSys)(x,FWINT,grampc->param->Nhor,t,x,u,grampc,&Wsys);
    /* final conditions for costates */
    if (!strncmp(grampc->opt->FinalCost,"on",VALUE_ONOFF)) {
      if (!strncmp(grampc->opt->ScaleProblem,"on",VALUE_ONOFF)) {
        for (i = 0; i <= grampc->param->Nx-1; i++) {
          x_unscaled[i] = x[grampc->param->Nx*(grampc->param->Nhor-1)+i]*grampc->param->xScale[i] + grampc->param->xOffset[i];
        }
        fcostjacx(adj+grampc->param->Nx*(grampc->param->Nhor-1),t[grampc->param->Nhor-1],x_unscaled,grampc->param->xdes,grampc->param->pCost);
        for (i = 0; i <= grampc->param->Nx-1; i++) {
          adj[grampc->param->Nx*(grampc->param->Nhor-1)+i] = grampc->param->xScale[i]*adj[grampc->param->Nx*(grampc->param->Nhor-1)+i];
        }
      }
      else {
        fcostjacx(adj+grampc->param->Nx*(grampc->param->Nhor-1),t[grampc->param->Nhor-1],x+grampc->param->Nx*(grampc->param->Nhor-1),grampc->param->xdes,grampc->param->pCost);
      }
    }
    /* integration of adjoint system in reverse time */
    (*pIntSys)(adj+grampc->param->Nx*(grampc->param->Nhor-1),BWINT,grampc->param->Nhor,t+(grampc->param->Nhor-1),x+grampc->param->Nx*(grampc->param->Nhor-1),u+grampc->param->Nu*(grampc->param->Nhor-1),grampc,&Wadjsys);
    /* residual in dHdu */
    dHdufct(grampc);
    /* adaptive line search */
    if (!strncmp(grampc->opt->LineSearchType,"adaptive",VALUE_LSTYPE)) {
      /* line search: adaption of interval */
      if (lsAdapt[NLS]>=lsAdapt[0]+(1-grampc->opt->LineSearchIntervalTol)*(lsAdapt[NLS-1]-lsAdapt[0]) && lsAdapt[NLS-1]<=grampc->opt->LineSearchMax) {
	for (i = 0; i <= NLS-1; i++) {
	  lsAdapt[i] = lsAdapt[i]*grampc->opt->LineSearchAdaptFactor;
	}
      }
      else if (lsAdapt[NLS]<=lsAdapt[0]+grampc->opt->LineSearchIntervalTol*(lsAdapt[NLS-1]-lsAdapt[0]) && lsAdapt[0]>=grampc->opt->LineSearchMin) {
	for (i = 0; i <= NLS-1; i++) {
	  lsAdapt[i] = lsAdapt[i]/grampc->opt->LineSearchAdaptFactor;
	}
      }      
      for (i = 0; i <= NLS-1; i++) {
        for (j = 0; j <= grampc->param->Nhor-1; j++) {
          for (k = 0; k <= grampc->param->Nu-1; k++) {
            uls[j*grampc->param->Nu+k] = u[j*grampc->param->Nu+k]-lsAdapt[i]*dHdu[j*grampc->param->Nu+k];
          }
        }
        inputproj(uls,grampc);
        /* grampc->rws->u = uls; */
        (*pIntSys)(x,FWINT,grampc->param->Nhor,t,x,uls,grampc,&Wsys);
        /* grampc->rws->u = u; */
        (*pIntCost)(lsAdapt+i+NLS+1,t,x,uls,grampc);
      }
      /* curve fitting */
      lsearch_fit2(lsAdapt+NLS,lsAdapt+2*NLS+1,lsAdapt,lsAdapt+NLS+1);
      /* Next trajectory */
      for (j = 0; j <= grampc->param->Nhor-1; j++) {
        for (k = 0; k <= grampc->param->Nu-1; k++) {
          u[j*grampc->param->Nu+k] = u[j*grampc->param->Nu+k]-lsAdapt[NLS]*dHdu[j*grampc->param->Nu+k];
        }
      }
      for (i = 0; i <= 2*(NLS+1)-1; i++) {
	grampc->rws->lsAdapt[i+igrad*NLS_adapt] = lsAdapt[i];
      }
    }
    /* explicit line search */
    else {
      lsExplicit[0] = 0;
      lsExplicit[1] = 0;
      for (j = 0; j <= grampc->param->Nhor-1; j++) {
        for (k = 0; k <= grampc->param->Nu-1; k++) {
          if (!strncmp(grampc->opt->LineSearchType,"explicit1",VALUE_LSTYPE)) {
            /* Formula 1 */
            lsExplicit[0] = lsExplicit[0] + (u[j*grampc->param->Nu+k]-uprev[j*grampc->param->Nu+k])*(u[j*grampc->param->Nu+k]-uprev[j*grampc->param->Nu+k]);
            lsExplicit[1] = lsExplicit[1] + (u[j*grampc->param->Nu+k]-uprev[j*grampc->param->Nu+k])*(dHdu[j*grampc->param->Nu+k]-dHduprev[j*grampc->param->Nu+k]);
          }
          else {
            /* Formula #2 */
            lsExplicit[0] = lsExplicit[0] + (u[j*grampc->param->Nu+k]-uprev[j*grampc->param->Nu+k])*(dHdu[j*grampc->param->Nu+k]-dHduprev[j*grampc->param->Nu+k]);
            lsExplicit[1] = lsExplicit[1] + (dHdu[j*grampc->param->Nu+k]-dHduprev[j*grampc->param->Nu+k])*(dHdu[j*grampc->param->Nu+k]-dHduprev[j*grampc->param->Nu+k]);
          }
        }
      }
      if (lsExplicit[0] > 0 && lsExplicit[1] > 0) {
        lsExplicit[2] = lsExplicit[0]/lsExplicit[1];
      }
      else {
        lsExplicit[2] = grampc->opt->LineSearchInit;
      }
      if (lsExplicit[2] > grampc->opt->LineSearchMax) {
        lsExplicit[2] = grampc->opt->LineSearchMax;
      }
      /* Next trajectory */
      for (j = 0; j <= grampc->param->Nhor-1; j++) {
        for (k = 0; k <= grampc->param->Nu-1; k++) {
          uprev[j*grampc->param->Nu+k]    = u[j*grampc->param->Nu+k];
          dHduprev[j*grampc->param->Nu+k] = dHdu[j*grampc->param->Nu+k];
          u[j*grampc->param->Nu+k]        = u[j*grampc->param->Nu+k]-lsExplicit[2]*dHdu[j*grampc->param->Nu+k];
        }
      }
    }
    inputproj(u,grampc);
  }
  /* END GRADIENT LOOP ***********************************************************/
  if (!strncmp(grampc->opt->ScaleProblem,"on",VALUE_ONOFF)) {
    for (i = 0; i <= grampc->param->Nu-1; i++) {
      grampc->sol->unext[i] = u[i]*grampc->param->uScale[i] + grampc->param->uOffset[i];
    }
  }
  else {
    for (i = 0; i <= grampc->param->Nu-1; i++) {
      grampc->sol->unext[i] = u[i];
    }
  }

  /* calculation of cost value */
  (*pIntSys)(x,FWINT,grampc->param->Nhor,t,x,u,grampc,&Wsys);
  (*pIntCost)(grampc->sol->J,t,x,u,grampc);
  
  /* calculation of xnext */
  i = 0;
  while (t[i] < grampc->param->dt) {
      i++;
  }
  interplin(grampc->sol->xnext,t,x,grampc->param->dt,grampc->param->Nx,i+1,1);
  if (!strncmp(grampc->opt->ScaleProblem,"on",VALUE_ONOFF)) {
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      grampc->sol->xnext[i] = grampc->sol->xnext[i]*grampc->param->xScale[i] + grampc->param->xOffset[i];
    }
  }
}


void dHdufct(typeGRAMPC *grampc)
{
  typeInt i, j;

  typeRNum *x_unscaled    = NULL;
  typeRNum *u_unscaled    = NULL;
  typeRNum *xadj_unscaled = NULL;

  typeRNum *dLdu = grampc->rws->rwsGradient;
  typeRNum *dfdu = dLdu + grampc->param->Nu;
  typeRNum *s    = dfdu + grampc->param->Nu*grampc->param->Nx;
  if (!strncmp(grampc->opt->ScaleProblem,"on",VALUE_ONOFF)) {
    x_unscaled   = grampc->rws->rwsScale;
    u_unscaled   = x_unscaled + grampc->param->Nx;
    xadj_unscaled = u_unscaled + grampc->param->Nu;
    for (i = 0; i <= grampc->param->Nhor-1; i++) {
      for (j = 0; j <= grampc->param->Nx-1; j++) {
        x_unscaled[j]   = grampc->rws->x[i*grampc->param->Nx+j]*grampc->param->xScale[j] + grampc->param->xOffset[j];
        xadj_unscaled[j] = grampc->rws->adj[i*grampc->param->Nx+j]/grampc->param->xScale[j];;
      }
      for (j = 0; j <= grampc->param->Nu-1; j++) {
        u_unscaled[j] = grampc->rws->u[i*grampc->param->Nu+j]*grampc->param->uScale[j] + grampc->param->uOffset[j];
      }
      if (!strncmp(grampc->opt->IntegralCost,"on",VALUE_ONOFF)) {
        icostjacu(dLdu,grampc->rws->t[i],x_unscaled,u_unscaled,grampc->param->xdes,grampc->param->udes,grampc->param->pCost);
      }
      if (!strncmp(grampc->opt->JacobianU,"sysjacuadj",VALUE_JACOBIANU)) {
        sysjacuadj(s,grampc->rws->t[i]+grampc->param->tk,x_unscaled,xadj_unscaled,u_unscaled,grampc->param->pSys);
      }
      else {
        sysjacu(dfdu,grampc->rws->t[i]+grampc->param->tk,x_unscaled,u_unscaled,grampc->param->pSys);
        MatMult(s,xadj_unscaled,dfdu,1,grampc->param->Nx,grampc->param->Nu);
      }
      for (j = 0; j <= grampc->param->Nu-1; j++) {
        grampc->rws->dHdu[i*grampc->param->Nu+j] = grampc->param->uScale[j]*dLdu[j] + s[j]*grampc->param->uScale[j];
      }
    }
  }
  else {
    for (i = 0; i <= grampc->param->Nhor-1; i++) {
      if (!strncmp(grampc->opt->IntegralCost,"on",VALUE_ONOFF)) {
        icostjacu(dLdu,grampc->rws->t[i],grampc->rws->x+i*grampc->param->Nx,grampc->rws->u+i*grampc->param->Nu,grampc->param->xdes,grampc->param->udes,grampc->param->pCost);
      }
      if (!strncmp(grampc->opt->JacobianU,"sysjacuadj",VALUE_JACOBIANU)) {
        sysjacuadj(s,grampc->rws->t[i]+grampc->param->tk,grampc->rws->x+i*grampc->param->Nx,grampc->rws->adj+i*grampc->param->Nx,grampc->rws->u+i*grampc->param->Nu,grampc->param->pSys);
      }
      else {
        sysjacu(dfdu,grampc->rws->t[i]+grampc->param->tk,grampc->rws->x+i*grampc->param->Nx,grampc->rws->u+i*grampc->param->Nu,grampc->param->pSys);
        MatMult(s,grampc->rws->adj+i*grampc->param->Nx,dfdu,1,grampc->param->Nx,grampc->param->Nu);
      }
      for (j = 0; j <= grampc->param->Nu-1; j++) {
        grampc->rws->dHdu[i*grampc->param->Nu+j] = dLdu[j] + s[j];
      }
    }
  }
}

void inputproj(typeRNum *u, typeGRAMPC *grampc)
{
  typeInt i, j;

  if (!strncmp(grampc->opt->ScaleProblem,"on",VALUE_ONOFF)) {
    for (i = 0; i <= grampc->param->Nhor-1; i++) {
      for (j = 0; j <= grampc->param->Nu-1; j++) {
        /* lower bound */
        if (u[j+i*grampc->param->Nu] < (grampc->param->umin[j]-grampc->param->uOffset[j])/grampc->param->uScale[j]) {
          u[j+i*grampc->param->Nu] = (grampc->param->umin[j]-grampc->param->uOffset[j])/grampc->param->uScale[j];
        }
        /* upper bound */
        else if (u[j+i*grampc->param->Nu] > (grampc->param->umax[j]-grampc->param->uOffset[j])/grampc->param->uScale[j]) {
          u[j+i*grampc->param->Nu] = (grampc->param->umax[j]-grampc->param->uOffset[j])/grampc->param->uScale[j];
        }
      }
    }
  }
  else {
    for (i = 0; i <= grampc->param->Nhor-1; i++) {
      for (j = 0; j <= grampc->param->Nu-1; j++) {
        /* lower bound */
        if (u[j+i*grampc->param->Nu] < grampc->param->umin[j]) {
          u[j+i*grampc->param->Nu] = grampc->param->umin[j];
        }
        /* upper bound */
        else if (u[j+i*grampc->param->Nu] > grampc->param->umax[j]) {
          u[j+i*grampc->param->Nu] = grampc->param->umax[j];
        }
      }
    }
  }
}

void lsearch_fit2(typeRNum *kfit, typeRNum *Jfit, typeRNum *k, typeRNum *J)
{
  typeRNum k02 = k[0]*k[0];
  typeRNum k12 = k[1]*k[1];
  typeRNum k22 = k[2]*k[2];
  typeRNum a0 = (J[2]*k[0]*(k[0]-k[1])*k[1]+k[2]*(J[0]*k[1]*(k[1]-k[2])+\
						  J[1]*k[0]*(-k[0]+k[2])))/((k[0]-k[1])*(k[0]-k[2])*(k[1]-k[2]));
  typeRNum a1 = (J[2]*(-k02+k12)+J[1]*(k02-k22)+J[0]*\
		 (-k12+k22))/((k[0]-k[1])*(k[0]-k[2])*(k[1]-k[2]));
  typeRNum a2 = ((J[0]-J[2])/(k[0]-k[2])+(-J[1]+J[2])/(k[1]-k[2]))/(k[0]-k[1]);

  /* minimum? */
  if (a2 >= aEPS) {
    kfit[0] = -a1/(2*a2);
    Jfit[0] = a0 + a1*kfit[0] + a2*kfit[0]*kfit[0];
  }
  /* smallest J */
  if (a2 < aEPS || kfit[0] < k[0] || kfit[0] > k[2]) {
    if (J[0] <= J[1]-JEPS && J[0] <= J[2]-JEPS) {
      kfit[0] = k[0];
      Jfit[0] = J[0];
    }
    else if (J[2] <= J[0]-JEPS && J[2] <= J[1]-JEPS) {
      kfit[0] = k[2];
      Jfit[0] = J[2];
    }
    else {
      kfit[0] = k[1];
      Jfit[0] = J[1];
    }
  }
}

void interplin(typeRNum *varint, typeRNum *tvec, typeRNum *varvec, typeRNum tint,
	       typeInt Nvar, typeInt Nvec, typeInt searchdir)
{
  /* option to determine position ioff in time vector such that tvec[ioff] <= tint <= tvec[ioff+1]
   * searchdir= 1: search in forward direction starting at ioff=0									 *
   * searchdir=-1: search in backward direction starting at ioff=Nvec-1
   */
  typeInt i;
  typeInt ioff;
  typeRNum dtratio;
  typeRNum *var0, *var1;

  if (tint < tvec[0]) {
    for (i = 0; i <= Nvar-1; i++) {
      varint[i] = varvec[i];
    }
  }
  else if (tint > tvec[Nvec-1]) {
    var0 = varvec + (Nvec-1)*Nvar;
    for (i = 0; i <= Nvar-1; i++)
      varint[i] = var0[i];
  }
  else {
    if (searchdir == 1) {
      ioff = 0;
      while (tvec[ioff] < tint) {
        ioff+=1;
      }
      ioff-= 1;
    }
    else {
      ioff = Nvec-2;
      while (tvec[ioff] > tint) {
        ioff-=1;
      }
    }
    dtratio = (tint - tvec[ioff])/(tvec[ioff+1] - tvec[ioff]);
    var0 = varvec + ioff*Nvar;
    var1 = var0 + Nvar;
    for (i = 0; i <= Nvar-1; i++) {
      varint[i] = var0[i] + dtratio*(var1[i] - var0[i]);
    }
  }
}

void MatAdd(typeRNum *C, typeRNum *A, typeRNum *B, typeInt n1, typeInt n2)
{
  /* matrix summation C = A+B *
   * A,B,C: (n1 x n2)         */
  typeInt i;
  for (i = 0; i <= n1*n2-1; i++) {
    C[i] = A[i] + B[i];
  }
}

void MatMult(typeRNum *C, typeRNum *A, typeRNum *B, typeInt n1, typeInt n2, typeInt n3)
{
  /* matrix multiplication C = A*B			   *
   * A: (n1 x n2)    B: (n2 x n3)   C: (n1 x n3) */
  typeInt i,j,k;
  typeRNum sigma;
  for (i = 0; i <= n1-1; i++) {
    for (j = 0; j <= n3-1; j++) {
      sigma = 0;
      for (k = 0; k <= n2-1; k++) {
        sigma += A[i*n2 + k] * B[k*n3 + j];
      }
      C[i*n3 + j] = sigma;
    }
  }
}

void minfct(typeRNum *amin, typeInt *amini, typeRNum *a, typeInt Na)
{
  typeInt i;
  amin[0] = a[0];
  amini[0] = 0;
  for (i = 1; i <= Na-1; i++) {
    if (a[i] < amin[0]) {
      amin[0]  = a[i];
      amini[0] = i;
    }
  }
}

void maxfct(typeRNum *amax, typeInt *amaxi, typeRNum *a, typeInt Na)
{
  typeInt i;
  amax[0] = a[0];
  amaxi[0] = 1;
  for (i = 1; i <= Na-1; i++) {
    if (a[i] > amax[0]) {
      amax[0]  = a[i];
      amaxi[0] = i;
    }
  }
}

void Wadjsys(typeRNum *s, typeRNum *adj, typeRNum *t, typeRNum *x, typeRNum *u, typeGRAMPC *grampc)
{
  typeRNum *x_   = x;
  typeRNum *u_   = u;
  typeRNum *adj_ = adj;
  typeRNum *dLdx = grampc->rws->rwsAdjIntegration;
  typeRNum *dfdx = dLdx + grampc->param->Nx;
  typeInt i;

  if (!strncmp(grampc->opt->ScaleProblem,"on",VALUE_ONOFF)) {
    x_   = grampc->rws->rwsScale;
    u_   = x_ + grampc->param->Nx;
    adj_ = u_ + grampc->param->Nu;
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      x_[i]    = grampc->param->xScale[i]*x[i] + grampc->param->xOffset[i];
      adj_[i] = adj[i]/grampc->param->xScale[i];
    }
    for (i = 0; i <= grampc->param->Nu-1; i++) {
      u_[i] = grampc->param->uScale[i]*u[i] + grampc->param->uOffset[i];
    }
  }

  if (!strncmp(grampc->opt->IntegralCost,"on",VALUE_ONOFF)) {
    icostjacx(dLdx,t[0],x_,u_,grampc->param->xdes,grampc->param->udes,grampc->param->pCost);
  }
  /*
   * No else branch is necessary since dLdx is initialized with zeros
   */
  if (!strncmp(grampc->opt->JacobianX,"sysjacxadj",VALUE_JACOBIANX)) {
    sysjacxadj(s,t[0]+grampc->param->tk,x_,adj_,u_,grampc->param->pSys);
  }
  else {
    sysjacx(dfdx,t[0]+grampc->param->tk,x_,u_,grampc->param->pSys);
    MatMult(s,adj_,dfdx,1,grampc->param->Nx,grampc->param->Nx);
  }
  if (!strncmp(grampc->opt->ScaleProblem,"on",VALUE_ONOFF)) {
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      s[i] = -dLdx[i]*grampc->param->xScale[i] - s[i]*grampc->param->xScale[i];
    }
  }
  else {
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      s[i] = -dLdx[i] - s[i];
    }
  }
}

void Wsys(typeRNum *s, typeRNum *x, typeRNum *t, typeRNum *dummy, typeRNum *u, typeGRAMPC *grampc)
{
  typeRNum *x_unscaled = NULL;
  typeRNum *u_unscaled = NULL;
  typeInt i;

  if (!strncmp(grampc->opt->ScaleProblem,"on",VALUE_ONOFF)) {
    x_unscaled    = grampc->rws->rwsScale;
    u_unscaled    = x_unscaled + grampc->param->Nx;
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      x_unscaled[i] = grampc->param->xScale[i]*x[i] + grampc->param->xOffset[i];
    }
    for (i = 0; i <= grampc->param->Nu-1; i++) {
      u_unscaled[i] = grampc->param->uScale[i]*u[i] + grampc->param->uOffset[i];
    }
    sysfct(s,t[0]+grampc->param->tk,x_unscaled,u_unscaled,grampc->param->pSys);
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      s[i] = s[i]/grampc->param->xScale[i];
    }
  }
  else {
    sysfct(s,t[0]+grampc->param->tk,x,u,grampc->param->pSys);
  }
}

void intCostTrapezodial(typeRNum *s, typeRNum *t, typeRNum *x, typeRNum *u, typeGRAMPC *grampc)
{
  typeInt i, j;
  typeRNum h;
  typeRNum *s1 = grampc->rws->rwsCostIntegration;
  
  typeRNum *x_ = NULL;
  typeRNum *u_ = NULL;
  
  s[0]  = 0;
  
  /* Integral Cost */
  if (!strncmp(grampc->opt->IntegralCost, "on", VALUE_ONOFF)) {
    for (i = 0; i <= grampc->param->Nhor - 1; i++) {
      
      /* Scaling */
      if (!strncmp(grampc->opt->ScaleProblem,"on",VALUE_ONOFF)) {
        x_ = grampc->rws->rwsScale;
        u_ = x_ + grampc->param->Nx;
        for (j = 0; j <= grampc->param->Nx-1; j++) {
          x_[j] = grampc->param->xScale[j] * x[j+i*grampc->param->Nx] + grampc->param->xOffset[j];
        }
        for (j = 0; j <= grampc->param->Nu-1; j++) {
          u_[j] = grampc->param->uScale[j] * u[j+i*grampc->param->Nu] + grampc->param->uOffset[j];
        }
      }
      else {
        x_ = x + i * grampc->param->Nx;
        u_ = u + i * grampc->param->Nu;
      }
      
      s1[0] = 0;
      icostfct(s1, t[i], x_, u_, grampc->param->xdes, grampc->param->udes, grampc->param->pCost);
      
      /* Integration */
      if(i == 0) {
        h = 0.5 * (t[i+1] - t[i]);
      }
      else if(i <= grampc->param->Nhor - 2) {
        h = 0.5 * (t[i+1] - t[i-1]);
      }
      else {
        h = 0.5 * (t[i] - t[i-1]);
      }
      s[0] = s[0] + h * s1[0];
    }
  }
  
  /* Final Cost */
  if (!strncmp(grampc->opt->FinalCost, "on", VALUE_ONOFF)) {
    /* Scaling */
    if (!strncmp(grampc->opt->ScaleProblem,"on",VALUE_ONOFF)) {
      x_ = grampc->rws->rwsScale;
      for (j = 0; j <= grampc->param->Nx-1; j++) {
        x_[j] = grampc->param->xScale[j] * x[j+(grampc->param->Nhor-1)*grampc->param->Nx] + grampc->param->xOffset[j];
      }
    }
    else {
      x_ = x + (grampc->param->Nhor - 1) * grampc->param->Nx;
    }
    fcostfct(s1, t[grampc->param->Nhor - 1], x_, grampc->param->xdes, grampc->param->pCost);
	s[0] = s[0] + s1[0];
  }
}

void intCostSimpson(typeRNum *s, typeRNum *t, typeRNum *x, typeRNum *u, typeGRAMPC *grampc)
{
  typeInt i, j;
  typeRNum h;
    
  typeRNum *s1 = grampc->rws->rwsCostIntegration;
  typeRNum *s2 = s1 + 1;
  typeRNum *ts = s2 + 1;
  typeRNum *xs = ts + 1;
  typeRNum *us = xs + grampc->param->Nx;  

  typeRNum *x1_ = NULL;
  typeRNum *u1_ = NULL;
  typeRNum *x2_ = NULL;
  typeRNum *u2_ = NULL;
      
  s[0]  = 0;
  
  /* Integral Cost */
  if (!strncmp(grampc->opt->IntegralCost,"on",VALUE_ONOFF)) {
    for (i = 0; i <= grampc->param->Nhor - 1; i++) {
      
      /* Scaling */
      if (!strncmp(grampc->opt->ScaleProblem,"on",VALUE_ONOFF)) {
        x1_ = grampc->rws->rwsScale;
        u1_ = x1_ + grampc->param->Nx;
        for (j = 0; j <= grampc->param->Nx-1; j++) {
          x1_[j] = grampc->param->xScale[j] * x[j+i*grampc->param->Nx] + grampc->param->xOffset[j];
        }
        for (j = 0; j <= grampc->param->Nu-1; j++) {
          u1_[j] = grampc->param->uScale[j] * u[j+i*grampc->param->Nu] + grampc->param->uOffset[j];
        }
      }
      else {
        x1_ = x + i * grampc->param->Nx;
        u1_ = u + i * grampc->param->Nu;
      }
      
      s1[0] = 0;
      icostfct(s1, t[i], x1_, u1_, grampc->param->xdes, grampc->param->udes, grampc->param->pCost);
      
      /* Integration */
      if(i == 0) {
        h = (t[i+1] - t[i]) / 6;
      }
      else if(i <= grampc->param->Nhor - 2) {
        h = (t[i+1] - t[i-1]) / 6;
      }
      else {
        h = (t[i] - t[i-1]) / 6;
      }
      s[0] = s[0] + h * s1[0];
      
      /* Cost for intermediate points */
      if (i <= grampc->param->Nhor - 2) {
        
        /* Scaling */
        if (!strncmp(grampc->opt->ScaleProblem,"on",VALUE_ONOFF)) {
          x2_ = u1_ + grampc->param->Nu;
          u2_ = x2_ + grampc->param->Nx;
          for (j = 0; j <= grampc->param->Nx-1; j++) {
            x2_[j] = grampc->param->xScale[j] * x[j+(i+1)*grampc->param->Nx] + grampc->param->xOffset[j];
          }
          for (j = 0; j <= grampc->param->Nu-1; j++) {
            u2_[j] = grampc->param->uScale[j] * u[j+(i+1)*grampc->param->Nu] + grampc->param->uOffset[j];
          }
        }
        else {
          x2_ = x + (i+1) * grampc->param->Nx;
          u2_ = u + (i+1) * grampc->param->Nu;
        }
        
        /* Interpolation */
        ts[0] = 0.5 * (t[i] + t[i+1]);
        for (j = 0; j <= grampc->param->Nx - 1; j++) {
          xs[j] = 0.5 * (x1_[j] + x2_[j]);
        }
        for (j = 0; j <= grampc->param->Nu - 1; j++) {
          us[j] = 0.5 * (u1_[j] + u2_[j]);
        }
        
        s2[0] = 0;
        icostfct(s2, ts[0], xs, us, grampc->param->xdes, grampc->param->udes, grampc->param->pCost);
        
        /* Integration */
        h = 4 * (t[i+1] - t[i]) / 6;
        s[0] = s[0] + h * s2[0];
      }
    }
  }
	/* Final Cost */
  if (!strncmp(grampc->opt->FinalCost, "on", VALUE_ONOFF)) {
    /* Scaling */
    if (!strncmp(grampc->opt->ScaleProblem,"on",VALUE_ONOFF)) {
      x1_ = grampc->rws->rwsScale;
      for (j = 0; j <= grampc->param->Nx-1; j++) {
        x1_[j] = grampc->param->xScale[j] * x[j+(grampc->param->Nhor-1)*grampc->param->Nx] + grampc->param->xOffset[j];
      }
    }
    else {
      x1_ = x + (grampc->param->Nhor - 1) * grampc->param->Nx;
    }
    fcostfct(s1, t[grampc->param->Nhor - 1], x1_, grampc->param->xdes, grampc->param->pCost);
    s[0] = s[0] + s1[0];
  }
}
