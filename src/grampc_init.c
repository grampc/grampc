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
 * File: grampc_init.c
 * Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Date: February 2014
 * Version: v1.0
 *
 * Initialization file for GRAMPC.
 *
 */


#include "grampc_init.h"
#include "grampc_mess.h"
#include "probfct.h"

void grampc_init(typeGRAMPC **grampc)
{
  typeInt i;
  size_t Nrws;

  /* MEMORY ALLOCATION OF GRAMPC STRUCTURE ***********************************/
  *grampc = (typeGRAMPC *)myCalloc(1,sizeof(**grampc));
  if (*grampc == NULL) {
    grampc_error(GRAMPC_ALLOC_FAILED);
  }
  /* END OF GRAMPC ALLOCATION ************************************************/

  /* MEMORY ALLOCATION OF PARAM STRUCTURE ************************************/
  (*grampc)->param = (typeGRAMPCparam *)myCalloc(1,sizeof(*(*grampc)->param));
  if ((*grampc)->param == NULL) {
    grampc_error(PARAM_ALLOC_FAILED);
  }
  /* END OF PARAM ALLOCATION *************************************************/

  /* MEMORY ALLOCATION OF MPC SOLUTION STRUCTURE *****************************/
  (*grampc)->sol = (typeGRAMPCsol *)myCalloc(1,sizeof(*(*grampc)->sol));
  if ((*grampc)->sol == NULL) {
    grampc_error(SOL_ALLOC_FAILED);
  }
  /* END OF MPC SOLUTION ALLOCATION ******************************************/

  /* MEMORY ALLOCATION OF WORKSPACE ARRAY ************************************/
  (*grampc)->rws = (typeGRAMPCrws *)myCalloc(1,sizeof(*(*grampc)->rws));
  if ((*grampc)->rws == NULL) {
    grampc_error(RWS_ALLOC_FAILED);
  }
  /* END WORKSPACE ALLOCATION ************************************************/

  /* MEMORY ALLOCATION OF OPTIONS ARRAY **************************************/
  (*grampc)->opt = (typeGRAMPCopt *)myCalloc(1,sizeof(*(*grampc)->opt));
  if ((*grampc)->opt == NULL) {
    grampc_error(OPT_ALLOC_FAILED);
  }
  /* END OPTIONS ALLOCATION **************************************************/

  /* Initialization of MPC options (default values) **************************/
  strncpy((*grampc)->opt->ShiftControl,"on",VALUE_ONOFF);
  strncpy((*grampc)->opt->ScaleProblem,"off",VALUE_ONOFF);
  strncpy((*grampc)->opt->CostIntegrator,"trapezodial",VALUE_COSTINTMETHOD);
  strncpy((*grampc)->opt->Integrator,"heun",VALUE_INTEGRATOR);
  strncpy((*grampc)->opt->LineSearchType,"adaptive",VALUE_LSTYPE);
  strncpy((*grampc)->opt->JacobianX,"sysjacxadj",VALUE_JACOBIANX);
  strncpy((*grampc)->opt->JacobianU,"sysjacuadj",VALUE_JACOBIANU);
  strncpy((*grampc)->opt->IntegralCost,"on",VALUE_ONOFF);
  strncpy((*grampc)->opt->FinalCost,"on",VALUE_ONOFF);

  (*grampc)->opt->MaxIter                  = 2;
  (*grampc)->opt->IntegratorRelTol         = 1e-6;
  (*grampc)->opt->IntegratorAbsTol         = 1e-8;
  (*grampc)->opt->LineSearchMax            = 0.75;
  (*grampc)->opt->LineSearchMin            = 1e-5;
  (*grampc)->opt->LineSearchInit           = 5.5e-4;
  (*grampc)->opt->LineSearchIntervalFactor = 0.85;
  (*grampc)->opt->LineSearchAdaptFactor    = 3.0/2.0;
  (*grampc)->opt->LineSearchIntervalTol    = 1e-1;
  /* End of initialization of MPC options (default values) *******************/

  /* Initialization of MPC structure (default values ) ***********************/
  /* Nx, Nu */
  sysdim(&(*grampc)->param->Nx,&(*grampc)->param->Nu);
  /* check number of states and control */
  if ((*grampc)->param->Nx <= 0) {
    grampc_error(INVALID_NX);
  }
  if ((*grampc)->param->Nu <= 0) {
    grampc_error(INVALID_NU);
  }
  /* Thor, dt, tk, Nhor, NpSys, NpCost (default values) */
  (*grampc)->param->Thor   = -1.0;
  (*grampc)->param->Nhor   = 30;
  (*grampc)->param->dt     = -1.0; 
  (*grampc)->param->tk     = 0.0;
  (*grampc)->param->NpSys  = 0;
  (*grampc)->param->NpCost = 0;
  /* xk and u0*/
  (*grampc)->param->xk = NULL;
  (*grampc)->param->u0 = NULL;
  /* xdes and udes */
  (*grampc)->param->xdes = NULL;
  (*grampc)->param->udes = NULL;
  /* Parameters for system dynamics and cost functional */
  (*grampc)->param->pSys  = NULL;
  (*grampc)->param->pCost = NULL;
  /* umax and umin */
  (*grampc)->param->umax = (typeRNum *)myCalloc((*grampc)->param->Nu,sizeof(*(*grampc)->param->umax));
  if ((*grampc)->param->umax == NULL) {
    grampc_error(PARAM_ALLOC_FAILED);
  }
  (*grampc)->param->umin = (typeRNum *)myCalloc((*grampc)->param->Nu,sizeof(*(*grampc)->param->umin));
  if ((*grampc)->param->umin == NULL) {
    grampc_error(PARAM_ALLOC_FAILED);
  }
  for (i = 0; i <= (*grampc)->param->Nu-1; i++) {
    (*grampc)->param->umax[i] = +INF;
    (*grampc)->param->umin[i] = -INF;
  }
  /* xScale and xOffset */
  (*grampc)->param->xScale = (typeRNum *)myCalloc((*grampc)->param->Nx,sizeof(*(*grampc)->param->xScale));
  if ((*grampc)->param->xScale == NULL) {
    grampc_error(PARAM_ALLOC_FAILED);
  }
  (*grampc)->param->xOffset = (typeRNum *)myCalloc((*grampc)->param->Nx,sizeof(*(*grampc)->param->xOffset));
  if ((*grampc)->param->xOffset == NULL) {
    grampc_error(PARAM_ALLOC_FAILED);
  }
  for (i = 0; i <= (*grampc)->param->Nx-1; i++) {
    (*grampc)->param->xScale[i]  = 1;
    (*grampc)->param->xOffset[i] = 0;
  }
  /* uScale and uOffset */
  (*grampc)->param->uScale = (typeRNum *)myCalloc((*grampc)->param->Nu,sizeof(*(*grampc)->param->uScale));
  if ((*grampc)->param->uScale == NULL) {
    grampc_error(PARAM_ALLOC_FAILED);
  }
  (*grampc)->param->uOffset = (typeRNum *)myCalloc((*grampc)->param->Nu,sizeof(*(*grampc)->param->uOffset));
  if ((*grampc)->param->uOffset == NULL) {
    grampc_error(PARAM_ALLOC_FAILED);
  }
  for (i = 0; i <= (*grampc)->param->Nu-1; i++) {
    (*grampc)->param->uScale[i]  = 1;
    (*grampc)->param->uOffset[i] = 0;
  }
  /* End of initialization of MPC structure (default values) *****************/

  /* Initialization of MPC solution structure ********************************/
  (*grampc)->sol->xnext = (typeRNum *)myCalloc((*grampc)->param->Nx,sizeof(*(*grampc)->sol->xnext));
  if ((*grampc)->sol->xnext == NULL) {
    grampc_error(SOL_ALLOC_FAILED);
  }
  (*grampc)->sol->unext = (typeRNum *)myCalloc((*grampc)->param->Nu,sizeof(*(*grampc)->sol->unext));
  if ((*grampc)->sol->unext == NULL) {
    grampc_error(SOL_ALLOC_FAILED);
  }
  (*grampc)->sol->J = (typeRNum *)myCalloc(1,sizeof(*(*grampc)->sol->J));
  if ((*grampc)->sol->J == NULL) {
    grampc_error(SOL_ALLOC_FAILED);
  }
  /* End Initialization of MPC solution structure ****************************/

  /* Initialization of workspace array rws ***********************************/
  /* t */
  (*grampc)->rws->t = (typeRNum *)myCalloc((*grampc)->param->Nhor,sizeof(*(*grampc)->rws->t));
  if ((*grampc)->rws->t == NULL) {
    grampc_error(RWS_ELEMENT_ALLOC_FAILED);
  }
  /* x */
  (*grampc)->rws->x = (typeRNum *)myCalloc((*grampc)->param->Nhor * (*grampc)->param->Nx,sizeof(*(*grampc)->rws->x));
  if ((*grampc)->rws->x == NULL) {
    grampc_error(RWS_ELEMENT_ALLOC_FAILED);
  }
  /* adj */
  (*grampc)->rws->adj = (typeRNum *)myCalloc((*grampc)->param->Nhor * (*grampc)->param->Nx,sizeof(*(*grampc)->rws->adj));
  if ((*grampc)->rws->adj == NULL) {
    grampc_error(RWS_ELEMENT_ALLOC_FAILED);
  }
  /* u */
  (*grampc)->rws->u = (typeRNum *)myCalloc((*grampc)->param->Nhor * (*grampc)->param->Nu,sizeof(*(*grampc)->rws->u));
  if ((*grampc)->rws->u == NULL) {
    grampc_error(RWS_ELEMENT_ALLOC_FAILED);
  }
  /* dHdu */
  (*grampc)->rws->dHdu = (typeRNum *)myCalloc((*grampc)->param->Nhor * (*grampc)->param->Nu,sizeof(*(*grampc)->rws->dHdu));
  if ((*grampc)->rws->dHdu == NULL) {
    grampc_error(RWS_ELEMENT_ALLOC_FAILED);
  }
  /* lsAdapt */
  (*grampc)->rws->lsAdapt = (typeRNum *)myCalloc(2*(NLS + 1)*(1+(*grampc)->opt->MaxIter),sizeof(*(*grampc)->rws->lsAdapt));
  if ((*grampc)->rws->lsAdapt == NULL) {
    grampc_error(RWS_ELEMENT_ALLOC_FAILED);
  }
  /* uls */
  (*grampc)->rws->uls = (typeRNum *)myCalloc((*grampc)->param->Nhor * (*grampc)->param->Nu,sizeof(*(*grampc)->rws->uls));
  if ((*grampc)->rws->uls == NULL) {
    grampc_error(RWS_ELEMENT_ALLOC_FAILED);
  }
  /* lsExplicit */
  (*grampc)->rws->lsExplicit = (typeRNum *)myCalloc(NLS,sizeof(*(*grampc)->rws->lsExplicit));
  if ((*grampc)->rws->lsExplicit == NULL) {
    grampc_error(RWS_ELEMENT_ALLOC_FAILED);
  }
  /* uprev */
  (*grampc)->rws->uprev = (typeRNum *)myCalloc((*grampc)->param->Nhor * (*grampc)->param->Nu,sizeof(*(*grampc)->rws->uprev));
  if ((*grampc)->rws->uprev == NULL) {
    grampc_error(RWS_ELEMENT_ALLOC_FAILED);
  }
  /* dHduprev */
  (*grampc)->rws->dHduprev = (typeRNum *)myCalloc((*grampc)->param->Nhor * (*grampc)->param->Nu,sizeof(*(*grampc)->rws->dHduprev));
  if ((*grampc)->rws->dHduprev == NULL) {
    grampc_error(RWS_ELEMENT_ALLOC_FAILED);
  }
  /* J */
  (*grampc)->rws->J = (typeRNum *)myCalloc((*grampc)->param->Nhor,sizeof(*(*grampc)->rws->J));
  if ((*grampc)->rws->J == NULL) {
    grampc_error(RWS_ELEMENT_ALLOC_FAILED);
  }
  /* rwsScale */
  Nrws = 2*(*grampc)->param->Nx + 2*(*grampc)->param->Nu;
  (*grampc)->rws->rwsScale = (typeRNum *)myCalloc(Nrws,sizeof(*(*grampc)->rws->rwsScale));
  if ((*grampc)->rws->rwsScale == NULL) {
    grampc_error(RWS_ELEMENT_ALLOC_FAILED);
  }
  /* rwsGradient */
  Nrws = (*grampc)->param->Nu*((*grampc)->param->Nx + 2);
  (*grampc)->rws->rwsGradient = (typeRNum *)myCalloc(Nrws,sizeof(*(*grampc)->rws->rwsGradient));
  if ((*grampc)->rws->rwsGradient == NULL) {
    grampc_error(RWS_ELEMENT_ALLOC_FAILED);
  }
  /* rwsCostIntegration */
  Nrws = 3 + (*grampc)->param->Nx + (*grampc)->param->Nu;
  (*grampc)->rws->rwsCostIntegration = (typeRNum *)myCalloc(Nrws,sizeof(*(*grampc)->rws->rwsCostIntegration));
  if ((*grampc)->rws->rwsCostIntegration == NULL) {
    grampc_error(RWS_ELEMENT_ALLOC_FAILED);
  }
  /* rwsAdjIntegration*/
  Nrws = (*grampc)->param->Nx*((*grampc)->param->Nx + 1);
  (*grampc)->rws->rwsAdjIntegration = (typeRNum *)myCalloc(Nrws,sizeof(*(*grampc)->rws->rwsAdjIntegration));
  if ((*grampc)->rws->rwsAdjIntegration == NULL) {
    grampc_error(RWS_ELEMENT_ALLOC_FAILED);
  }
  /* rwsIntegration */
  Nrws = 17*(*grampc)->param->Nx + (*grampc)->param->Nu;
  (*grampc)->rws->rwsIntegration = (typeRNum *)myCalloc(Nrws,sizeof(*(*grampc)->rws->rwsIntegration));
  if ((*grampc)->rws->rwsIntegration == NULL) {
    grampc_error(RWS_ELEMENT_ALLOC_FAILED);
  }
  /* End initialization of workspace array (*grampc)->rws ****************************/

  /* Initialization **********************************************************/
  /* time vector */
  for (i = 0; i <= (*grampc)->param->Nhor-1; i++) {
    (*grampc)->rws->t[i] = (*grampc)->param->Thor/((*grampc)->param->Nhor-1)*i;
  }
  /* linesearch */
  /* adaptive line search */
  for (i = 0; i <= (*grampc)->opt->MaxIter; i++) {
    (*grampc)->rws->lsAdapt[0+i*2*(NLS+1)] = (*grampc)->opt->LineSearchInit*(1-(*grampc)->opt->LineSearchIntervalFactor);
    (*grampc)->rws->lsAdapt[1+i*2*(NLS+1)] = (*grampc)->opt->LineSearchInit;
    (*grampc)->rws->lsAdapt[2+i*2*(NLS+1)] = (*grampc)->opt->LineSearchInit*(1+(*grampc)->opt->LineSearchIntervalFactor);
    (*grampc)->rws->lsAdapt[3+i*2*(NLS+1)] = (*grampc)->opt->LineSearchInit;
  }
  /* explicit line search */
  (*grampc)->rws->lsExplicit[2] = (*grampc)->opt->LineSearchInit;
  /* End of initialization ***************************************************/
}

void grampc_free(typeGRAMPC **grampc)
{
  myFree((*grampc)->rws->t);
  myFree((*grampc)->rws->x);
  myFree((*grampc)->rws->adj);
  myFree((*grampc)->rws->u);
  myFree((*grampc)->rws->dHdu);
  myFree((*grampc)->rws->lsAdapt);
  myFree((*grampc)->rws->J);
  myFree((*grampc)->rws->uls);
  myFree((*grampc)->rws->lsExplicit);
  myFree((*grampc)->rws->uprev);
  myFree((*grampc)->rws->dHduprev);
  myFree((*grampc)->rws->rwsScale);
  myFree((*grampc)->rws->rwsGradient);
  myFree((*grampc)->rws->rwsCostIntegration);
  myFree((*grampc)->rws->rwsAdjIntegration);
  myFree((*grampc)->rws->rwsIntegration);
  myFree((*grampc)->rws);

  myFree((*grampc)->param->xk);
  myFree((*grampc)->param->u0);
  myFree((*grampc)->param->xdes);
  myFree((*grampc)->param->udes);
  myFree((*grampc)->param->umax);
  myFree((*grampc)->param->umin);
  myFree((*grampc)->param->xScale);
  myFree((*grampc)->param->xOffset);
  myFree((*grampc)->param->uScale);
  myFree((*grampc)->param->uOffset);
  myFree((*grampc)->param->pCost);
  myFree((*grampc)->param->pSys);
  myFree((*grampc)->param);

  myFree((*grampc)->sol->xnext);
  myFree((*grampc)->sol->unext);
  myFree((*grampc)->sol->J);
  myFree((*grampc)->sol);

  myFree((*grampc)->opt);

  myFree(*grampc);
}
