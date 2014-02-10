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
 * File: grampc_setparam.c
 * Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Date: February 2014
 * Version: v1.0
 *
 * Parameter setting file for GRAMPC.
 *
 */


#include "grampc_setparam.h"
#include "grampc_mess.h"

void grampc_setparam_real(typeGRAMPC *grampc, typeChar paramName[], typeRNum paramValue)
{
  int i;
  /* Set prediction horizon */
  if (!strncmp(paramName,"Thor",NAME_THOR)) {
    grampc->param->Thor = paramValue;
    if (grampc->param->Thor <= 0) {
      grampc_error_addstring(INVALID_PARAM_VALUE,paramName);
    }
    for (i = 0; i <= grampc->param->Nhor-1; i++) {
      grampc->rws->t[i] = grampc->param->Thor/(grampc->param->Nhor-1)*i;
    }
  }
  /* Set sampling time */
  else if (!strncmp(paramName,"dt",NAME_DT)) {
    grampc->param->dt = paramValue;
    if (grampc->param->dt <= 0) {
      grampc_error_addstring(INVALID_PARAM_VALUE,paramName);
    }
  }
  /* Set current time */
  else if (!strncmp(paramName,"tk",NAME_TK)) {
    grampc->param->tk = paramValue;
  }
  /* Undefined parameter */
  else {
    grampc_error_addstring(INVALID_PARAM_NAME,paramName);
  }
}

void grampc_setparam_int(typeGRAMPC *grampc, typeChar paramName[], typeInt paramValue)
{
  int i, j;
  /* Set discretization points */
  if (!strncmp(paramName,"Nhor",NAME_NHOR)) {
    grampc->param->Nhor = paramValue;
    if (grampc->param->Nhor <= 0) {
      grampc_error_addstring(INVALID_PARAM_VALUE,paramName);
    }
    /* Reallocation and -definition */
    /* t */
    /* myFree(grampc->rws->t); */
    /* grampc->rws->t = (typeRNum *)myCalloc(grampc->param->Nhor,sizeof(*(grampc->rws)->t)); */
    grampc->rws->t = (typeRNum *)myRealloc(grampc->rws->t,grampc->param->Nhor * sizeof(*(grampc->rws)->t));
    if (grampc->rws->t == NULL) {
      grampc_error(RWS_ELEMENT_ALLOC_FAILED);
    }
    /* x */
    /* myFree(grampc->rws->x); */
    /* grampc->rws->x = (typeRNum *)myCalloc(grampc->param->Nhor * grampc->param->Nx,sizeof(*(grampc->rws)->x)); */
    grampc->rws->x = (typeRNum *)myRealloc(grampc->rws->x,grampc->param->Nhor * grampc->param->Nx * sizeof(*(grampc->rws)->x));
    if (grampc->rws->x == NULL) {
      grampc_error(RWS_ELEMENT_ALLOC_FAILED);
    }
    /* adj */
    /* myFree(grampc->rws->adj); */
    /* grampc->rws->adj = (typeRNum *)myCalloc(grampc->param->Nhor * grampc->param->Nx,sizeof(*(grampc->rws)->adj)); */
    grampc->rws->adj = (typeRNum *)myRealloc(grampc->rws->adj,grampc->param->Nhor * grampc->param->Nx * sizeof(*(grampc->rws)->adj));
    if (grampc->rws->adj == NULL) {
      grampc_error(RWS_ELEMENT_ALLOC_FAILED);
    }
    /* u */
    /* myFree(grampc->rws->u); */
    /* grampc->rws->u = (typeRNum *)myCalloc(grampc->param->Nhor * grampc->param->Nu,sizeof(*(grampc->rws)->u)); */
    grampc->rws->u = (typeRNum *)myRealloc(grampc->rws->u,grampc->param->Nhor * grampc->param->Nu * sizeof(*(grampc->rws)->u));
    if (grampc->rws->u == NULL) {
      grampc_error(RWS_ELEMENT_ALLOC_FAILED);
    }
    /* dHdu */
    /* myFree(grampc->rws->dHdu); */
    /* grampc->rws->dHdu = (typeRNum *)myCalloc(grampc->param->Nhor * grampc->param->Nu,sizeof(*(grampc->rws)->dHdu)); */
    grampc->rws->dHdu = (typeRNum *)myRealloc(grampc->rws->dHdu,grampc->param->Nhor * grampc->param->Nu * sizeof(*(grampc->rws)->dHdu));
    if (grampc->rws->dHdu == NULL) {
      grampc_error(RWS_ELEMENT_ALLOC_FAILED);
    }
    /* uls */
    /* myFree(grampc->rws->uls); */
    /* grampc->rws->uls = (typeRNum *)myCalloc(grampc->param->Nhor * grampc->param->Nu,sizeof(*(grampc->rws)->uls)); */
    grampc->rws->uls = (typeRNum *)myRealloc(grampc->rws->uls,grampc->param->Nhor * grampc->param->Nu * sizeof(*(grampc->rws)->uls));
    if (grampc->rws->uls == NULL) {
      grampc_error(RWS_ELEMENT_ALLOC_FAILED);
    }
    /* uprev */
    /* myFree(grampc->rws->uprev); */
    /* grampc->rws->uprev = (typeRNum *)myCalloc(grampc->param->Nhor * grampc->param->Nu,sizeof(*(grampc->rws)->uprev)); */
    grampc->rws->uprev = (typeRNum *)myRealloc(grampc->rws->uprev,grampc->param->Nhor * grampc->param->Nu * sizeof(*(grampc->rws)->uprev));
    if (grampc->rws->uprev == NULL) {
      grampc_error(RWS_ELEMENT_ALLOC_FAILED);
    }
    /* dHduprev */
    /* myFree(grampc->rws->dHduprev); */
    /* grampc->rws->dHduprev = (typeRNum *)myCalloc(grampc->param->Nhor * grampc->param->Nu,sizeof(*(grampc->rws)->dHduprev)); */
    grampc->rws->dHduprev = (typeRNum *)myRealloc(grampc->rws->dHduprev,grampc->param->Nhor * grampc->param->Nu * sizeof(*(grampc->rws)->dHduprev));
    if (grampc->rws->dHduprev == NULL) {
      grampc_error(RWS_ELEMENT_ALLOC_FAILED);
    }
    /* J */
    /* myFree(grampc->rws->J); */
    /* grampc->rws->J = (typeRNum *)myCalloc(grampc->param->Nhor,sizeof(*(grampc->rws)->J)); */
    grampc->rws->J = (typeRNum *)myRealloc(grampc->rws->J,grampc->param->Nhor * sizeof(*(grampc->rws)->J));
    if (grampc->rws->J == NULL) {
      grampc_error(RWS_ELEMENT_ALLOC_FAILED);
    }
    /* time vector */
    for (i = 0; i <= grampc->param->Nhor-1; i++) {
      grampc->rws->t[i] = grampc->param->Thor/(grampc->param->Nhor-1)*i;
    }
    /* control vector */
    if (grampc->param->u0 != NULL) {
      for (i = 0; i <= grampc->param->Nhor-1; i++) {
        for (j = 0; j <= grampc->param->Nu-1; j++) {
          if (!strncmp(grampc->opt->ScaleProblem,"on",VALUE_ONOFF)) {
            grampc->rws->u[j+i*grampc->param->Nu] = (grampc->param->u0[j]-grampc->param->uOffset[j])/grampc->param->uScale[j];
          }
          else {
            grampc->rws->u[j+i*grampc->param->Nu] = grampc->param->u0[j];
          }
        }
      }
    }
  }
  /* NpCost */
  else if (!strncmp(paramName,"NpCost",NAME_NPCOST)) {
    grampc->param->NpCost = paramValue;
    if (grampc->param->NpCost < 0) {
      grampc_error_addstring(INVALID_PARAM_VALUE,paramName);
    }
  }
  /* NpSys */
  else if (!strncmp(paramName,"NpSys",NAME_NPSYS)) {
    grampc->param->NpSys = paramValue;
    if (grampc->param->NpSys < 0) {
      grampc_error_addstring(INVALID_PARAM_VALUE,paramName);
    }
  }
  /* Undefined parameter */
  else {
    grampc_error_addstring(INVALID_PARAM_NAME,paramName);
  }
}

void grampc_setparam_vector(typeGRAMPC *grampc, typeChar paramName[], typeRNum *paramValue)
{
  int i,j;
  /* xk */
  if (!strncmp(paramName,"xk",NAME_XK)) {
    if (grampc->param->xk == NULL) {
      grampc->param->xk = (typeRNum *)myCalloc(grampc->param->Nx,sizeof(*grampc->param->xk));
      if (grampc->param->xk == NULL) {
        grampc_error(PARAM_ALLOC_FAILED);
      }
    }
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      grampc->param->xk[i] = paramValue[i];
    }
  }
  /* u0 */
  else if (!strncmp(paramName,"u0",NAME_U0)) {
    if (grampc->param->u0 == NULL) {
      grampc->param->u0 = (typeRNum *)myCalloc(grampc->param->Nu,sizeof(*grampc->param->u0));
      if (grampc->param->u0 == NULL) {
        grampc_error(PARAM_ALLOC_FAILED);
      }
    }
    if (!strncmp(grampc->opt->ScaleProblem,"on",VALUE_ONOFF)) {
      for (j = 0; j <= grampc->param->Nu-1; j++) {
	grampc->param->u0[j] = paramValue[j];
	for (i = 0; i <= grampc->param->Nhor-1; i++) {
	  grampc->rws->u[j+i*grampc->param->Nu] = (grampc->param->u0[j]-grampc->param->uOffset[j])/grampc->param->uScale[j];
	}
      }
    }
    else {
      for (j = 0; j <= grampc->param->Nu-1; j++) {
	grampc->param->u0[j] = paramValue[j];
	for (i = 0; i <= grampc->param->Nhor-1; i++) {
	  grampc->rws->u[j+i*grampc->param->Nu] = grampc->param->u0[j];
	}
      }
    }
  }
  /* xdes */
  else if (!strncmp(paramName,"xdes",NAME_XDES)) {
    if (grampc->param->xdes == NULL) {
      grampc->param->xdes = (typeRNum *)myCalloc(grampc->param->Nx,sizeof(*grampc->param->xdes));
      if (grampc->param->xdes == NULL) {
        grampc_error(PARAM_ALLOC_FAILED);
      }
    }
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      grampc->param->xdes[i] = paramValue[i];
    }
  }
  /* udes */
  else if (!strncmp(paramName,"udes",NAME_UDES)) {
    if (grampc->param->udes == NULL) {
      grampc->param->udes = (typeRNum *)myCalloc(grampc->param->Nu,sizeof(*grampc->param->udes));
      if (grampc->param->udes == NULL) {
        grampc_error(PARAM_ALLOC_FAILED);
      }
    }
    for (i = 0; i <= grampc->param->Nu-1; i++) {
      grampc->param->udes[i] = paramValue[i];
    }
  }
  /* pCost */
  else if (!strncmp(paramName,"pCost",NAME_PCOST)) {
    if (grampc->param->NpCost == 0) {
      grampc_error(NPCOST_ZERO);
    }
    myFree(grampc->param->pCost);
    grampc->param->pCost = (typeRNum *)myCalloc(grampc->param->NpCost,sizeof(*grampc->param->pCost));
    if (grampc->param->pCost == NULL) {
      grampc_error(PARAM_ALLOC_FAILED);
    }
    for (i = 0; i <= grampc->param->NpCost-1; i++) {
      grampc->param->pCost[i] = paramValue[i];
    }
  }
  /* pSys */
  else if (!strncmp(paramName,"pSys",NAME_PSYS)) {
    if (grampc->param->NpSys == 0) {
      grampc_error(NPSYS_ZERO);
    }
    myFree(grampc->param->pSys);
    grampc->param->pSys = (typeRNum *)myCalloc(grampc->param->NpSys,sizeof(*grampc->param->pSys));
    if (grampc->param->pSys == NULL) {
      grampc_error(PARAM_ALLOC_FAILED);
    }
    for (i = 0; i <= grampc->param->NpSys-1; i++) {
      grampc->param->pSys[i] = paramValue[i];
    }
  }
  /* umax */
  else if (!strncmp(paramName,"umax",NAME_UMAX)) {
    for (i = 0; i <= grampc->param->Nu-1; i++) {
      grampc->param->umax[i] = paramValue[i];
    }
  }
  /* umin */
  else if (!strncmp(paramName,"umin",NAME_UMIN)) {
    for (i = 0; i <= grampc->param->Nu-1; i++) {
      grampc->param->umin[i] = paramValue[i];
    }
  }
  /* xScale */
  else if (!strncmp(paramName,"xScale",NAME_XSCALE)) {
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      if (grampc->param->xScale[i] == 0) {
        grampc_error(STATESCALE_VALUE_ZERO);
      }
      grampc->param->xScale[i] = paramValue[i];
    }
  }
  /* xOffset */
  else if (!strncmp(paramName,"xOffset",NAME_XOFFSET)) {
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      grampc->param->xOffset[i] = paramValue[i];
    }
  }
  /* uScale */
  else if (!strncmp(paramName,"uScale",NAME_USCALE)) { 
    if (!strncmp(grampc->opt->ScaleProblem,"on",VALUE_ONOFF) && grampc->param->u0 != NULL) {
      for (j = 0; j <= grampc->param->Nu-1; j++) {
	grampc->param->uScale[j] = paramValue[j];
	if (grampc->param->uScale[j] == 0) {
	  grampc_error(CONTROLSCALE_VALUE_ZERO);  
	}      
      	for (i = 0; i <= grampc->param->Nhor-1; i++) {
      	  grampc->rws->u[j+i*grampc->param->Nu] = (grampc->param->u0[j]-grampc->param->uOffset[j])/grampc->param->uScale[j];
      	}
      }
    }
    else {
      for (j = 0; j <= grampc->param->Nu-1; j++) {
	grampc->param->uScale[j] = paramValue[j];
	if (grampc->param->uScale[j] == 0) {
	  grampc_error(CONTROLSCALE_VALUE_ZERO);  
	}    
      }
    }
  }
  /* uOffset */
  else if (!strncmp(paramName,"uOffset",NAME_UOFFSET)) {
    if (!strncmp(grampc->opt->ScaleProblem,"on",VALUE_ONOFF) && grampc->param->u0 != NULL) {
      for (j = 0; j <= grampc->param->Nu-1; j++) {
	grampc->param->uOffset[j] = paramValue[j];
      	for (i = 0; i <= grampc->param->Nhor-1; i++) {
      	  grampc->rws->u[j+i*grampc->param->Nu] = (grampc->param->u0[j]-grampc->param->uOffset[j])/grampc->param->uScale[j];
      	}
      }
    }
    else {
      for (j = 0; j <= grampc->param->Nu-1; j++) {
	grampc->param->uOffset[j] = paramValue[j];
      }
    }
  }
  /* Undefined parameter */
  else {
    grampc_error_addstring(INVALID_PARAM_NAME,paramName);
  }
}

void grampc_printparam(typeGRAMPC *grampc)
{
  int i;
  myPrint("%s","-- MPC PARAMETERS --\n");
  myPrint("     Nx: %d\n",   grampc->param->Nx);
  myPrint("     Nu: %d\n",   grampc->param->Nu);

  if (grampc->param->xk == NULL) {
    myPrint("%s","     xk: []\n");
  }
  else {
    myPrint("%s","     xk: [");
    for (i = 0; i <= grampc->param->Nx-2; i++) {
      myPrint("%.2f,", grampc->param->xk[i]);
    }
    myPrint("%.2f]\n", grampc->param->xk[grampc->param->Nx-1]);
  }

  if (grampc->param->u0 == NULL) {
    myPrint("%s","     u0: []\n");
  }
  else {
    myPrint("%s","     u0: [");
    for (i = 0; i <= grampc->param->Nu-2; i++) {
      myPrint("%.2f,", grampc->param->u0[i]);
    }
    myPrint("%.2f]\n", grampc->param->u0[grampc->param->Nu-1]);
  }

  if (grampc->param->xdes == NULL) {
    myPrint("%s","   xdes: []\n");
  }
  else {
    myPrint("%s","   xdes: [");
    for (i = 0; i <= grampc->param->Nx-2; i++) {
      myPrint("%.2f,", grampc->param->xdes[i]);
    }
    myPrint("%.2f]\n", grampc->param->xdes[grampc->param->Nx-1]);
  }

  if (grampc->param->udes == NULL) {
    myPrint("%s","   udes: []\n");
  }
  else {
    myPrint("%s","   udes: [");
    for (i = 0; i <= grampc->param->Nu-2; i++) {
      myPrint("%.2f,", grampc->param->udes[i]);
    }
    myPrint("%.2f]\n", grampc->param->udes[grampc->param->Nu-1]);
  }

  myPrint("   Thor: %.2f\n", grampc->param->Thor);
  myPrint("     dt: %.4f\n", grampc->param->dt);
  myPrint("     tk: %.4f\n", grampc->param->tk);
  myPrint("   Nhor: %d\n",   grampc->param->Nhor);

  if (grampc->param->pCost == NULL) {
    myPrint("%s","  pCost: []\n");
  }
  else {
    myPrint("%s","  pCost: [");
    for (i = 0; i <= grampc->param->NpCost-2; i++) {
      myPrint("%.2f,", grampc->param->pCost[i]);
    }
    myPrint("%.2f]\n", grampc->param->pCost[grampc->param->NpCost-1]);
  }
  myPrint(" NpCost: %d\n", grampc->param->NpCost);

  if (grampc->param->pSys == NULL) {
    myPrint("%s","   pSys: []\n");
  }
  else {
    myPrint("%s","   pSys: [");
    for (i = 0; i <= grampc->param->NpSys-2; i++) {
      myPrint("%.2f,", grampc->param->pSys[i]);
    }
    myPrint("%.2f]\n", grampc->param->pSys[grampc->param->NpSys-1]);
  }
  myPrint("  NpSys: %d\n", grampc->param->NpSys);

  myPrint("%s","   umax: [");
  for (i = 0; i <= grampc->param->Nu-2; i++) {
    if (grampc->param->umax[i] == INF) {
      myPrint("%s","inf,");
    }
    else {
      myPrint("%.2f,", grampc->param->umax[i]);
    }
  }
  if (grampc->param->umax[grampc->param->Nu-1] == INF) {
    myPrint("%s","inf]\n");
  }
  else {
    myPrint("%.2f]\n", grampc->param->umax[grampc->param->Nu-1]);
  }

  myPrint("%s","   umin: [");
  for (i = 0; i <= grampc->param->Nu-2; i++) {
    if (grampc->param->umin[i] == -INF) {
      myPrint("%s","-inf,");
    }
    else {
      myPrint("%.2f,", grampc->param->umin[i]);
    }
  }
  if (grampc->param->umin[grampc->param->Nu-1] == -INF) {
    myPrint("%s","-inf]\n");
  }
  else {
    myPrint("%.2f]\n", grampc->param->umin[grampc->param->Nu-1]);
  }

  myPrint("%s"," xScale: [");
  for (i = 0; i <= grampc->param->Nx-2; i++) {
    myPrint("%.2f,", grampc->param->xScale[i]);
  }
  myPrint("%.2f]\n", grampc->param->xScale[grampc->param->Nx-1]);

  myPrint("%s","xOffset: [");
  for (i = 0; i <= grampc->param->Nx-2; i++) {
    myPrint("%.2f,", grampc->param->xOffset[i]);
  }
  myPrint("%.2f]\n", grampc->param->xOffset[grampc->param->Nx-1]);

  myPrint("%s"," uScale: [");
  for (i = 0; i <= grampc->param->Nu-2; i++) {
    myPrint("%.2f,", grampc->param->uScale[i]);
  }
  myPrint("%.2f]\n", grampc->param->uScale[grampc->param->Nu-1]);

  myPrint("%s","uOffset: [");
  for (i = 0; i <= grampc->param->Nu-2; i++) {
    myPrint("%.2f,", grampc->param->uOffset[i]);
  }
  myPrint("%.2f]\n", grampc->param->uOffset[grampc->param->Nu-1]);
}
