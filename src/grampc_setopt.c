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
 * File: grampc_setopt.c
 * Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Date: February 2014
 * Version: v1.0
 *
 * Option setting file for GRAMPC.
 *
 */


#include "grampc_setopt.h"
#include "grampc_mess.h"

void grampc_setopt_real(typeGRAMPC *grampc, typeChar optName[], typeRNum optValue)
{
  typeInt i;
  /* Integrator relative tolerance */
  if (!strncmp(optName,"IntegratorRelTol",NAME_INTRELTOL)) {
    grampc->opt->IntegratorRelTol = optValue;
    if (grampc->opt->IntegratorRelTol <= 0.0) {
      grampc_error_addstring(INVALID_OPTION_VALUE,optName);
    }
  }
  /* Integrator absolute tolerance */
  else if (!strncmp(optName,"IntegratorAbsTol",NAME_INTABSTOL)) {
    grampc->opt->IntegratorAbsTol = optValue;
    if (grampc->opt->IntegratorAbsTol <= 0.0) {
      grampc_error_addstring(INVALID_OPTION_VALUE,optName);
    }
  }
  /* Max value linesearch */
  else if (!strncmp(optName,"LineSearchMax",NAME_LSMAX)) {
    grampc->opt->LineSearchMax = optValue;
    if (grampc->opt->LineSearchMax < 0.0) {
      grampc_error_addstring(INVALID_OPTION_VALUE,optName);
    }
  }
  /* Min value linesearch */
  else if (!strncmp(optName,"LineSearchMin",NAME_LSMIN)) {
    grampc->opt->LineSearchMin = optValue;
    if (grampc->opt->LineSearchMin < 0.0) {
      grampc_error_addstring(INVALID_OPTION_VALUE,optName);
    }
  }
  /* Init value linesearch */
  else if (!strncmp(optName,"LineSearchInit",NAME_LSINIT)) {
    grampc->opt->LineSearchInit = optValue;
    if (grampc->opt->LineSearchInit < 0.0) {
      grampc_error_addstring(INVALID_OPTION_VALUE,optName);
    }
    for (i = 0; i <= grampc->opt->MaxIter; i++) {
      grampc->rws->lsAdapt[0+i*2*(NLS+1)] = grampc->opt->LineSearchInit*(1-grampc->opt->LineSearchIntervalFactor);
      grampc->rws->lsAdapt[1+i*2*(NLS+1)] = grampc->opt->LineSearchInit;
      grampc->rws->lsAdapt[2+i*2*(NLS+1)] = grampc->opt->LineSearchInit*(1+grampc->opt->LineSearchIntervalFactor);
      grampc->rws->lsAdapt[3+i*2*(NLS+1)] = grampc->opt->LineSearchInit;
    }
    grampc->rws->lsExplicit[2] = grampc->opt->LineSearchInit;
  }
  /* LineSearchIntervalFactor */
  else if (!strncmp(optName,"LineSearchIntervalFactor",NAME_LSINTFACTOR)) {
    grampc->opt->LineSearchIntervalFactor = optValue;
    if (grampc->opt->LineSearchIntervalFactor <= 0.0 || grampc->opt->LineSearchIntervalFactor >= 1.0) {
      grampc_error_addstring(INVALID_OPTION_VALUE,optName);
    }
    for (i = 0; i <= grampc->opt->MaxIter; i++) {
      grampc->rws->lsAdapt[0+i*2*(NLS+1)] = grampc->opt->LineSearchInit*(1-grampc->opt->LineSearchIntervalFactor);
      grampc->rws->lsAdapt[1+i*2*(NLS+1)] = grampc->opt->LineSearchInit;
      grampc->rws->lsAdapt[2+i*2*(NLS+1)] = grampc->opt->LineSearchInit*(1+grampc->opt->LineSearchIntervalFactor);
      grampc->rws->lsAdapt[3+i*2*(NLS+1)] = grampc->opt->LineSearchInit;
    }
  }
  /* LineSearchAdaptFactor */
  else if (!strncmp(optName,"LineSearchAdaptFactor",NAME_LSADAPTFACTOR)) {
    grampc->opt->LineSearchAdaptFactor = optValue;
    if (grampc->opt->LineSearchAdaptFactor <= 1.0) {
      grampc_error_addstring(INVALID_OPTION_VALUE,optName);
    }
  }
  /* LineSearchIntervalTol */
  else if (!strncmp(optName,"LineSearchIntervalTol",NAME_LSINTTOL)) {
    grampc->opt->LineSearchIntervalTol = optValue;
    if (grampc->opt->LineSearchIntervalTol <= 0.0 || grampc->opt->LineSearchIntervalTol >= 1.0) {
      grampc_error_addstring(INVALID_OPTION_VALUE,optName);
    }
  }
  /* Undefined optName */
  else {
    grampc_error_addstring(INVALID_OPTION_NAME,optName);
  }
}

void grampc_setopt_int(typeGRAMPC *grampc, typeChar optName[], typeInt optValue)
{
  typeInt i;
  /* MaxIter */
  if (!strncmp(optName,"MaxIter",NAME_MAXITER)) {
    grampc->opt->MaxIter = optValue;
    if (grampc->opt->MaxIter <= 0) {
      grampc_error_addstring(INVALID_OPTION_VALUE,optName);
    }
    /* Reallocaton of lsAdapt */
    /* myFree(grampc->rws->lsAdapt); */
    /* grampc->rws->lsAdapt = (typeRNum *)myCalloc(2*(NLS + 1)*(1+grampc->opt->MaxIter),sizeof(*(grampc->rws)->lsAdapt)); */
    grampc->rws->lsAdapt = (typeRNum *)myRealloc(grampc->rws->lsAdapt,2*(NLS + 1)*(1+grampc->opt->MaxIter)*sizeof(*(grampc->rws)->lsAdapt));
    if (grampc->rws->lsAdapt == NULL) {
        grampc_error(RWS_ELEMENT_ALLOC_FAILED);
    }
    for (i = 0; i <= grampc->opt->MaxIter; i++) {
      grampc->rws->lsAdapt[0+i*2*(NLS+1)] = grampc->opt->LineSearchInit*(1-grampc->opt->LineSearchIntervalFactor);
      grampc->rws->lsAdapt[1+i*2*(NLS+1)] = grampc->opt->LineSearchInit;
      grampc->rws->lsAdapt[2+i*2*(NLS+1)] = grampc->opt->LineSearchInit*(1+grampc->opt->LineSearchIntervalFactor);
      grampc->rws->lsAdapt[3+i*2*(NLS+1)] = grampc->opt->LineSearchInit;
    }
  }
  /* Undefined optName */
  else {
    grampc_error_addstring(INVALID_OPTION_NAME,optName);
  }
}

void grampc_setopt_string(typeGRAMPC *grampc, typeChar optName[], typeChar optValue[])
{
  /* ControlShift */
  if (!strncmp(optName,"ShiftControl",NAME_SHIFTCONTROL)) {
    if (!strncmp(optValue,"on",VALUE_ONOFF)) {
      strncpy(grampc->opt->ShiftControl,optValue,VALUE_ONOFF);
    }
    else if (!strncmp(optValue,"off",VALUE_ONOFF)) {
      strncpy(grampc->opt->ShiftControl,optValue,VALUE_ONOFF);
    }
    else {
      grampc_error_addstring(INVALID_OPTION_VALUE,optName);
    }
  }
  /* ProblemScale */
  else if (!strncmp(optName,"ScaleProblem",NAME_SCALEPROBLEM)) {
    if (!strncmp(optValue,"on",VALUE_ONOFF)) {
      strncpy(grampc->opt->ScaleProblem,optValue,VALUE_ONOFF);
    }
    else if (!strncmp(optValue,"off",VALUE_ONOFF)) {
      strncpy(grampc->opt->ScaleProblem,optValue,VALUE_ONOFF);
    }
    else {
      grampc_error_addstring(INVALID_OPTION_VALUE,optName);
    }
  }
  /* CostIntegrator */
  else if (!strncmp(optName,"CostIntegrator",NAME_COSTINTMETHOD)) {
    if (!strncmp(optValue,"trapezodial",VALUE_COSTINTMETHOD)) {
      strncpy(grampc->opt->CostIntegrator,optValue,VALUE_COSTINTMETHOD);
    }
    else if (!strncmp(optValue,"simpson",VALUE_COSTINTMETHOD)) {
      strncpy(grampc->opt->CostIntegrator,optValue,VALUE_COSTINTMETHOD);
    }
    else {
      grampc_error_addstring(INVALID_OPTION_VALUE,optName);
    }
  }
  /* Integrator type */
  else if (!strncmp(optName,"Integrator",NAME_INTEGRATOR)) {
    if (!strncmp(optValue,"euler",VALUE_INTEGRATOR)) {
      strncpy(grampc->opt->Integrator,optValue,VALUE_INTEGRATOR);
    }
    else if (!strncmp(optValue,"modeuler",VALUE_INTEGRATOR)) {
      strncpy(grampc->opt->Integrator,optValue,VALUE_INTEGRATOR);
    }
    else if (!strncmp(optValue,"heun",VALUE_INTEGRATOR)) {
      strncpy(grampc->opt->Integrator,optValue,VALUE_INTEGRATOR);
    }
    else if (!strncmp(optValue,"ruku45",VALUE_INTEGRATOR)) {
      strncpy(grampc->opt->Integrator,optValue,VALUE_INTEGRATOR);
    }
    else {
      grampc_error_addstring(INVALID_OPTION_VALUE,optName);
    }
  }
  /* LineSearchType */
  else if (!strncmp(optName,"LineSearchType",NAME_LSTYPE)) {
    if (!strncmp(optValue,"adaptive",VALUE_LSTYPE)) {
      strncpy(grampc->opt->LineSearchType,optValue,VALUE_LSTYPE);
    }
    else if (!strncmp(optValue,"explicit1",VALUE_LSTYPE)) {
      strncpy(grampc->opt->LineSearchType,optValue,VALUE_LSTYPE);
    }
    else if (!strncmp(optValue,"explicit2",VALUE_LSTYPE)) {
      strncpy(grampc->opt->LineSearchType,optValue,VALUE_LSTYPE);
    }
    else {
      grampc_error_addstring(INVALID_OPTION_VALUE,optName);
    }
  }
  /* activate/deactivate JacobianX */
  else if (!strncmp(optName,"JacobianX",NAME_JACOBIANX)) {
    if (!strncmp(optValue,"sysjacx",VALUE_JACOBIANX)) {
      strncpy(grampc->opt->JacobianX,optValue,VALUE_JACOBIANX);
    }
    else if (!strncmp(optValue,"sysjacxadj",VALUE_JACOBIANX)) {
      strncpy(grampc->opt->JacobianX,optValue,VALUE_JACOBIANX);
    }
    else {
      grampc_error_addstring(INVALID_OPTION_VALUE,optName);
    }
  }
  /* activate/deactivate JacobianU */
  else if (!strncmp(optName,"JacobianU",NAME_JACOBIANU)) {
    if (!strncmp(optValue,"sysjacu",VALUE_JACOBIANU)) {
      strncpy(grampc->opt->JacobianU,optValue,VALUE_JACOBIANU);
    }
    else if (!strncmp(optValue,"sysjacuadj",VALUE_JACOBIANU)) {
      strncpy(grampc->opt->JacobianU,optValue,VALUE_JACOBIANU);
    }
    else {
      grampc_error_addstring(INVALID_OPTION_VALUE,optName);
    }
  }
  /* activate/deactivate IntegralCost */
  else if (!strncmp(optName,"IntegralCost",NAME_INTEGRALCOST)) {
    if (!strncmp(optValue,"on",VALUE_ONOFF)) {
      strncpy(grampc->opt->IntegralCost,optValue,VALUE_ONOFF);
    }
    else if (!strncmp(optValue,"off",VALUE_ONOFF)) {
      strncpy(grampc->opt->IntegralCost,optValue,VALUE_ONOFF);
    }
    else {
      grampc_error_addstring(INVALID_OPTION_VALUE,optName);
    }
  }
  /* activate/deactivate FinalCost */
  else if (!strncmp(optName,"FinalCost",NAME_FINALCOST)) {
    if (!strncmp(optValue,"on",VALUE_ONOFF)) {
      strncpy(grampc->opt->FinalCost,optValue,VALUE_ONOFF);
    }
    else if (!strncmp(optValue,"off",VALUE_ONOFF)) {
      strncpy(grampc->opt->FinalCost,optValue,VALUE_ONOFF);
    }
    else {
      grampc_error_addstring(INVALID_OPTION_VALUE,optName);
    }
  }
  /* Undefined optName */
  else {
    grampc_error_addstring(INVALID_OPTION_NAME,optName);
  }
}

void grampc_printopt(typeGRAMPC *grampc)
{
  myPrint("%s","-- MPC OPTIONS --\n"                                              );
  myPrint("                 MaxIter: %d\n",  grampc->opt->MaxIter                 );
  myPrint("            ShiftControl: %s\n",  grampc->opt->ShiftControl            );
  myPrint("            ScaleProblem: %s\n",  grampc->opt->ScaleProblem            );
  myPrint("          CostIntegrator: %s\n",  grampc->opt->CostIntegrator          );
  myPrint("              Integrator: %s\n",  grampc->opt->Integrator              );
  myPrint("        IntegratorRelTol: %.2e\n",grampc->opt->IntegratorRelTol        );
  myPrint("        IntegratorAbsTol: %.2e\n",grampc->opt->IntegratorAbsTol        );
  myPrint("          LineSearchType: %s\n",  grampc->opt->LineSearchType          );
  myPrint("           LineSearchMax: %f\n",  grampc->opt->LineSearchMax           );
  myPrint("           LineSearchMin: %f\n",  grampc->opt->LineSearchMin           );
  myPrint("          LineSearchInit: %f\n",  grampc->opt->LineSearchInit          );
  myPrint("LineSearchIntervalFactor: %f\n",  grampc->opt->LineSearchIntervalFactor);
  myPrint("   LineSearchAdaptFactor: %f\n",  grampc->opt->LineSearchAdaptFactor   );
  myPrint("   LineSearchIntervalTol: %f\n",  grampc->opt->LineSearchIntervalTol   );
  myPrint("               JacobianX: %s\n",  grampc->opt->JacobianX               );
  myPrint("               JacobianU: %s\n",  grampc->opt->JacobianU               );
  myPrint("            IntegralCost: %s\n",  grampc->opt->IntegralCost            );
  myPrint("               FinalCost: %s\n",  grampc->opt->FinalCost               );
}
