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
 * File: grampc_init_Cmex.c
 * Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Date: February 2014
 * Version: v1.0
 * 
 * Mex function interface to use the function grampc_init 
 * for initialization of GRAMPC.
 *
 */

#define nGRAMPCparam      	20
#define nGRAMPCrws      	16
#define nGRAMPCopt		18
#define nGRAMPC                 3

/* include headers and functions */
#include "mex.h"
#include "matrix.h"
#include "grampc.h"


/* global variable */
const mxArray *mxMPC = NULL;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  typeInt i;
  typeInt idx;

  typeGRAMPC *grampc;

  typeInt Nrws;

  mxArray *p = NULL;
  mxArray *paramOut[nGRAMPCparam];
  mxArray *rwsOut[nGRAMPCrws];
  mxArray *optOut[nGRAMPCopt];

  const char *grampcOutName[nGRAMPC] = {"param","rws","opt"};
  mxArray *grampcOutStructure[nGRAMPC];

  const char *paramOutName[nGRAMPCparam]       = {"Nx","Nu","xk","u0","xdes","udes","Thor","dt","tk","Nhor","pCost","pSys","umax","umin","xScale","xOffset","uScale","uOffset","NpCost","NpSys"};

  const char *rwsOutName[nGRAMPCrws]       = {"t","x","adj","u","dHdu","lsAdapt","uls","lsExplicit","uprev","dHduprev","J","rwsScale","rwsGradient","rwsCostIntegration","rwsAdjIntegration","rwsIntegration"};

  const char *optOutName[nGRAMPCopt] = {"MaxIter","ShiftControl","ScaleProblem","CostIntegrator","Integrator","IntegratorRelTol","IntegratorAbsTol","LineSearchType","LineSearchMax","LineSearchMin","LineSearchInit","LineSearchIntervalFactor","LineSearchAdaptFactor","LineSearchIntervalTol","JacobianX","JacobianU","IntegralCost","FinalCost"};

  /* check proper number of input arguments */
  if (nrhs != 0) {
    mexErrMsgTxt("Function has no input arguments.");
  }

  if (nlhs > 3) {
    mexErrMsgTxt("Too many output arguments.");
  }

  /* Memory allocation and initialization */
  grampc_init(&grampc);

  /* generation of outputs */
  plhs[0] = mxCreateStructMatrix(1,1,nGRAMPC,grampcOutName);

  grampcOutStructure[0] = mxCreateStructMatrix(1,1,nGRAMPCparam,paramOutName);
  grampcOutStructure[1] = mxCreateStructMatrix(1,1,nGRAMPCrws,rwsOutName);
  grampcOutStructure[2] = mxCreateStructMatrix(1,1,nGRAMPCopt,optOutName);

  /* param STRUCTURE **********************************************************/
  /* Nx */
  idx = 0;
  paramOut[idx] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
  *((typeInt *)mxGetData(paramOut[idx])) = grampc->param->Nx;
  /* Nu */
  idx += 1;
  paramOut[idx] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
  *((typeInt *)mxGetData(paramOut[idx])) = grampc->param->Nu;
  /* xk */
  idx += 1;
  paramOut[idx] = NULL;
  /* paramOut[idx] = mxCreateDoubleMatrix(1,1,mxREAL);  */
  /* *(mxGetPr(paramOut[idx])) = mxGetNaN();  */
  /* u0 */
  idx += 1;
  paramOut[idx] = NULL;
  /* paramOut[idx] = mxCreateDoubleMatrix(1,1,mxREAL); */
  /* *(mxGetPr(paramOut[idx])) = mxGetNaN(); */
  /* xdes */
  idx += 1;
  paramOut[idx] = NULL;
  /* paramOut[idx] = mxCreateDoubleMatrix(1,1,mxREAL); */
  /* *(mxGetPr(paramOut[idx])) = mxGetNaN(); */
  /* udes */
  idx += 1;
  paramOut[idx] = NULL;
  /* paramOut[idx] = mxCreateDoubleMatrix(1,1,mxREAL); */
  /* *(mxGetPr(paramOut[idx])) = mxGetNaN(); */
  /* Thor */
  idx += 1;
  paramOut[idx] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(paramOut[idx])) = grampc->param->Thor;
  /* dt */
  idx += 1;
  paramOut[idx] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(paramOut[idx])) = grampc->param->dt;
  /* tk */
  idx += 1;
  paramOut[idx] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(paramOut[idx])) = grampc->param->tk;
  /* Nhor */
  idx += 1;
  paramOut[idx] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
  *((typeInt *)mxGetData(paramOut[idx])) = grampc->param->Nhor;
  /* pCost */
  idx += 1;
  paramOut[idx] = NULL;
  /* paramOut[idx] = mxCreateDoubleMatrix(1,1,mxREAL); */
  /* *(mxGetPr(paramOut[idx])) = mxGetNaN(); */
  /* pSys */
  idx += 1;
  paramOut[idx] = NULL;
  /* paramOut[idx] = mxCreateDoubleMatrix(1,1,mxREAL); */
  /* *(mxGetPr(paramOut[idx])) = mxGetNaN(); */
  /* umax */
  idx += 1;
  paramOut[idx] = mxCreateDoubleMatrix(1,grampc->param->Nu,mxREAL);
  for (i = 0; i <= grampc->param->Nu-1; i++) {
    *(mxGetPr(paramOut[idx])+i) = grampc->param->umax[i];
  }
  /* umin */
  idx += 1;
  paramOut[idx] = mxCreateDoubleMatrix(1,grampc->param->Nu,mxREAL);
  for (i = 0; i <= grampc->param->Nu-1; i++) {
    *(mxGetPr(paramOut[idx])+i) = grampc->param->umin[i];
  }
  /* xScale and xOffset*/
  idx += 1;
  paramOut[idx]   = mxCreateDoubleMatrix(1,grampc->param->Nx,mxREAL);
  paramOut[idx+1] = mxCreateDoubleMatrix(1,grampc->param->Nx,mxREAL);
  for (i = 0; i <= grampc->param->Nx-1; i++) {
    *(mxGetPr(paramOut[idx])+i)   = grampc->param->xScale[i];
    *(mxGetPr(paramOut[idx+1])+i) = grampc->param->xOffset[i];
  }
  /* uScale and uOffset */
  idx += 2;
  paramOut[idx]   = mxCreateDoubleMatrix(1,grampc->param->Nu,mxREAL);
  paramOut[idx+1] = mxCreateDoubleMatrix(1,grampc->param->Nu,mxREAL);
  for (i = 0; i <= grampc->param->Nu-1; i++) {
    *(mxGetPr(paramOut[idx])+i)   = grampc->param->uScale[i];
    *(mxGetPr(paramOut[idx+1])+i) = grampc->param->uOffset[i];
  }
  /* NpCost */
  idx += 2;
  paramOut[idx] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
  *((typeInt *)mxGetData(paramOut[idx])) = grampc->param->NpCost;
  /* NpSys */
  idx += 1;
  paramOut[idx] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
  *((typeInt *)mxGetData(paramOut[idx])) = grampc->param->NpSys;
  /* Set param Structure */
  for (i = 0; i <= nGRAMPCparam-1; i++) {
    mxSetFieldByNumber(grampcOutStructure[0],0,i,paramOut[i]);
  }
  /**************************************************************************/

  /* RWS STRUCTURE **********************************************************/
  /* t */
  idx = 0;
  rwsOut[idx] = mxCreateDoubleMatrix(1,grampc->param->Nhor,mxREAL);
  for (i = 0; i <= grampc->param->Nhor-1; i++) {
    *(mxGetPr(rwsOut[idx])+i) = grampc->rws->t[i];
  }
  /* x */
  idx += 1;
  rwsOut[idx] = mxCreateDoubleMatrix(grampc->param->Nx,grampc->param->Nhor,mxREAL);
  for (i = 0; i <= grampc->param->Nhor*grampc->param->Nx-1; i++) {
    *(mxGetPr(rwsOut[idx])+i) = grampc->rws->x[i];
  }
  /* adj */
  idx += 1;
  rwsOut[idx] = mxCreateDoubleMatrix(grampc->param->Nx,grampc->param->Nhor,mxREAL);
  for (i = 0; i <= grampc->param->Nhor*grampc->param->Nx-1; i++) {
    *(mxGetPr(rwsOut[idx])+i) = grampc->rws->adj[i];
  }
  /* u */
  idx += 1;
  rwsOut[idx] = mxCreateDoubleMatrix(grampc->param->Nu,grampc->param->Nhor,mxREAL);
  for (i = 0; i <= grampc->param->Nhor*grampc->param->Nu-1; i++) {
    *(mxGetPr(rwsOut[idx])+i) = grampc->rws->u[i];
  }
  /* dHdu */
  idx += 1;
  rwsOut[idx] = mxCreateDoubleMatrix(grampc->param->Nu,grampc->param->Nhor,mxREAL);
  for (i = 0; i <= grampc->param->Nhor*grampc->param->Nu-1; i++) {
    *(mxGetPr(rwsOut[idx])+i) = grampc->rws->dHdu[i];
  } 
  /* lsAdapt */
  idx += 1;
  rwsOut[idx] = mxCreateDoubleMatrix(1,2*(NLS+1)*(1+grampc->opt->MaxIter),mxREAL);
  for (i = 0; i <= 2*(NLS+1)*(1+grampc->opt->MaxIter)-1; i++) {
    *(mxGetPr(rwsOut[idx])+i) = grampc->rws->lsAdapt[i];
  }
  /* uls */
  idx += 1;
  rwsOut[idx] = mxCreateDoubleMatrix(grampc->param->Nu,grampc->param->Nhor,mxREAL);
  for (i = 0; i <= grampc->param->Nhor*grampc->param->Nu-1; i++) {
    *(mxGetPr(rwsOut[idx])+i) = grampc->rws->uls[i];
  }
  /* lsExplicit */
  idx += 1;
  rwsOut[idx] = mxCreateDoubleMatrix(1,NLS,mxREAL);
  for (i = 0; i <= NLS-1; i++) {
    *(mxGetPr(rwsOut[idx])+i) = grampc->rws->lsExplicit[i];
  }
  /* uprev */
  idx += 1;
  rwsOut[idx] = mxCreateDoubleMatrix(grampc->param->Nu,grampc->param->Nhor,mxREAL);
  for (i = 0; i <= grampc->param->Nhor*grampc->param->Nu-1; i++) {
    *(mxGetPr(rwsOut[idx])+i) = grampc->rws->uprev[i];
  }
  /* dHduprev */
  idx += 1;
  rwsOut[idx] = mxCreateDoubleMatrix(grampc->param->Nu,grampc->param->Nhor,mxREAL);
  for (i = 0; i <= grampc->param->Nhor*grampc->param->Nu-1; i++) {
    *(mxGetPr(rwsOut[idx])+i) = grampc->rws->dHduprev[i];
  }
  /* J */
  idx += 1;
  rwsOut[idx] = mxCreateDoubleMatrix(1,grampc->param->Nhor,mxREAL);
  for (i = 0; i <= grampc->param->Nhor-1; i++) {
    *(mxGetPr(rwsOut[idx])+i) = grampc->rws->J[i];
  }
  /* rwsScaling */
  Nrws = 2*grampc->param->Nx + 2*grampc->param->Nu;
  idx += 1;
  rwsOut[idx] = mxCreateDoubleMatrix(1,Nrws,mxREAL);
  for (i = 0; i <= Nrws-1; i++) {
    *(mxGetPr(rwsOut[idx])+i) = grampc->rws->rwsScale[i];
  }
  /* rwsGradient */
  Nrws = grampc->param->Nu*(grampc->param->Nx + 2);
  idx += 1;
  rwsOut[idx] = mxCreateDoubleMatrix(1,Nrws,mxREAL);
  for (i = 0; i <= Nrws-1; i++) {
    *(mxGetPr(rwsOut[idx])+i) = grampc->rws->rwsGradient[i];
  }
  /* rwsCostIntegration */
  Nrws = 4 + grampc->param->Nx + grampc->param->Nu;
  idx += 1;
  rwsOut[idx] = mxCreateDoubleMatrix(1,Nrws,mxREAL);
  for (i = 0; i <= Nrws-1; i++) {
    *(mxGetPr(rwsOut[idx])+i) = grampc->rws->rwsCostIntegration[i];
  }
  /* rwsAdjIntegration */
  Nrws = grampc->param->Nx*(grampc->param->Nx + 1);
  idx += 1;
  rwsOut[idx] = mxCreateDoubleMatrix(1,Nrws,mxREAL);
  for (i = 0; i <= Nrws-1; i++) {
    *(mxGetPr(rwsOut[idx])+i) = grampc->rws->rwsAdjIntegration[i];
  }
  /* rwsIntegration */
  Nrws = 17*grampc->param->Nx + grampc->param->Nu;
  idx += 1;
  rwsOut[idx] = mxCreateDoubleMatrix(1,Nrws,mxREAL);
  for (i = 0; i <= Nrws-1; i++) {
    *(mxGetPr(rwsOut[idx])+i) = grampc->rws->rwsIntegration[i];
  }
  /* Set Structure */
  for (i = 0; i <= nGRAMPCrws-1; i++) {
    mxSetFieldByNumber(grampcOutStructure[1],0,i,rwsOut[i]);
  }
  /**************************************************************************/

  /* OPTIONS STRUCTURE ******************************************************/
  /* MaxIter */
  idx = 0;
  optOut[idx] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
  *((typeInt *)mxGetData(optOut[idx])) = grampc->opt->MaxIter;
  /* ShiftControl */
  idx += 1;
  optOut[idx] = mxCreateString(grampc->opt->ShiftControl);
  /* ScaleProblem */
  idx += 1;
  optOut[idx] = mxCreateString(grampc->opt->ScaleProblem);
  /* CostIntegrator */
  idx += 1;
  optOut[idx] = mxCreateString(grampc->opt->CostIntegrator);
  /* Integrator */
  idx += 1;
  optOut[idx] = mxCreateString(grampc->opt->Integrator);
  /* IntegratorRelTol */
  idx += 1;
  optOut[idx] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(optOut[idx])) = grampc->opt->IntegratorRelTol;
  /* IntegratorAbsTol */
  idx += 1;
  optOut[idx] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(optOut[idx])) = grampc->opt->IntegratorAbsTol;
  /* LineSearchType */
  idx += 1;
  optOut[idx] = mxCreateString(grampc->opt->LineSearchType);
  /* LineSearchMax */
  idx += 1;
  optOut[idx] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(optOut[idx])) = grampc->opt->LineSearchMax;
  /* LineSearchMin */
  idx += 1;
  optOut[idx] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(optOut[idx])) = grampc->opt->LineSearchMin;
  /* LineSearchInit */
  idx += 1;
  optOut[idx] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(optOut[idx])) = grampc->opt->LineSearchInit;
  /* LineSearchIntervalFactor */
  idx += 1;
  optOut[idx] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(optOut[idx])) = grampc->opt->LineSearchIntervalFactor;
  /* LineSearchAdaptFactor */
  idx += 1;
  optOut[idx] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(optOut[idx])) = grampc->opt->LineSearchAdaptFactor;
  /* LineSearchIntervalTol */
  idx += 1;
  optOut[idx] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(optOut[idx])) = grampc->opt->LineSearchIntervalTol;
  /* JacobianX */
  idx += 1;
  optOut[idx] = mxCreateString(grampc->opt->JacobianX);
  /* JacobianU */
  idx += 1;
  optOut[idx] = mxCreateString(grampc->opt->JacobianU);
  /* IntegralCost */
  idx += 1;
  optOut[idx] = mxCreateString(grampc->opt->IntegralCost);
  /* fCostFct */
  idx += 1;
  optOut[idx] = mxCreateString(grampc->opt->FinalCost);
  /* Set Structure */
  for (i = 0; i <= nGRAMPCopt-1; i++) {
    mxSetFieldByNumber(grampcOutStructure[2],0,i,optOut[i]);
  }
  /**************************************************************************/

  /* GRAMPC STRUCTURE *******************************************************/
  for (i = 0; i <= nGRAMPC-1; i++) {
    mxSetFieldByNumber(plhs[0],0,i,grampcOutStructure[i]);
  }
  /**************************************************************************/

  /* free allocated memory */
  grampc_free(&grampc);
}
