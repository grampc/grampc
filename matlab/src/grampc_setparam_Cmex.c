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
 * File: grampc_setparam_Cmex.c
 * Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Date: February 2014
 * Version: v1.0
 * 
 * Mex function interface to use the function grampc_setparam for setting
 * MPC parameters of GRAMPC.
 *
 */
#define nrwsOut                  9

#include "mex.h"
#include "matrix.h"

#include "grampc_init.h"
#include "grampc_setparam.h"
#include "grampc_mess.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  typeInt i, idx;
  typeRNum *p      = NULL;
  mxArray *mxopt   = NULL;
  mxArray *mxparam = NULL;
  mxArray *mxrws   = NULL;

  typeChar *paramname;

  mxArray *paramOut;
  mxArray *rwsOut[nrwsOut];

	const char *rwsOutName[nrwsOut] = { "t", "x", "adj", "u", "dHdu", "uls", "uprev", "dHduprev", "J" };

  typeGRAMPC grampc;
  typeGRAMPCopt opt;
  typeGRAMPCparam param;
  typeGRAMPCrws rws;

  /* check proper number of input arguments */
  if (nrhs < 3) {
    mexErrMsgTxt("Not enough input arguments.");
  }
  if (nrhs > 3) {
    mexErrMsgTxt("Too many input arguments.");
  }
  /* first argument: param structure */
  if (!mxIsStruct(prhs[0])) {
    mexErrMsgTxt("First argument (grampc structure) must be a struct.");
  }
  /* Second argument: option name */
  if (!mxIsChar(prhs[1])) {
    mexErrMsgTxt("Second argument (option name) must be a string.");
  }
  /* Third argument: option value */
  if (!mxIsNumeric(prhs[2])) {
    mexErrMsgTxt("Third argument (option value) must be a number.");
  }
  /* check proper number of output arguments */
  if (nlhs != 0) {
    mexErrMsgTxt("No output arguments.");
  }

  mxopt   = mxGetField(prhs[0],0,"opt");
  mxparam = mxGetField(prhs[0],0,"param");
  mxrws   = mxGetField(prhs[0],0,"rws");

  paramname = mxArrayToString(prhs[1]);

  param.Nx = (typeInt)mxGetScalar(mxGetField(mxparam,0,"Nx"));
  param.Nu = (typeInt)mxGetScalar(mxGetField(mxparam,0,"Nu"));
  
  if (mxGetField(mxparam,0,"xk") == NULL) {
    param.xk = NULL;
  }
  else {
    param.xk = mxGetPr(mxGetField(mxparam,0,"xk"));
  }
  
  if (mxGetField(mxparam,0,"u0") == NULL) {
    param.u0 = NULL;
  }
  else {
    param.u0 = mxGetPr(mxGetField(mxparam,0,"u0"));
  }
  
  if (mxGetField(mxparam,0,"xdes") == NULL) {
    param.xdes = NULL;
  }
  else {
    param.xdes = mxGetPr(mxGetField(mxparam,0,"xdes"));
  }
  
  if (mxGetField(mxparam,0,"udes") == NULL) {
    param.udes = NULL;
  }
  else {
    param.udes = mxGetPr(mxGetField(mxparam,0,"udes"));
  }
  param.Thor    = (typeRNum)mxGetScalar(mxGetField(mxparam,0,"Thor"));
  param.dt      = (typeRNum)mxGetScalar(mxGetField(mxparam,0,"dt"));
  param.tk      = (typeRNum)mxGetScalar(mxGetField(mxparam,0,"tk"));
  param.Nhor    = (typeInt)mxGetScalar(mxGetField(mxparam,0,"Nhor"));
  param.umax    = mxGetPr(mxGetField(mxparam,0,"umax"));
  param.umin    = mxGetPr(mxGetField(mxparam,0,"umin"));
  param.xScale  = mxGetPr(mxGetField(mxparam,0,"xScale"));
  param.xOffset = mxGetPr(mxGetField(mxparam,0,"xOffset"));
  param.uScale  = mxGetPr(mxGetField(mxparam,0,"uScale"));
  param.uOffset = mxGetPr(mxGetField(mxparam,0,"uOffset"));
  param.NpCost	 = (typeInt)mxGetScalar(mxGetField(mxparam,0,"NpCost"));
  param.NpSys   = (typeInt)mxGetScalar(mxGetField(mxparam,0,"NpSys"));
  
  if (mxGetField(mxparam,0,"pCost") == NULL ||
      !strncmp(paramname,"pCost",NAME_PCOST)) {
    param.pCost = NULL;
  }
  else {
    param.pCost = mxGetPr(mxGetField(mxparam,0,"pCost"));
  }
  if (mxGetField(mxparam,0,"pSys") == NULL ||
      !strncmp(paramname,"pSys",NAME_PSYS)) {
    param.pSys = NULL;
  }
  else {
    param.pSys = mxGetPr(mxGetField(mxparam,0,"pSys"));
  }

  /* option structure */
  opt.MaxIter 		       = (typeInt)mxGetScalar(mxGetField(mxopt,0,"MaxIter"));
  opt.IntegratorRelTol 	       = mxGetScalar(mxGetField(mxopt,0,"IntegratorRelTol"));
  opt.IntegratorAbsTol 	       = mxGetScalar(mxGetField(mxopt,0,"IntegratorAbsTol"));
  opt.LineSearchMax 	       = mxGetScalar(mxGetField(mxopt,0,"LineSearchMax"));
  opt.LineSearchMin 	       = mxGetScalar(mxGetField(mxopt,0,"LineSearchMin"));
  opt.LineSearchInit 	       = mxGetScalar(mxGetField(mxopt,0,"LineSearchInit"));
  opt.LineSearchIntervalFactor = mxGetScalar(mxGetField(mxopt,0,"LineSearchIntervalFactor"));
  opt.LineSearchAdaptFactor    = mxGetScalar(mxGetField(mxopt,0,"LineSearchAdaptFactor"));
  opt.LineSearchIntervalTol    = mxGetScalar(mxGetField(mxopt,0,"LineSearchIntervalTol"));

  mxGetString(mxGetField(mxopt, 0, "ShiftControl"), opt.ShiftControl, VALUE_ONOFF);
  mxGetString(mxGetField(mxopt, 0, "ScaleProblem"), opt.ScaleProblem, VALUE_ONOFF);
  mxGetString(mxGetField(mxopt, 0, "CostIntegrator"), opt.CostIntegrator, VALUE_COSTINTMETHOD);
  mxGetString(mxGetField(mxopt, 0, "Integrator"), opt.Integrator, VALUE_INTEGRATOR);
  mxGetString(mxGetField(mxopt, 0, "LineSearchType"), opt.LineSearchType, VALUE_LSTYPE);
  mxGetString(mxGetField(mxopt, 0, "JacobianX"), opt.JacobianX, VALUE_JACOBIANX);
  mxGetString(mxGetField(mxopt, 0, "JacobianU"), opt.JacobianU, VALUE_JACOBIANU);
  mxGetString(mxGetField(mxopt, 0, "IntegralCost"), opt.IntegralCost, VALUE_ONOFF);
  mxGetString(mxGetField(mxopt, 0, "FinalCost"), opt.FinalCost, VALUE_ONOFF);

  /* rws structure */
  if (!strncmp(paramname,"Nhor",NAME_NHOR)) {
    rws.t        = NULL;
    rws.x        = NULL;
    rws.adj      = NULL;
    rws.u        = NULL;
    rws.dHdu     = NULL;
    rws.uls      = NULL;
    rws.uprev    = NULL;
    rws.dHduprev = NULL;
    rws.J        = NULL;
  }
  if (!strncmp(paramname,"Thor",NAME_THOR)) {
    rws.t = mxGetPr(mxGetField(mxrws,0,"t"));
  }
  if (!strncmp(paramname,"u0",NAME_U0) ||
      !strncmp(paramname,"uScale",NAME_USCALE) ||
      !strncmp(paramname,"uOffset",NAME_UOFFSET)) {
    rws.u = mxGetPr(mxGetField(mxrws,0,"u"));
  }

  /* check array length of pointer supposed to be changed */
  if ( (!strncmp(paramname,"xk",NAME_XK)) ||
       (!strncmp(paramname,"xdes",NAME_XDES)) ||
       (!strncmp(paramname,"xScale",NAME_XSCALE)) ||
       (!strncmp(paramname,"xOffset",NAME_XOFFSET)) ) {
    if (mxGetNumberOfElements(prhs[2]) != param.Nx) {
      mexErrMsgTxt("No. of elements of input vector does not correspond to Nx.");
    }
  }
  if ( (!strncmp(paramname,"u0",NAME_U0)) ||
       (!strncmp(paramname,"udes",NAME_UDES)) ||
       (!strncmp(paramname,"uScale",NAME_USCALE)) ||
       (!strncmp(paramname,"uOffset",NAME_UOFFSET)) ||
       (!strncmp(paramname,"umax",NAME_UMAX)) ||
       (!strncmp(paramname,"umin",NAME_UMIN)) ) {
    if (mxGetNumberOfElements(prhs[2]) != param.Nu) {
      mexErrMsgTxt("No. of elements of input vector does not correspond to Nu.");
    }
  }
  if (!strncmp(paramname,"pSys",NAME_PSYS)) {
    if (mxGetNumberOfElements(prhs[2]) != param.NpSys) {
      mexErrMsgTxt("No. of elements of pSys does not correspond to NpSys.");
    }
  }
  if (!strncmp(paramname,"pCost",NAME_PCOST)) {
    if (mxGetNumberOfElements(prhs[2]) != param.NpCost) {
      mexErrMsgTxt("No. of elements of pCost does not correspond to NpCost.");
    }
  }
  
  grampc.param = &param;
  grampc.opt   = &opt;
  grampc.rws   = &rws;

  /* Set parameters */
  if ( (!strncmp(paramname,"Nhor",NAME_NHOR)) || (!strncmp(paramname,"NpCost",NAME_NPCOST)) || (!strncmp(paramname,"NpSys",NAME_NPSYS)) ) {
    grampc_setparam_int(&grampc,paramname,(typeInt)mxGetScalar(prhs[2]));
  }
  else if ( (!strncmp(paramname,"Thor",NAME_THOR)) || (!strncmp(paramname,"dt",NAME_DT)) || (!strncmp(paramname,"tk",NAME_TK)) ){
    grampc_setparam_real(&grampc,paramname,mxGetScalar(prhs[2]));
  }
  else {
    grampc_setparam_vector(&grampc,paramname,mxGetPr(prhs[2]));
  }

  /* param STRUCTURE **********************************************************/
  /* xk */
  if (!strncmp(paramname,"xk",NAME_XK)) {
    if (param.xk == NULL) {
      paramOut = NULL;
    }
    else {
      paramOut = mxCreateDoubleMatrix(1,param.Nx,mxREAL);
      for (i = 0; i <= param.Nx-1; i++) {
        *(mxGetPr(paramOut)+i) = param.xk[i];
      }
    }
  }
  /* u0 */
  else if (!strncmp(paramname,"u0",NAME_U0)) {
    if (param.u0 == NULL) {
      paramOut = NULL;
    }
    else {
      paramOut = mxCreateDoubleMatrix(1,param.Nu,mxREAL);
      for (i = 0; i <= param.Nu-1; i++) {
        *(mxGetPr(paramOut)+i) = param.u0[i];
      }
    }
  }
  /* xdes */
  else if (!strncmp(paramname,"xdes",NAME_XDES)) {
    if (param.xdes == NULL) {
      paramOut = NULL;
    }
    else {
      paramOut = mxCreateDoubleMatrix(1,param.Nx,mxREAL);
      for (i = 0; i <= param.Nx-1; i++) {
        *(mxGetPr(paramOut)+i) = param.xdes[i];
      }
    }
  }
  /* udes */
  else if (!strncmp(paramname,"udes",NAME_UDES)) {
    if (param.udes == NULL) {
      paramOut = NULL;
    }
    else {
      paramOut = mxCreateDoubleMatrix(1,param.Nu,mxREAL);
      for (i = 0; i <= param.Nu-1; i++) {
        *(mxGetPr(paramOut)+i) = param.udes[i];
      }
    }
  }
  /* Thor */
  else if (!strncmp(paramname,"Thor",NAME_THOR)) {
    paramOut = mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(paramOut)) = param.Thor;
  }
  /* dt */
  else if (!strncmp(paramname,"dt",NAME_DT)) {
    paramOut = mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(paramOut)) = param.dt;
  }
  /* tk */
  else if (!strncmp(paramname,"tk",NAME_TK)) {
    paramOut = mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(paramOut)) = param.tk;
  }
  /* Nhor */
  else if (!strncmp(paramname,"Nhor",NAME_NHOR)) {
    paramOut = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    *((typeInt *)mxGetData(paramOut)) = param.Nhor;
  }
  /* pCost */
  else if (!strncmp(paramname,"pCost",NAME_PCOST)) {
    if (param.pCost == NULL) {
      paramOut = NULL;
    }
    else {
      paramOut = mxCreateDoubleMatrix(1,param.NpCost,mxREAL);
      for (i = 0; i <= param.NpCost-1; i++) {
        *(mxGetPr(paramOut)+i) = param.pCost[i];
      }
	  free(param.pCost);
    }
  }
  /* pSys */
  else if (!strncmp(paramname,"pSys",NAME_PSYS)) {
    if (param.pSys == NULL) {
      paramOut = NULL;
    }
    else {
      paramOut = mxCreateDoubleMatrix(1,param.NpSys,mxREAL);
      for (i = 0; i <= param.NpSys-1; i++) {
        *(mxGetPr(paramOut)+i) = param.pSys[i];
      }
	  free(param.pSys);
    }
  }
  /* umax */
  else if (!strncmp(paramname,"umax",NAME_UMAX)) {
    paramOut = mxCreateDoubleMatrix(1,param.Nu,mxREAL);
    for (i = 0; i <= param.Nu-1; i++) {
      *(mxGetPr(paramOut)+i) = param.umax[i];
    }
  }
  /* umin */
  else if (!strncmp(paramname,"umin",NAME_UMIN)) {
    paramOut = mxCreateDoubleMatrix(1,param.Nu,mxREAL);
    for (i = 0; i <= param.Nu-1; i++) {
      *(mxGetPr(paramOut)+i) = param.umin[i];
    }
  }
  /* xScale */
  else if (!strncmp(paramname,"xScale",NAME_XSCALE)) {
    paramOut = mxCreateDoubleMatrix(1,param.Nx,mxREAL);
    for (i = 0; i <= param.Nx-1; i++) {
      *(mxGetPr(paramOut)+i) = param.xScale[i];
    }
  }
  /* xOffset */
  else if (!strncmp(paramname,"xOffset",NAME_XOFFSET)) {
    paramOut = mxCreateDoubleMatrix(1,param.Nx,mxREAL);
    for (i = 0; i <= param.Nx-1; i++) {
      *(mxGetPr(paramOut)+i) = param.xOffset[i];
    }
  }
  /* uScale */
  else if (!strncmp(paramname,"uScale",NAME_USCALE)) {
    paramOut = mxCreateDoubleMatrix(1,param.Nu,mxREAL);
    for (i = 0; i <= param.Nu-1; i++) {
      *(mxGetPr(paramOut)+i) = param.uScale[i];
    }
  }
  /* uOffset */
  else if (!strncmp(paramname,"uOffset",NAME_UOFFSET)) {
    paramOut = mxCreateDoubleMatrix(1,param.Nu,mxREAL);
    for (i = 0; i <= param.Nu-1; i++) {
      *(mxGetPr(paramOut)+i) = param.uOffset[i];
    }
  }
  /* NpCost */
  else if (!strncmp(paramname,"NpCost",NAME_NPCOST)) {
    paramOut = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    *((typeInt *)mxGetData(paramOut)) = param.NpCost;
  }
  /* NpSys */
  else if (!strncmp(paramname,"NpSys",NAME_NPSYS)) {
    paramOut = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    *((typeInt *)mxGetData(paramOut)) = param.NpSys;
  }
  /* Set field */
  mxDestroyArray(mxGetField(mxparam, 0, paramname));
  mxSetField(mxparam,0,paramname,paramOut);
  /**************************************************************************/

  /* RWS STRUCTURE **********************************************************/
  if (!strncmp(paramname,"Nhor",NAME_NHOR)) {
    /* t */
    idx = 0;
    rwsOut[idx] = mxCreateDoubleMatrix(1,param.Nhor,mxREAL);
    for (i = 0; i <= param.Nhor-1; i++) {
      *(mxGetPr(rwsOut[idx])+i) = rws.t[i];
    }
	free(rws.t);
    /* x */
    idx += 1;
    rwsOut[idx] = mxCreateDoubleMatrix(param.Nx,param.Nhor,mxREAL);
    for (i = 0; i <= param.Nhor*param.Nx-1; i++) {
      *(mxGetPr(rwsOut[idx])+i) = rws.x[i];
    }
	free(rws.x);
    /* adj */
    idx += 1;
    rwsOut[idx] = mxCreateDoubleMatrix(param.Nx,param.Nhor,mxREAL);
    for (i = 0; i <= param.Nhor*param.Nx-1; i++) {
      *(mxGetPr(rwsOut[idx])+i) = rws.adj[i];
    }
	free(rws.adj);
    /* u */
    idx += 1;
    rwsOut[idx] = mxCreateDoubleMatrix(param.Nu,param.Nhor,mxREAL);
    for (i = 0; i <= param.Nhor*param.Nu-1; i++) {
      *(mxGetPr(rwsOut[idx])+i) = rws.u[i];
    }
	free(rws.u);
    /* dHdu */
    idx += 1;
    rwsOut[idx] = mxCreateDoubleMatrix(param.Nu,param.Nhor,mxREAL);
    for (i = 0; i <= param.Nhor*param.Nu-1; i++) {
      *(mxGetPr(rwsOut[idx])+i) = rws.dHdu[i];
    }
	free(rws.dHdu);
    /* uls */
    idx += 1;
    rwsOut[idx] = mxCreateDoubleMatrix(param.Nu,param.Nhor,mxREAL);
    for (i = 0; i <= param.Nhor*param.Nu-1; i++) {
      *(mxGetPr(rwsOut[idx])+i) = rws.uls[i];
    }
	free(rws.uls);
    /* uprev */
    idx += 1;
    rwsOut[idx] = mxCreateDoubleMatrix(param.Nu,param.Nhor,mxREAL);
    for (i = 0; i <= param.Nhor*param.Nu-1; i++) {
      *(mxGetPr(rwsOut[idx])+i) = rws.uprev[i];
    }
	free(rws.uprev);
    /* dHduprev */
    idx += 1;
    rwsOut[idx] = mxCreateDoubleMatrix(param.Nu,param.Nhor,mxREAL);
    for (i = 0; i <= param.Nhor*param.Nu-1; i++) {
      *(mxGetPr(rwsOut[idx])+i) = rws.dHduprev[i];
    }
	free(rws.dHduprev);
    /* J */
    idx += 1;
    rwsOut[idx] = mxCreateDoubleMatrix(1,param.Nhor,mxREAL);
    for (i = 0; i <= param.Nhor-1; i++) {
      *(mxGetPr(rwsOut[idx])+i) = rws.J[i];
    }
	free(rws.J);
    /* Set field */
    for (i = 0; i <= nrwsOut-1; i++) {
	  mxDestroyArray(mxGetField(mxrws, 0, rwsOutName[i]));
      mxSetField(mxrws,0,rwsOutName[i],rwsOut[i]);
	  
    }
  }
  mxFree(paramname);
  /**************************************************************************/

}
