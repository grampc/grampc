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
 * File: grampc_setopt_Cmex.c
 * Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Date: February 2014
 * Version: v1.0
 * 
 * Mex function interface to use the function grampc_setopt for setting
 * MPC options of GRAMPC.
 *
 */

#define nrwsOut                  2

#include "mex.h"
#include "matrix.h"

#include "grampc_init.h"
#include "grampc_setopt.h"
#include "grampc_mess.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  typeInt i;
  mxArray *mxopt = NULL;
  mxArray *mxrws = NULL;

  mxArray *optOut;
  mxArray *rwsOut[nrwsOut];

  typeGRAMPC grampc;
  typeGRAMPCopt opt;
  typeGRAMPCrws rws;

  typeChar *optname;
  typeChar *optvalue;

  /* check proper number of input arguments */
  if (nrhs < 3) {
    mexErrMsgTxt("Not enough input arguments.");
  }
  if (nrhs > 3) {
    mexErrMsgTxt("Too many input arguments.");
  }
  /* first argument: grampc structure  */
  if (!mxIsStruct(prhs[0])) {
    mexErrMsgTxt("First argument (grampc structure) must be a struct.");
  }
  /* second argument: MPC structure */
  if (!mxIsChar(prhs[1])) {
    mexErrMsgTxt("Second argument (option name) must be a string.");
  }
  /* third argument: option structure */
  if (!mxIsNumeric(prhs[2]) && !mxIsChar(prhs[2])) {
    mexErrMsgTxt("Third argument (option value) must be a string or a number.");
  }
  /* check proper number of output arguments */
  if (nlhs > 0) {
    mexErrMsgTxt("No output arguments.");
  }
  optname = mxArrayToString(prhs[1]);

  mxopt = mxGetField(prhs[0],0,"opt");

  opt.MaxIter 	               = (typeInt)mxGetScalar(mxGetField(mxopt,0,"MaxIter"));
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
  mxGetString(mxGetField(mxopt, 0, "CostIntegrator"),opt.CostIntegrator,  VALUE_COSTINTMETHOD);
  mxGetString(mxGetField(mxopt,0,"Integrator"),opt.Integrator,VALUE_INTEGRATOR);
  mxGetString(mxGetField(mxopt,0,"LineSearchType"),opt.LineSearchType,VALUE_LSTYPE);
  mxGetString(mxGetField(mxopt,0,"JacobianX"),opt.JacobianX,VALUE_JACOBIANX);
  mxGetString(mxGetField(mxopt,0,"JacobianU"),opt.JacobianU,VALUE_JACOBIANU);
  mxGetString(mxGetField(mxopt,0,"IntegralCost"),opt.IntegralCost,VALUE_ONOFF);
  mxGetString(mxGetField(mxopt,0,"FinalCost"),opt.FinalCost,VALUE_ONOFF);

  /* rws structure (only necessary if MaxIter is modified) */
  mxrws = mxGetField(prhs[0],0,"rws");
  if (!strncmp(optname,"MaxIter",NAME_MAXITER)) {
    rws.lsAdapt = NULL;
	//rws.lsAdapt = mxGetPr(mxGetField(mxrws, 0, "lsAdapt"));
  }
  if (!strncmp(optname,"LineSearchInit",NAME_LSINIT) ||
      !strncmp(optname,"LineSearchIntervalFactor",NAME_LSINTFACTOR)) {
    rws.lsAdapt    = mxGetPr(mxGetField(mxrws,0,"lsAdapt"));
    rws.lsExplicit = mxGetPr(mxGetField(mxrws,0,"lsExplicit"));
  }

  grampc.opt = &opt;
  grampc.rws = &rws;

  if (mxIsChar(prhs[2])) {
	  optvalue = mxArrayToString(prhs[2]);
	  grampc_setopt_string(&grampc, optname, optvalue);
	  mxFree(optvalue);
  }
  else if (!strncmp(optname,"MaxIter",NAME_MAXITER)) {
    grampc_setopt_int(&grampc,optname,(typeInt)mxGetScalar(prhs[2]));
  }
  else {
    grampc_setopt_real(&grampc,optname,mxGetScalar(prhs[2]));
  }

  /* MaxIter */
  if (!strncmp(optname,"MaxIter",NAME_MAXITER)) {
    optOut = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    *((typeInt *)mxGetData(optOut)) = opt.MaxIter;
	rwsOut[0] = mxCreateDoubleMatrix(1,2*(NLS+1)*(1+opt.MaxIter),mxREAL);
    for (i = 0; i <= 2*(NLS+1)*(1+opt.MaxIter)-1; i++) {
      *(mxGetPr(rwsOut[0])+i) = rws.lsAdapt[i];
    }
	free(rws.lsAdapt);
	mxDestroyArray(mxGetField(mxrws, 0, "lsAdapt"));
    mxSetField(mxrws,0,"lsAdapt",rwsOut[0]);
  }
  /* ShiftControl */
  else if (!strncmp(optname,"ShiftControl",NAME_SHIFTCONTROL)) {
    optOut = mxCreateString(opt.ShiftControl);
  }
  /* ScaleProblem */
  else if (!strncmp(optname,"ScaleProblem",NAME_SCALEPROBLEM)) {
    optOut = mxCreateString(opt.ScaleProblem);
  }
  /* CostIntegrator */
  else if (!strncmp(optname,"CostIntegrator",NAME_COSTINTMETHOD)) {
    optOut = mxCreateString(opt.CostIntegrator);
  }
  /* Integrator */
  else if (!strncmp(optname,"Integrator",NAME_INTEGRATOR)) {
    optOut = mxCreateString(opt.Integrator);
  }
  /* IntegratorRelTol */
  else if (!strncmp(optname,"IntegratorRelTol",NAME_INTRELTOL)) {
    optOut = mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(optOut)) = opt.IntegratorRelTol;
  }
  /* IntegratorAbsTol */
  else if (!strncmp(optname,"IntegratorAbsTol",NAME_INTABSTOL)) {
    optOut = mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(optOut)) = opt.IntegratorAbsTol;
  }
  /* LineSearchType */
  else if (!strncmp(optname,"LineSearchType",NAME_LSTYPE)) {
    optOut = mxCreateString(opt.LineSearchType);
  }
  /* LineSearchMax */
  else if (!strncmp(optname,"LineSearchMax",NAME_LSMAX)) {
    optOut = mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(optOut)) = opt.LineSearchMax;
  }
  /* LineSearchMin */
  else if (!strncmp(optname,"LineSearchMin",NAME_LSMIN)) {
    optOut = mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(optOut)) = opt.LineSearchMin;
  }
  /* LineSearchInit */
  else if (!strncmp(optname,"LineSearchInit",NAME_LSINIT)) {
    optOut = mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(optOut)) = opt.LineSearchInit;
  }
  /* LineSearchIntervalFactor */
  else if (!strncmp(optname,"LineSearchIntervalFactor",NAME_LSINTFACTOR)) {
    optOut = mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(optOut)) = opt.LineSearchIntervalFactor;
  }
  /* LineSearchAdaptFactor */
  else if (!strncmp(optname,"LineSearchAdaptFactor",NAME_LSADAPTFACTOR)) {
    optOut = mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(optOut)) = opt.LineSearchAdaptFactor;
  }
  /* LineSearchIntervalTol */
  else if (!strncmp(optname,"LineSearchIntervalTol",NAME_LSINTTOL)) {
    optOut = mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(optOut)) = opt.LineSearchIntervalTol;
  }
  /* JacobianX */
  else if (!strncmp(optname,"JacobianX",NAME_JACOBIANX)) {
    optOut = mxCreateString(opt.JacobianX);
  }
  /* JacobianU */
  else if (!strncmp(optname,"JacobianU",NAME_JACOBIANU)) {
    optOut = mxCreateString(opt.JacobianU);
  }
  /* IntegralCost */
  else if (!strncmp(optname,"IntegralCost",NAME_INTEGRALCOST)) {
    optOut = mxCreateString(opt.IntegralCost);
  }
  /* FinalCost */
  else if (!strncmp(optname,"FinalCost",NAME_FINALCOST)) {
    optOut = mxCreateString(opt.FinalCost);
  }
  /* Set field */
  mxDestroyArray(mxGetField(mxopt, 0, optname));
  mxSetField(mxopt,0,optname,optOut);

  mxFree(optname);
}
