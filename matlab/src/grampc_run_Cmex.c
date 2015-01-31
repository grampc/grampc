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
 * File: grampc_run_Cmex.c
 * Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Date: February 2014
 * Version: v1.0
 * 
 * Mex function interface to use the function grampc_run 
 * for running GRAMPC.
 *
 */

#include "mex.h"
#include "matrix.h"


#include "grampc.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  typeGRAMPC grampc;
  typeGRAMPCparam param;
  typeGRAMPCsol sol;
  typeGRAMPCrws rws;
  typeGRAMPCopt opt;

  mxArray *mxopt   = NULL;
  mxArray *mxparam = NULL;
  mxArray *mxrws   = NULL;

  /* check proper number of input arguments */
  if (nrhs < 1) {
    mexErrMsgTxt("Not enough input arguments.");
  }
  if (nrhs > 1) {
    mexErrMsgTxt("Too many input arguments.");
  }
  /* first argument: param structure */
  if (!mxIsStruct(prhs[0])) {
    mexErrMsgTxt("First argument (grampc structure) must be a struct.");
  }

  /* check proper number of output arguments */
  if (nlhs > 3) {
    mexErrMsgTxt("Too many output arguments.");
  }

  mxopt   = mxGetField(prhs[0],0,"opt");
  mxparam = mxGetField(prhs[0],0,"param");
  mxrws   = mxGetField(prhs[0],0,"rws");

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
  param.Thor = mxGetScalar(mxGetField(mxparam,0,"Thor"));
  param.dt   = mxGetScalar(mxGetField(mxparam,0,"dt"));
  param.tk   = mxGetScalar(mxGetField(mxparam,0,"tk"));
  param.Nhor = (typeInt)mxGetScalar(mxGetField(mxparam,0,"Nhor"));
  if (mxGetField(mxparam,0,"pCost") == NULL) {
    param.pCost = NULL;
  }
  else {
    param.pCost = mxGetPr(mxGetField(mxparam,0,"pCost"));
  }
  if (mxGetField(mxparam,0,"pSys") == NULL) {
    param.pSys = NULL;
  }
  else {
    param.pSys = mxGetPr(mxGetField(mxparam,0,"pSys"));
  }
  param.umax    = mxGetPr(mxGetField(mxparam,0,"umax"));
  param.umin    = mxGetPr(mxGetField(mxparam,0,"umin"));
  param.xScale  = mxGetPr(mxGetField(mxparam,0,"xScale"));
  param.xOffset = mxGetPr(mxGetField(mxparam,0,"xOffset"));
  param.uScale  = mxGetPr(mxGetField(mxparam,0,"uScale"));
  param.uOffset = mxGetPr(mxGetField(mxparam,0,"uOffset"));
  param.NpCost  = (typeInt)mxGetScalar(mxGetField(mxparam,0,"NpCost"));
  param.NpSys   = (typeInt)mxGetScalar(mxGetField(mxparam,0,"NpSys"));

  rws.t                  = mxGetPr(mxGetField(mxrws,0,"t"));
  rws.x                  = mxGetPr(mxGetField(mxrws,0,"x"));
  rws.adj                = mxGetPr(mxGetField(mxrws,0,"adj"));
  rws.u                  = mxGetPr(mxGetField(mxrws,0,"u"));
  rws.dHdu               = mxGetPr(mxGetField(mxrws,0,"dHdu"));
  rws.lsAdapt            = mxGetPr(mxGetField(mxrws,0,"lsAdapt"));
  rws.uls                = mxGetPr(mxGetField(mxrws,0,"uls"));
  rws.lsExplicit         = mxGetPr(mxGetField(mxrws,0,"lsExplicit"));
  rws.uprev              = mxGetPr(mxGetField(mxrws,0,"uprev"));
  rws.dHduprev           = mxGetPr(mxGetField(mxrws,0,"dHduprev"));
  rws.J                  = mxGetPr(mxGetField(mxrws,0,"J"));
  rws.rwsScale           = mxGetPr(mxGetField(mxrws,0,"rwsScale"));
  rws.rwsGradient        = mxGetPr(mxGetField(mxrws,0,"rwsGradient"));
  rws.rwsCostIntegration = mxGetPr(mxGetField(mxrws,0,"rwsCostIntegration"));
  rws.rwsAdjIntegration  = mxGetPr(mxGetField(mxrws,0,"rwsAdjIntegration"));
  rws.rwsIntegration     = mxGetPr(mxGetField(mxrws,0,"rwsIntegration"));

  opt.MaxIter                  = (typeInt)mxGetScalar(mxGetField(mxopt,0,"MaxIter"));
  opt.IntegratorRelTol         = mxGetScalar(mxGetField(mxopt,0,"IntegratorRelTol"));
  opt.IntegratorAbsTol         = mxGetScalar(mxGetField(mxopt,0,"IntegratorAbsTol"));
  opt.LineSearchMax            = mxGetScalar(mxGetField(mxopt,0,"LineSearchMax"));
  opt.LineSearchMin            = mxGetScalar(mxGetField(mxopt,0,"LineSearchMin"));
  opt.LineSearchInit           = mxGetScalar(mxGetField(mxopt,0,"LineSearchInit"));
  opt.LineSearchIntervalFactor = mxGetScalar(mxGetField(mxopt,0,"LineSearchIntervalFactor"));
  opt.LineSearchAdaptFactor    = mxGetScalar(mxGetField(mxopt,0,"LineSearchAdaptFactor"));
  opt.LineSearchIntervalTol    = mxGetScalar(mxGetField(mxopt,0,"LineSearchIntervalTol"));

  strcpy(opt.ShiftControl,mxArrayToString(mxGetField(mxopt,0,"ShiftControl")));
  strcpy(opt.ScaleProblem,mxArrayToString(mxGetField(mxopt,0,"ScaleProblem")));
  strcpy(opt.CostIntegrator,mxArrayToString(mxGetField(mxopt,0,"CostIntegrator")));
  strcpy(opt.Integrator,mxArrayToString(mxGetField(mxopt,0,"Integrator")));
  strcpy(opt.LineSearchType,mxArrayToString(mxGetField(mxopt,0,"LineSearchType")));
  strcpy(opt.JacobianX,mxArrayToString(mxGetField(mxopt,0,"JacobianX")));
  strcpy(opt.JacobianU,mxArrayToString(mxGetField(mxopt,0,"JacobianU")));
  strcpy(opt.IntegralCost,mxArrayToString(mxGetField(mxopt,0,"IntegralCost")));
  strcpy(opt.FinalCost,mxArrayToString(mxGetField(mxopt,0,"FinalCost")));

  plhs[0] = mxCreateDoubleMatrix(param.Nx,1,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(param.Nu,1,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);

  sol.xnext = mxGetPr(plhs[0]);
  sol.unext = mxGetPr(plhs[1]);
  sol.J     = mxGetPr(plhs[2]);

  grampc.param = &param;
  grampc.opt   = &opt;
  grampc.rws   = &rws;
  grampc.sol   = &sol;

  grampc_run(&grampc);
}
