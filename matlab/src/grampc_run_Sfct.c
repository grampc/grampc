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
 * File: grampc_run_Sfct.c
 * Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Date: February 2014
 * Version: v1.0
 * 
 * Level 2 S-function to use GRAMPC in MATLAB/Simulink.
 * The function requires a variable of datatype typeGRAMPC
 * as input parameter for the MATLAB/Simulink block.
 *
 */

#define S_FUNCTION_NAME  grampc_run_Sfct
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include <math.h>

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

/*include functions for gradient algorithm*/
#include "grampc.h"

static void mdlInitializeSizes(SimStruct *S)
{
  int_T  i;                                   	/* loop variable */

  typeInt Nx, Nu;

  sysdim(&Nx,&Nu);

  ssSetNumSFcnParams(S, 1);			

  /* parameter tuning not allowed during simulation */
  for (i=0; i<ssGetSFcnParamsCount(S); i++) {
    ssSetSFcnParamTunable(S, i, 1);
  }

  /*Input ports*/
  if(!ssSetNumInputPorts(S, 4)) return;       	/* number of input ports */

  ssSetInputPortWidth(S, 0, Nx);		/* no. of inputs at input port 0 (xk)   */
  ssSetInputPortWidth(S, 1,  1);		/* no. of inputs at input port 1 (tk)   */
  ssSetInputPortWidth(S, 2, Nx);	        /* no. of inputs at input port 2 (xdes) */
  ssSetInputPortWidth(S, 3, Nu);                /* no. of inputs at input port 3 (udes) */

  ssSetInputPortDirectFeedThrough(S, 0, 1);    	/* direct input feedthrough */
  ssSetInputPortDirectFeedThrough(S, 1, 1);    	/* direct input feedthrough */
  ssSetInputPortDirectFeedThrough(S, 2, 1);    	/* direct input feedthrough */
  ssSetInputPortDirectFeedThrough(S, 3, 1);    	/* direct input feedthrough */

  /*Output ports*/
  if(!ssSetNumOutputPorts(S, 3)) return;      	/* number of output ports                  */
  ssSetOutputPortWidth(S, 0, Nx);     	      	/* no. of outputs at output port 0 (xnext) */
  ssSetOutputPortWidth(S, 1, Nu);		/* no. of outputs at output port 1 (unext) */
  ssSetOutputPortWidth(S, 2, 1);             	/* no. of outputs at output port 2 (J)     */

  /*Other*/
  ssSetNumSampleTimes(S, 0);   	             	/* number of sample times         */
  ssSetNumRWork(      S, 0);	    	        /* number of real work vectors    */
  ssSetNumIWork(      S, 0);                	/* number of integer work vectors */
  ssSetNumPWork(      S, 1);                	/* number of pointer work vectors */

  ssSetNumModes(        S, 0);                 	/* number of modes to switch between    */
  ssSetNumNonsampledZCs(S, 0);                 	/* number of non-sampled zero-crossings */
}


/***************************************************************************
 * mdlInitializeSampleTimes                                                *
 ***************************************************************************/

static void mdlInitializeSampleTimes(SimStruct *S)
{
  mxArray *mxparam = mxGetField(ssGetSFcnParam(S,0),0,"param");
  ssSetSampleTime(S, 0, mxGetScalar(mxGetField(mxparam,0,"dt")));
  ssSetOffsetTime(S, 0, 0.0);
}


/***************************************************************************
 * mdlStart                                                                *
 ***************************************************************************/

#define MDL_START 
#if defined (MDL_START)
static void mdlStart(SimStruct *S)
{
  typeGRAMPC *grampc;
  void **ptr;
  mxArray *mxopt   = mxGetField(ssGetSFcnParam(S,0),0,"opt");
  mxArray *mxparam = mxGetField(ssGetSFcnParam(S,0),0,"param");

  ptr = ssGetPWork(S);

  /* Memory allocation and initialization */
  grampc_init(&grampc);

  /* adapt option structure */
  grampc_setopt_int(grampc,"MaxIter",(typeInt)mxGetScalar(mxGetField(mxopt,0,"MaxIter")));
  grampc_setopt_string(grampc,"ShiftControl",mxArrayToString(mxGetField(mxopt,0,"ShiftControl")));
  grampc_setopt_string(grampc,"ScaleProblem",mxArrayToString(mxGetField(mxopt,0,"ScaleProblem")));
  grampc_setopt_string(grampc,"CostIntegrator",mxArrayToString(mxGetField(mxopt,0,"CostIntegrator")));
  grampc_setopt_string(grampc,"Integrator",mxArrayToString(mxGetField(mxopt,0,"Integrator")));  
  grampc_setopt_real(grampc,"IntegratorRelTol",mxGetScalar(mxGetField(mxopt,0,"IntegratorRelTol")));
  grampc_setopt_real(grampc,"IntegratorAbsTol",mxGetScalar(mxGetField(mxopt,0,"IntegratorAbsTol")));
  grampc_setopt_string(grampc,"LineSearchType",mxArrayToString(mxGetField(mxopt,0,"LineSearchType")));  
  grampc_setopt_real(grampc,"LineSearchMax",mxGetScalar(mxGetField(mxopt,0,"LineSearchMax")));
  grampc_setopt_real(grampc,"LineSearchMin",mxGetScalar(mxGetField(mxopt,0,"LineSearchMin")));
  grampc_setopt_real(grampc,"LineSearchInit",mxGetScalar(mxGetField(mxopt,0,"LineSearchInit")));
  grampc_setopt_real(grampc,"LineSearchIntervalFactor",mxGetScalar(mxGetField(mxopt,0,"LineSearchIntervalFactor")));
  grampc_setopt_real(grampc,"LineSearchAdaptFactor",mxGetScalar(mxGetField(mxopt,0,"LineSearchAdaptFactor")));
  grampc_setopt_real(grampc,"LineSearchIntervalTol",mxGetScalar(mxGetField(mxopt,0,"LineSearchIntervalTol")));
  grampc_setopt_string(grampc,"JacobianX",mxArrayToString(mxGetField(mxopt,0,"JacobianX")));
  grampc_setopt_string(grampc,"JacobianU",mxArrayToString(mxGetField(mxopt,0,"JacobianU")));
  grampc_setopt_string(grampc,"IntegralCost",mxArrayToString(mxGetField(mxopt,0,"IntegralCost")));
  grampc_setopt_string(grampc,"FinalCost",mxArrayToString(mxGetField(mxopt,0,"FinalCost")));

  /* adapt param structure */
  grampc_setparam_int(grampc,"Nhor",(typeInt)mxGetScalar(mxGetField(mxparam,0,"Nhor")));
  grampc_setparam_vector(grampc,"umax",mxGetPr(mxGetField(mxparam,0,"umax")));
  grampc_setparam_vector(grampc,"umin",mxGetPr(mxGetField(mxparam,0,"umin")));
  grampc_setparam_vector(grampc,"xScale",mxGetPr(mxGetField(mxparam,0,"xScale")));
  grampc_setparam_vector(grampc,"xOffset",mxGetPr(mxGetField(mxparam,0,"xOffset")));
  grampc_setparam_vector(grampc,"uScale",mxGetPr(mxGetField(mxparam,0,"uScale")));
  grampc_setparam_vector(grampc,"uOffset",mxGetPr(mxGetField(mxparam,0,"uOffset")));
  grampc_setparam_vector(grampc,"xk",mxGetPr(mxGetField(mxparam,0,"xk")));
  grampc_setparam_vector(grampc,"u0",mxGetPr(mxGetField(mxparam,0,"u0")));
  grampc_setparam_vector(grampc,"xdes",mxGetPr(mxGetField(mxparam,0,"xdes")));
  grampc_setparam_vector(grampc,"udes",mxGetPr(mxGetField(mxparam,0,"udes")));
  grampc_setparam_real(grampc,"Thor",mxGetScalar(mxGetField(mxparam,0,"Thor")));
  grampc_setparam_real(grampc,"dt",mxGetScalar(mxGetField(mxparam,0,"dt")));
  grampc_setparam_real(grampc,"tk",mxGetScalar(mxGetField(mxparam,0,"tk")));
  if ((typeInt)mxGetScalar(mxGetField(mxparam,0,"NpCost")) > 0) {
    grampc_setparam_int(grampc,"NpCost",(typeInt)mxGetScalar(mxGetField(mxparam,0,"NpCost")));
    grampc_setparam_vector(grampc,"pCost",mxGetPr(mxGetField(mxparam,0,"pCost")));
  }
  if ((typeInt)mxGetScalar(mxGetField(mxparam,0,"NpSys")) > 0) {
    grampc_setparam_int(grampc,"NpSys",(typeInt)mxGetScalar(mxGetField(mxparam,0,"NpSys")));
    grampc_setparam_vector(grampc,"pSys",mxGetPr(mxGetField(mxparam,0,"pSys")));
  }

  /* print options */
  grampc_printopt(grampc);

  /* print MPC parameter */
  grampc_printparam(grampc);

  ptr[0] = (void *)grampc;
}
#endif   /*MDL_START*/


/***************************************************************************
 * mdlOutputs = Enthält den Code für die Berechnung  der Ausgänge          *
 ***************************************************************************/

static void mdlOutputs(SimStruct *S, int_T tid)
{
  /*inputs*/
  typeRNum *xk   = (typeRNum *)ssGetInputPortRealSignalPtrs(S,0)[0];
  typeRNum tk    = (typeInt)ssGetInputPortRealSignalPtrs(S,1)[0];
  typeRNum *xdes = (typeRNum *)ssGetInputPortRealSignalPtrs(S,2)[0];
  typeRNum *udes = (typeRNum *)ssGetInputPortRealSignalPtrs(S,3)[0];

  /*outputs*/
  typeRNum *xnext  = (typeRNum *)ssGetOutputPortRealSignal(S,0);
  typeRNum *unext  = (typeRNum *)ssGetOutputPortRealSignal(S,1);
  typeRNum *cost   = (typeRNum *)ssGetOutputPortRealSignal(S,2);

  /*parameters*/
  typeInt i;

  typeGRAMPC *grampc = (typeGRAMPC *)ssGetPWorkValue(S,0);

  grampc_setparam_vector(grampc, "xk", xk);
  grampc_setparam_real(grampc, "tk", tk);
  grampc_setparam_vector(grampc, "xdes", xdes);
  grampc_setparam_vector(grampc, "udes", udes);

  grampc_run(grampc);

  for (i = 0; i <= grampc->param->Nx-1; i++) {
    xnext[i] = grampc->sol->xnext[i];
  }
  for (i = 0; i <= grampc->param->Nu-1; i++) {
    unext[i] = grampc->sol->unext[i];
  }
  cost[0] = grampc->sol->J[0];
}


/****************************************************************************
 * mdlTerminate                                                             *
 ****************************************************************************/
static void mdlTerminate(SimStruct *S)
{
  typeGRAMPC *grampc = (typeGRAMPC *)ssGetPWorkValue(S,0);

  /* free allocated memory */
  grampc_free(&grampc);
}


/****************************************************************************
 * Compiler dependent settings                                              *
 ****************************************************************************/

#ifdef	MATLAB_MEX_FILE 
#include "simulink.c"   
#else
#include "cg_sfun.h"    
#endif


/****************************************************************************
 * End of C-Code S-Function	                                            *
 ****************************************************************************/
