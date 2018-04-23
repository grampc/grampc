/* This file is part of GRAMPC - (https://sourceforge.net/projects/grampc/)
 *
 * GRAMPC -- A software framework for embedded nonlinear model predictive
 * control using a gradient-based augmented Lagrangian approach
 *
 * Copyright (C) 2014-2018 by Tobias Englert, Knut Graichen, Felix Mesmer,
 * Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 * Developed at the Institute of Measurement, Control, and Microtechnology,
 * Ulm University. All rights reserved.
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
 * License along with GRAMPC. If not, see <http://www.gnu.org/licenses/>
 *
 */

#define S_FUNCTION_NAME  grampc_ffct_Sfct
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

/* include system functions */
#include "probfct.h"
#include "grampc_init.h"

static void mdlInitializeSizes(SimStruct *S)
{
	int_T  i;                                   	/* loop variable */

	typeInt Nx, Nu, Np, Ng, Nh, NgT, NhT;
	ocp_dim(&Nx, &Nu, &Np, &Ng, &Nh, &NgT, &NhT, NULL);

	/* Parameters */
	ssSetNumSFcnParams(S, 1);
	/* parameter tuning not allowed during simulation */
	ssSetSFcnParamTunable(S, 0, SS_PRM_NOT_TUNABLE);

	/*Input ports*/
	if (!ssSetNumInputPorts(S, 4)) return;        /* number of input ports */
	ssSetInputPortWidth(S, 0, 1);                 /* no. of inputs at input port 0 (t0) */
	ssSetInputPortWidth(S, 1, MAX(Nu, 1));        /* no. of inputs at input port 1 (u) */
	ssSetInputPortWidth(S, 2, MAX(Np, 1));        /* no. of inputs at input port 2 (p) */
	ssSetInputPortWidth(S, 3, Nx);                /* no. of inputs at input port 3 (x0) */

	for (i = 0; i < 4; i++) {
		ssSetInputPortDirectFeedThrough(S, i, 1);   /* direct input feedthrough */
		ssSetInputPortRequiredContiguous(S, i, 1);  /* must be contiguous to use  ssGetInputPortRealSignal*/
		ssSetInputPortDataType(S, i, SS_TYPERNUM);
	}

	/*Output ports*/
	if (!ssSetNumOutputPorts(S, 1)) return;       /* number of output ports */
	ssSetOutputPortWidth(S, 0, Nx);               /* no. of outputs at output port 0 (dxdt) */

	ssSetOutputPortDataType(S, 0, SS_TYPERNUM);

	/*Other*/
	ssSetNumSampleTimes(S, 0);                    /* number of sample times         */
	ssSetNumRWork(S, 0);                          /* number of real work vectors    */
	ssSetNumIWork(S, 0);                          /* number of integer work vectors */
	ssSetNumPWork(S, 1);                          /* number of pointer work vectors */

	ssSetNumModes(S, 0);                          /* number of modes to switch between    */
	ssSetNumNonsampledZCs(S, 0);                  /* number of non-sampled zero-crossings */
}


/***************************************************************************
* mdlInitializeSampleTimes                                                *
***************************************************************************/

static void mdlInitializeSampleTimes(SimStruct *S)
{
	/* Attention: Read dt */
	ssSetSampleTime(S, 0, -1);
	ssSetOffsetTime(S, 0, 0.0);
}


/***************************************************************************
* mdlStart                                                                *
***************************************************************************/

#define MDL_START 
#if defined (MDL_START)
static void mdlStart(SimStruct *S)
{
	size_t size_userpram;
	int_T  i;
	typeUSERPARAM *userparam;
	typeRNum *real_userparam;
	void **ptr;

	const mxArray* const mxuserparam = ssGetSFcnParam(S, 0);
	const double *d_userparam = mxGetPr(mxuserparam);

	ptr = ssGetPWork(S);

	/* Memory allocation for userparam */
	size_userpram = mxGetM(mxuserparam)*mxGetN(mxuserparam);
	userparam = calloc(size_userpram * sizeof(typeRNum), 1);
	real_userparam = (typeRNum*)userparam;
	for (i = 0; i < size_userpram; i++) {
		real_userparam[i] = (typeRNum)d_userparam[i];
	}

	ptr[0] = (void*)real_userparam;

}
#endif   /*MDL_START*/


/***************************************************************************
* mdlOutputs                                                              *
***************************************************************************/

static void mdlOutputs(SimStruct *S, int_T tid)
{
	/* inputs */
	ctypeRNum *t = (typeRNum *)ssGetInputPortRealSignal(S, 0);
	ctypeRNum *u = (typeRNum *)ssGetInputPortRealSignal(S, 1);
	ctypeRNum *p = (typeRNum *)ssGetInputPortRealSignal(S, 2);
	ctypeRNum *x = (typeRNum *)ssGetInputPortRealSignal(S, 3);

	/* outputs */
	typeRNum *dxdt = (typeRNum *)ssGetOutputPortRealSignal(S, 0);

	/* workspace */
	typeUSERPARAM *userparam = (typeUSERPARAM *)ssGetPWorkValue(S, 0);

	/* remove warning unreferenced formal parameter */
	(void)(tid);

	/* set the current quantitites */
	ffct(dxdt, t[0], x, u, p, userparam);
}


/****************************************************************************
* mdlTerminate                                                             *
****************************************************************************/
static void mdlTerminate(SimStruct *S)
{
	free(ssGetPWorkValue(S, 0));
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
* End of C-Code S-Function	                                               *
****************************************************************************/
