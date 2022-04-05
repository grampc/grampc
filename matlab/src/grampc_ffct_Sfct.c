/* This file is part of GRAMPC - (https://sourceforge.net/projects/grampc/)
 *
 * GRAMPC -- A software framework for embedded nonlinear model predictive
 * control using a gradient-based augmented Lagrangian approach
 *
 * Copyright 2014-2019 by Tobias Englert, Knut Graichen, Felix Mesmer,
 * Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 * All rights reserved.
 *
 * GRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
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

	/* no Parameters */
	ssSetNumSFcnParams(S, 0);

	/*Input ports*/
	if (!ssSetNumInputPorts(S, 5)) return;        /* number of input ports */
	ssSetInputPortWidth(S, 0, 1);                 /* no. of inputs at input port 0 (t0) */
	ssSetInputPortWidth(S, 1, MAX(Nu, 1));        /* no. of inputs at input port 1 (u) */
	ssSetInputPortWidth(S, 2, MAX(Np, 1));        /* no. of inputs at input port 2 (p) */
	ssSetInputPortWidth(S, 3, Nx);                /* no. of inputs at input port 3 (x0) */
	ssSetInputPortWidth(S, 4, DYNAMICALLY_SIZED); /* no. of inputs at input port 4 (userparam) */

	for (i = 0; i < 5; i++) {
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


#if defined(MATLAB_MEX_FILE)
#define MDL_SET_INPUT_PORT_WIDTH
static void mdlSetInputPortWidth(SimStruct *S, int_T port, int_T inputPortWidth)
{
	if ((ssGetInputPortWidth(S, port) != DYNAMICALLY_SIZED) && (ssGetInputPortWidth(S, port) != inputPortWidth))
	{
		ssSetErrorStatus(S, "something wrong");
		return;
	}
	ssSetInputPortWidth(S, port, inputPortWidth);
}

#define MDL_SET_OUTPUT_PORT_WIDTH
static void mdlSetOutputPortWidth(SimStruct *S, int_T port, int_T outputPortWidth)
{
	/* We have no dynamic outputs */
}

#define MDL_SET_DEFAULT_PORT_DIMENSION_INFO
static void mdlSetDefaultPortDimensionInfo(SimStruct *S)
{
	ssSetInputPortWidth(S, 0, 1); /* default to one */
}

#define MDL_SET_WORK_WIDTHS
static void mdlSetWorkWidths(SimStruct *S)
{
	/* We have no states */
}
#endif


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
	
	ptr = ssGetPWork(S);

	/* Memory allocation for userparam, use size of input port */
	const typeRNum *d_userparam = (typeRNum *)ssGetInputPortRealSignal(S, 4);
	size_userpram = ssGetInputPortWidth(S, 4);
	userparam = calloc(size_userpram * sizeof(typeRNum), 1);
	real_userparam = (typeRNum*)userparam;
	for (i = 0; i < size_userpram; i++) {
		real_userparam[i] = (typeRNum)d_userparam[i];
	}

	/* set userparam as work vector */
	ptr[0] = (void*)real_userparam;
}
#endif   /*MDL_START*/


/***************************************************************************
* mdlOutputs                                                              *
***************************************************************************/

static void mdlOutputs(SimStruct *S, int_T tid)
{
	int_T  i;

	/* inputs */
	ctypeRNum *t = (typeRNum *)ssGetInputPortRealSignal(S, 0);
	ctypeRNum *u = (typeRNum *)ssGetInputPortRealSignal(S, 1);
	ctypeRNum *p = (typeRNum *)ssGetInputPortRealSignal(S, 2);
	ctypeRNum *x = (typeRNum *)ssGetInputPortRealSignal(S, 3);

	/* outputs */
	typeRNum *dxdt = (typeRNum *)ssGetOutputPortRealSignal(S, 0);

	/* workspace */
	typeUSERPARAM *userparam = (typeUSERPARAM *)ssGetPWorkValue(S, 0);
	
	/* update userparm */
	const typeRNum *d_userparam = (typeRNum *)ssGetInputPortRealSignal(S, 4);
	size_t size_userpram = ssGetInputPortWidth(S, 4);
	typeRNum *loc_userparam = (typeRNum *)userparam;
	for (i = 0; i < size_userpram; i++) {
		loc_userparam[i] = (typeRNum)d_userparam[i];
	}
	
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
