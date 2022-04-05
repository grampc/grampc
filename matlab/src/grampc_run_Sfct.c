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

#define S_FUNCTION_NAME  grampc_run_Sfct
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

/*include functions for gradient algorithm*/
#include "grampc.h"

void sdata2typeGRAMPC(const typeGRAMPC *grampc, ctypeInt *intoptidx, ctypeInt *intopt, ctypeInt *numoptidx,
	ctypeRNum *numopt, ctypeInt *numparamidx, ctypeRNum *numparam) {
	/* The option and parameter sequences must correspond to the sequence in typeGRAMPC2mx */
	typeInt i;

	/* vector parameters */
	i = 0;
	grampc_setparam_real_vector(grampc, "x0", (numparam + numparamidx[i])); i++;
	grampc_setparam_real_vector(grampc, "xdes", (numparam + numparamidx[i])); i++;

	grampc_setparam_real_vector(grampc, "u0", (numparam + numparamidx[i])); i++;
	grampc_setparam_real_vector(grampc, "udes", (numparam + numparamidx[i])); i++;
	grampc_setparam_real_vector(grampc, "umax", (numparam + numparamidx[i])); i++;
	grampc_setparam_real_vector(grampc, "umin", (numparam + numparamidx[i])); i++;

	grampc_setparam_real_vector(grampc, "p0", (numparam + numparamidx[i])); i++;
	grampc_setparam_real_vector(grampc, "pmax", (numparam + numparamidx[i])); i++;
	grampc_setparam_real_vector(grampc, "pmin", (numparam + numparamidx[i])); i++;

	grampc_setparam_real(grampc, "Thor", numparam[numparamidx[i]]); i++;
	grampc_setparam_real(grampc, "Tmax", numparam[numparamidx[i]]); i++;
	grampc_setparam_real(grampc, "Tmin", numparam[numparamidx[i]]); i++;

	grampc_setparam_real(grampc, "dt", numparam[numparamidx[i]]); i++;
	grampc_setparam_real(grampc, "t0", numparam[numparamidx[i]]); i++;


	/* Integer options */
	i = 0;
	grampc_setopt_int(grampc, "Nhor", intopt[intoptidx[i]]); i++;
	grampc_setopt_int(grampc, "MaxGradIter", intopt[intoptidx[i]]); i++;
	grampc_setopt_int(grampc, "MaxMultIter", intopt[intoptidx[i]]); i++;
	grampc_setopt_string(grampc, "ShiftControl", intopt[intoptidx[i]] == INT_ON ? "on" : "off"); i++;

	grampc_setopt_string(grampc, "TimeDiscretization", intopt[intoptidx[i]] == INT_UNIFORM ? "uniform" : "nonuniform"); i++;

	grampc_setopt_string(grampc, "IntegralCost", intopt[intoptidx[i]] == INT_ON ? "on" : "off"); i++;
	grampc_setopt_string(grampc, "TerminalCost", intopt[intoptidx[i]] == INT_ON ? "on" : "off"); i++;
	grampc_setopt_string(grampc, "IntegratorCost", intopt[intoptidx[i]] == INT_TRAPZ ? "trapezodial" : "simpson"); i++;

	grampc_setopt_string(grampc, "Integrator", IntegratorInt2Str(intopt[intoptidx[i]])); i++;
	grampc_setopt_int(grampc, "IntegratorMaxSteps", intopt[intoptidx[i]]); i++;
	grampc_setopt_int_vector(grampc, "FlagsRodas", (intopt + intoptidx[i])); i++;

	grampc_setopt_string(grampc, "LineSearchType", LineSearchTypeInt2Str(intopt[intoptidx[i]])); i++;
	grampc_setopt_string(grampc, "LineSearchExpAutoFallback", intopt[intoptidx[i]] == INT_ON ? "on" : "off"); i++;

	grampc_setopt_string(grampc, "OptimControl", intopt[intoptidx[i]] == INT_ON ? "on" : "off"); i++;
	grampc_setopt_string(grampc, "OptimParam", intopt[intoptidx[i]] == INT_ON ? "on" : "off"); i++;
	grampc_setopt_string(grampc, "OptimTime", intopt[intoptidx[i]] == INT_ON ? "on" : "off"); i++;

	grampc_setopt_string(grampc, "ScaleProblem", intopt[intoptidx[i]] == INT_ON ? "on" : "off"); i++;

	grampc_setopt_string(grampc, "EqualityConstraints", intopt[intoptidx[i]] == INT_ON ? "on" : "off"); i++;
	grampc_setopt_string(grampc, "InequalityConstraints", intopt[intoptidx[i]] == INT_ON ? "on" : "off"); i++;
	grampc_setopt_string(grampc, "TerminalEqualityConstraints", intopt[intoptidx[i]] == INT_ON ? "on" : "off"); i++;
	grampc_setopt_string(grampc, "TerminalInequalityConstraints", intopt[intoptidx[i]] == INT_ON ? "on" : "off"); i++;
	grampc_setopt_string(grampc, "ConstraintsHandling", intopt[intoptidx[i]] == INT_EXTPEN ? "extpen" : "auglag"); i++;

	grampc_setopt_string(grampc, "ConvergenceCheck", intopt[intoptidx[i]] == INT_ON ? "on" : "off"); i++;

	/* real options */
	i = 0;
	grampc_setopt_real(grampc, "IntegratorRelTol", numopt[numoptidx[i]]); i++;
	grampc_setopt_real(grampc, "IntegratorAbsTol", numopt[numoptidx[i]]); i++;
	grampc_setopt_real(grampc, "IntegratorMinStepSize", numopt[numoptidx[i]]); i++;

	grampc_setopt_real(grampc, "LineSearchMax", numopt[numoptidx[i]]); i++;
	grampc_setopt_real(grampc, "LineSearchMin", numopt[numoptidx[i]]); i++;
	grampc_setopt_real(grampc, "LineSearchInit", numopt[numoptidx[i]]); i++;
	grampc_setopt_real(grampc, "LineSearchAdaptAbsTol", numopt[numoptidx[i]]); i++;
	grampc_setopt_real(grampc, "LineSearchAdaptFactor", numopt[numoptidx[i]]); i++;
	grampc_setopt_real(grampc, "LineSearchIntervalTol", numopt[numoptidx[i]]); i++;
	grampc_setopt_real(grampc, "LineSearchIntervalFactor", numopt[numoptidx[i]]); i++;

	grampc_setopt_real(grampc, "OptimParamLineSearchFactor", numopt[numoptidx[i]]); i++;
	grampc_setopt_real(grampc, "OptimTimeLineSearchFactor", numopt[numoptidx[i]]); i++;

	grampc_setopt_real_vector(grampc, "xScale", (numopt + numoptidx[i])); i++;
	grampc_setopt_real_vector(grampc, "xOffset", (numopt + numoptidx[i])); i++;
	grampc_setopt_real_vector(grampc, "uScale", (numopt + numoptidx[i])); i++;
	grampc_setopt_real_vector(grampc, "uOffset", (numopt + numoptidx[i])); i++;
	grampc_setopt_real_vector(grampc, "pScale", (numopt + numoptidx[i])); i++;
	grampc_setopt_real_vector(grampc, "pOffset", (numopt + numoptidx[i])); i++;
	grampc_setopt_real(grampc, "TScale", numopt[numoptidx[i]]); i++;
	grampc_setopt_real(grampc, "TOffset", numopt[numoptidx[i]]); i++;
	grampc_setopt_real(grampc, "JScale", numopt[numoptidx[i]]); i++;
	grampc_setopt_real_vector(grampc, "cScale", (numopt + numoptidx[i])); i++;

	grampc_setopt_real_vector(grampc, "ConstraintsAbsTol", (numopt + numoptidx[i])); i++;

	grampc_setopt_real(grampc, "MultiplierMax", numopt[numoptidx[i]]); i++;
	grampc_setopt_real(grampc, "MultiplierDampingFactor", numopt[numoptidx[i]]); i++;
	grampc_setopt_real(grampc, "PenaltyMax", numopt[numoptidx[i]]); i++;
	grampc_setopt_real(grampc, "PenaltyMin", numopt[numoptidx[i]]); i++;
	grampc_setopt_real(grampc, "PenaltyIncreaseFactor", numopt[numoptidx[i]]); i++;
	grampc_setopt_real(grampc, "PenaltyDecreaseFactor", numopt[numoptidx[i]]); i++;
	grampc_setopt_real(grampc, "PenaltyIncreaseThreshold", numopt[numoptidx[i]]); i++;
	grampc_setopt_real(grampc, "AugLagUpdateGradientRelTol", numopt[numoptidx[i]]); i++;

	grampc_setopt_real(grampc, "ConvergenceGradientRelTol", numopt[numoptidx[i]]); i++;
}

static void mdlInitializeSizes(SimStruct *S)
{
	int_T  i;                                   	/* loop variable */

	typeInt Nx, Nu, Np, Ng, Nh, NgT, NhT;

	ocp_dim(&Nx, &Nu, &Np, &Ng, &Nh, &NgT, &NhT, NULL);

	/*Parameters */
	ssSetNumSFcnParams(S, 6);
	/* parameter tuning not allowed during simulation */
	for (i = 0; i < ssGetSFcnParamsCount(S); i++) {
		ssSetSFcnParamTunable(S, i, SS_PRM_NOT_TUNABLE);
	}

	/*Input ports*/
	if (!ssSetNumInputPorts(S, 5)) return;        /* number of input ports */
	ssSetInputPortWidth(S, 0, Nx);                /* no. of inputs at input port 0 (x0)        */
	ssSetInputPortWidth(S, 1, 1);                 /* no. of inputs at input port 1 (t0)        */
	ssSetInputPortWidth(S, 2, Nx);                /* no. of inputs at input port 2 (xdes)      */
	ssSetInputPortWidth(S, 3, MAX(Nu, 1));        /* no. of inputs at input port 3 (udes)      */
	ssSetInputPortWidth(S, 4, DYNAMICALLY_SIZED); /* no. of inputs at input port 4 (userparam) */

	for (i = 0; i < 5; i++) {
		ssSetInputPortDirectFeedThrough(S, i, 1);   /* direct input feedthrough */
		ssSetInputPortRequiredContiguous(S, i, 1);  /* must be contiguous to use  ssGetInputPortRealSignal*/
		ssSetInputPortDataType(S, i, SS_TYPERNUM);
	}

	/*Output ports*/
	if (!ssSetNumOutputPorts(S, 8)) return;  /* number of output ports                     */
	ssSetOutputPortWidth(S, 0, Nx);          /* no. of outputs at output port 0 (xnext)    */
	ssSetOutputPortWidth(S, 1, MAX(Nu, 1));  /* no. of outputs at output port 2 (unext)    */
	ssSetOutputPortWidth(S, 2, MAX(Np, 1));  /* no. of outputs at output port 3 (pnext)    */
	ssSetOutputPortWidth(S, 3, 1);           /* no. of outputs at output port 4 (Tnext)    */
	ssSetOutputPortWidth(S, 4, 2);           /* no. of outputs at output port 5 (J)        */
	ssSetOutputPortWidth(S, 5, 1);           /* no. of outputs at output port 6 (||cfct||) */
	ssSetOutputPortWidth(S, 6, 1);           /* no. of outputs at output port 7 (||pen||)  */
	ssSetOutputPortWidth(S, 7, 1);           /* no. of outputs at output port 8 (status)   */
	/*ssSetOutputPortWidth(S, 8, 1);*/           /* no. of outputs at output port 8 (iter)     */

	for (i = 0; i < 7; i++) {
		ssSetOutputPortDataType(S, i, SS_TYPERNUM);
	}
	ssSetOutputPortDataType(S, 7, SS_INT32);

	/*Other*/
	ssSetNumSampleTimes(S, 0);    /* number of sample times         */
	ssSetNumRWork(S, 0);          /* number of real work vectors    */
	ssSetNumIWork(S, 0);          /* number of integer work vectors */
	ssSetNumPWork(S, 1);          /* number of pointer work vectors */

	ssSetNumModes(S, 0);          /* number of modes to switch between    */
	ssSetNumNonsampledZCs(S, 0);  /* number of non-sampled zero-crossings */
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
	ssSetInputPortWidth(S, 4, 1); /* default to one */
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
	typeGRAMPC *grampc;
	typeUSERPARAM *userparam;
	size_t size_userpram;
	typeInt *intoptidx, *numoptidx, *numparamidx, *intopt, i;
	typeRNum *numopt, *numparam, *real_userparam;
	void **ptr;

	const double *d_intoptidx = mxGetData(ssGetSFcnParam(S, 0));
	const double *d_intopt = mxGetData(ssGetSFcnParam(S, 1));
	const double *d_numoptidx = mxGetData(ssGetSFcnParam(S, 2));
	const double *d_numopt = mxGetData(ssGetSFcnParam(S, 3));
	const double *d_numparamidx = mxGetData(ssGetSFcnParam(S, 4));
	const double *d_numparam = mxGetData(ssGetSFcnParam(S, 5));

	ptr = ssGetPWork(S);

	/* Memory allocation for userparam, use size of input port */
	const typeRNum *d_userparam = (typeRNum *)ssGetInputPortRealSignal(S, 4);
	size_userpram = ssGetInputPortWidth(S, 4);
	userparam = calloc(size_userpram * sizeof(typeRNum), 1);
	real_userparam = (typeRNum*)userparam;
	for (i = 0; i < size_userpram; i++) {
		real_userparam[i] = (typeRNum)d_userparam[i];
	}

	/* Memory allocation and initialization of grampc */
	grampc_init(&grampc, userparam);

	/* cast parameters and options to the correct data type */
	CastDvec2Intvec(&intoptidx, d_intoptidx + 1, (size_t)d_intoptidx[0]);
	CastDvec2Intvec(&numoptidx, d_numoptidx + 1, (size_t)d_numoptidx[0]);
	CastDvec2Intvec(&numparamidx, d_numparamidx + 1, (size_t)d_numparamidx[0]);
	CastDvec2Intvec(&intopt, d_intopt + 1, (size_t)d_intopt[0]);
	CastDvec2Numvec(&numopt, d_numopt + 1, (size_t)d_numopt[0]);
	CastDvec2Numvec(&numparam, d_numparam + 1, (size_t)d_numparam[0]);

	/* Update the params and options */
	sdata2typeGRAMPC(grampc, intoptidx, intopt, numoptidx, numopt, numparamidx, numparam);

	/* free temporary used memory */
	free(intoptidx);
	free(numoptidx);
	free(numparamidx);
	free(intopt);
	free(numopt);
	free(numparam);

	/* print options and parameters */
	grampc_printopt(grampc);
	grampc_printparam(grampc);

	ptr[0] = (void *)grampc;
}
#endif   /*MDL_START*/


/***************************************************************************
 * mdlOutputs                                                              *
 ***************************************************************************/

static void mdlOutputs(SimStruct *S, int_T tid)
{
	int_T  i;

	/* inputs */
	ctypeRNum *x0 = (typeRNum *)ssGetInputPortRealSignal(S, 0);
	ctypeRNum *t0 = (typeRNum *)ssGetInputPortRealSignal(S, 1);
	ctypeRNum *xdes = (typeRNum *)ssGetInputPortRealSignal(S, 2);
	ctypeRNum *udes = (typeRNum *)ssGetInputPortRealSignal(S, 3);

	/* outputs */
	typeRNum *xnext = (typeRNum *)ssGetOutputPortRealSignal(S, 0);
	typeRNum *unext = (typeRNum *)ssGetOutputPortRealSignal(S, 1);
	typeRNum *pnext = (typeRNum *)ssGetOutputPortRealSignal(S, 2);
	typeRNum *Tnext = (typeRNum *)ssGetOutputPortRealSignal(S, 3);
	typeRNum *cost = (typeRNum *)ssGetOutputPortRealSignal(S, 4);
	typeRNum *cfct = (typeRNum *)ssGetOutputPortRealSignal(S, 5);
	typeRNum *pen = (typeRNum *)ssGetOutputPortRealSignal(S, 6);
	typeInt *status = (typeInt *)ssGetOutputPortRealSignal(S, 7);

	/* workspace */
	const typeGRAMPC *grampc = (typeGRAMPC *)ssGetPWorkValue(S, 0);

	/* remove warning unreferenced formal parameter */
	(void)(tid);

	/* set the current quantitites */
	grampc_setparam_real_vector(grampc, "x0", x0);
	grampc_setparam_real(grampc, "t0", t0[0]);
	grampc_setparam_real_vector(grampc, "xdes", xdes);
	grampc_setparam_real_vector(grampc, "udes", udes);
	
	/* update userparm */
	typeUSERPARAM *userparam = grampc->userparam;
	const typeRNum *d_userparam = (typeRNum *)ssGetInputPortRealSignal(S, 4);
	size_t size_userpram = ssGetInputPortWidth(S, 4);
	typeRNum *loc_userparam = (typeRNum *)userparam;
	for (i = 0; i < size_userpram; i++) {
		loc_userparam[i] = (typeRNum)d_userparam[i];
	}

	/* run mpc */
	grampc_run(grampc);

	/* Copy the solution */
	MatCopy(xnext, grampc->sol->xnext, 1, grampc->param->Nx);
	MatCopy(unext, grampc->sol->unext, 1, grampc->param->Nu);
	MatCopy(pnext, grampc->sol->pnext, 1, grampc->param->Np);
	Tnext[0] = grampc->sol->Tnext;
	MatCopy(cost, grampc->sol->J, 1, 2);
	cfct[0] = grampc->sol->cfct;
	pen[0] = grampc->sol->pen;
	status[0] = grampc->sol->status;
}


/****************************************************************************
 * mdlTerminate                                                             *
 ****************************************************************************/
static void mdlTerminate(SimStruct *S)
{
	typeGRAMPC *grampc = (typeGRAMPC *)ssGetPWorkValue(S, 0);

	/* free allocated memory */
	free(grampc->userparam);
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
	* End of C-Code S-Function	                                               *
	****************************************************************************/
