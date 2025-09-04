/* This file is part of GRAMPC - (https://github.com/grampc/grampc)
 *
 * GRAMPC -- A software framework for embedded nonlinear model predictive
 * control using a gradient-based augmented Lagrangian approach
 *
 * Copyright 2014-2025 by Knut Graichen, Andreas Voelz, Thore Wietzke,
 * Tobias Englert (<v2.3), Felix Mesmer (<v2.3), Soenke Rhein (<v2.3),
 * Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 * All rights reserved.
 *
 * GRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
 *
 */

#include "mex.h"
#include "grampc.h"
#include "grampc_conversion_Cmex.h"
#include "timing.h"

#ifndef N_TIMER
#define N_TIMER 1
#endif


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	typeGRAMPC* grampc;
	typeRNum* comptime;

	mxArray *mxuserparam = NULL;

	/* variables for runtime measurement */
    typeTime tic, toc;
#if N_TIMER > 1
    typeGRAMPC* grampcArray[N_TIMER-1];
#endif

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
	if (nlhs > 7) {
		mexErrMsgTxt("Too many output arguments.");
	}

	/* Allocate grampc struct and copy mxData */
	mxuserparam = mxDuplicateArray(mxGetField(prhs[0], 0, "userparam"));
	mx2typeGRAMPC(&grampc, prhs[0], mxuserparam);

	/* Allocate memory for the outputs */
	plhs[1] = mxCreateNumericMatrix(1, 1, mxtypeRNum_CLASS, mxREAL);
	comptime = (typeRNum*)mxGetPr(plhs[1]);
	*comptime = -1;

	/* runtime analysis */

#if N_TIMER > 1
	/* copy GRAMPC struct for runtime analysis */
	for (int i = 0; i < N_TIMER-1; i++) {
		mx2typeGRAMPC((grampcArray + i), prhs[0], mxuserparam);
	}
#endif

    timer_now(&tic);

	/* run GRAMPC while measuring time */
	grampc_run(grampc);

#if N_TIMER > 1
	for (int i = 0; i < N_TIMER-1; i++) {
		grampc_run(grampcArray[i]);
	}
#endif

    timer_now(&toc);
    *comptime = timer_diff_ms(&tic, &toc) / N_TIMER;

#if N_TIMER > 1
	/* free copies of GRAMPC struct */
	for (int i = 0; i < N_TIMER-1; i++) {
		grampc_free((grampcArray + i));
	}
#endif

	/* create mx structure */
	typeGRAMPC2mx(plhs, grampc, mxuserparam);

	/* free allocated memory */
	grampc_free(&grampc);
}
