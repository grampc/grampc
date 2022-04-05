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

#include "mex.h" 
#include "grampc_conversion_Cmex.h"
#include "grampc.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *run_grampc;
	typeRNum* status;
	typeGRAMPC *grampc;
	mxArray *mxuserparam = NULL;

	/* Userparam is required for ocp_dim */
	if (nrhs > 2) {
		mexErrMsgTxt("Wrong number of input arguments. Input arguments are grampc and a flag for evaluating grampc_run");
	}
	if (!mxIsStruct(prhs[0])) {
		mexErrMsgTxt("First argument (grampc structure) must be a struct.");
	}
	if (mxGetM(prhs[1]) > 1 || mxGetN(prhs[1]) > 1) {
		mexErrMsgTxt("Second argument (flag for evaluating grampc_run) must be a scalar");
	}
	/* check proper number of output arguments */
	if (nlhs != 2) {
		mexErrMsgTxt("Wrong number of output arguments. Output arguments are grampc and the estimated value of PenaltyMin");
	}

	/* Allocate grampc struct and copy mxData */
	mxuserparam = mxDuplicateArray(mxGetField(prhs[0], 0, "userparam"));
	mx2typeGRAMPC(&grampc, prhs[0], mxuserparam);
	run_grampc = mxGetPr(prhs[1]);

	/* Allocate memory for output argument */
	plhs[1] = mxCreateNumericMatrix(1, 1, mxtypeRNum_CLASS, mxREAL);

	/* pointer to output argument */
	status = (typeRNum*)mxGetPr(plhs[1]);
	
	/* function call */
	*status = (typeRNum)grampc_estim_penmin(grampc, (typeInt)(*run_grampc));

	/* create grampc mx structure */
	typeGRAMPC2mx(plhs, grampc, mxuserparam);

	/* free allocated memory */
	grampc_free(&grampc);
}
