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
#include "grampc_conversion_Cmex.h"
#include "grampc.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	typeGRAMPC *grampc;
	mxArray *mxuserparam = NULL;

	/* Userparam is required for ocp_dim */
	if (nrhs != 3) {
		mexErrMsgTxt("Wrong number of input arguments. Input argument is grampc, tolerance and step size");
	}
	if (!mxIsStruct(prhs[0])) {
		mexErrMsgTxt("First argument (grampc structure) must be a struct.");
	}
	if (!mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2])) {
		mexErrMsgTxt("Tolerance and step size must be a double.");
	}

	/* Allocate grampc struct and copy mxData */
	mxuserparam = mxDuplicateArray(mxGetField(prhs[0], 0, "userparam"));
	mx2typeGRAMPC(&grampc, prhs[0], mxuserparam);
	
	/* function call */
	grampc_check_gradients(grampc, (typeRNum)mxGetScalar(prhs[1]), (typeRNum)mxGetScalar(prhs[2]));

	/* create grampc mx structure */
	typeGRAMPC2mx(plhs, grampc, mxuserparam);

	/* free allocated memory */
	grampc_free(&grampc);
}
