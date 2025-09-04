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
    typeChar *config_file_path;
	typeGRAMPC *grampc;
	mxArray *mxuserparam = NULL;

	/* Userparam is required for ocp_dim */
	if (nrhs > 2) {
        mexErrMsgTxt("Wrong number of input arguments. Input arguments are grampc and path to config file.");
	}
	if (!mxIsStruct(prhs[0])) {
		mexErrMsgTxt("First argument (grampc structure) must be a struct.");
	}
    if (!mxIsChar(prhs[1])) {
        mexErrMsgTxt("Second argument (path to config file) must be a char array.");
    }
	/* check proper number of output arguments */
    if (nlhs != 1) {
        mexErrMsgTxt("Wrong number of output arguments. Output argument is grampc.");
	}

	/* Allocate grampc struct and copy mxData */
	mxuserparam = mxDuplicateArray(mxGetField(prhs[0], 0, "userparam"));
	mx2typeGRAMPC(&grampc, prhs[0], mxuserparam);
    config_file_path = mxArrayToString(prhs[1]);

	/* Allocate memory for output argument */
	plhs[1] = mxCreateNumericMatrix(1, 1, mxtypeRNum_CLASS, mxREAL);
	
	/* function call */
    grampc_get_config_from_file(grampc, config_file_path);

	/* create grampc mx structure */
	typeGRAMPC2mx(plhs, grampc, mxuserparam);

	/* free allocated memory */
	grampc_free(&grampc);
}
