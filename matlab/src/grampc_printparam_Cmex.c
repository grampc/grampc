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
#include "grampc_setparam.h"
#include "grampc_conversion_Cmex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mxArray *mxuserparam = NULL;
	typeGRAMPC *grampc;

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
	if (nlhs > 0) {
		mexErrMsgTxt("Too many output arguments.");
	}

	/* remove warning unreferenced formal parameter */
	(void)(plhs);

	mxuserparam = mxGetField(prhs[0], 0, "userparam");
	mx2typeGRAMPC(&grampc, prhs[0], mxuserparam);

	grampc_printparam(grampc);
	grampc_free(&grampc);
}
