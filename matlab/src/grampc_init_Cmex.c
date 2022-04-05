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


/* include headers and functions */
#include "mex.h"
#include "grampc.h"
#include "grampc_conversion_Cmex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	typeGRAMPC *grampc;
	typeUSERPARAM *userparam;
	mxArray* mxuserparam;

	/* check proper number of input arguments */
	if (nrhs != 1) {
		mexErrMsgTxt("Wrong number of input arguments. Only Input argument is userparam.");
	}
	/* check proper number of output arguments */
	if (nlhs > 1) {
		mexErrMsgTxt("Too many output arguments. Only Output argument is grampc structure");
	}

	/* Memory allocation and initialization */
	mxuserparam = mxDuplicateArray(prhs[0]);
	userparam = (typeUSERPARAM*)mxGetData(mxuserparam);

	/* init grampc */
	grampc_init(&grampc, userparam);

	/* create mx structure */
	typeGRAMPC2mx(plhs, grampc, mxuserparam);

	/* free allocated memory */
	grampc_free(&grampc);
}
