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
#include "grampc_mess.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	typeInt status, level;
	typeChar *levelname;
	typeRNum *printed;

	/* check proper number of input arguments */
	if (nrhs < 1) {
		mexErrMsgTxt("Not enough input arguments.");
	}
	if (nrhs > 2) {
		mexErrMsgTxt("Too many input arguments.");
	}
	if (mxGetM(prhs[0]) > 1 || mxGetN(prhs[0]) > 1) {
		mexErrMsgTxt("Status (input argument 1) must be a scalar.");
	}

	/* check proper number of output arguments */
	if (nlhs > 1) {
		mexErrMsgTxt("Too many output arguments.");
	}
	if (nrhs < 2) {
		level = STATUS_LEVEL_ERROR;
	}
	else {
		levelname = mxArrayToString(prhs[1]);
		if (!strcmp(levelname, "Error") || !strcmp(levelname, "error") || !strcmp(levelname, "ERROR")) {
			level = STATUS_LEVEL_ERROR;
		}
		else if (!strcmp(levelname, "Warn") || !strcmp(levelname, "warn") || !strcmp(levelname, "WARN")
			|| !strcmp(levelname, "Warning") || !strcmp(levelname, "warning") || !strcmp(levelname, "WARNING")) {
			level = STATUS_LEVEL_WARN;
		}
		else if (!strcmp(levelname, "Info") || !strcmp(levelname, "info") || !strcmp(levelname, "INFO")) {
			level = STATUS_LEVEL_INFO;
		}
		else if (!strcmp(levelname, "Debug") || !strcmp(levelname, "debug") || !strcmp(levelname, "DEBUG")) {
			level = STATUS_LEVEL_DEBUG;
		}
		mxFree(levelname);
	}

	/* create a matrix for the output argument */
	plhs[0] = mxCreateNumericMatrix(1, 1, mxtypeRNum_CLASS, mxREAL);

	/* pointer to output argument */
	printed = (typeRNum*)mxGetPr(plhs[0]);

	/* assign the input argument */
	if (mxGetClassID(prhs[0]) == mxDOUBLE_CLASS) {
		status = (typeInt)(*((double*)mxGetData(prhs[0])));
	}
	else if (mxGetClassID(prhs[0]) == mxtypeInt_CLASS) {
		status = *((typeInt*)mxGetData(prhs[0]));
	}
	else {
		mexErrMsgTxt("Invalid datatype of first input argument.");
	}

	/* print status */
	*printed = (typeRNum)grampc_printstatus(status, level);
}
