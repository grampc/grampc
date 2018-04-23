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
