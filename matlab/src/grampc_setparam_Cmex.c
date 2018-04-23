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

#include "grampc_init.h"
#include "grampc_setparam.h"
#include "grampc_mess.h"
#include "grampc_conversion_Cmex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mxArray *mxuserparam = NULL;
	typeGRAMPC *grampc;
	typeChar *paramname;
	typeRNum *realparamval;

	/* check proper number of input arguments */
	if (nrhs < 3) {
		mexErrMsgTxt("Not enough input arguments.");
	}
	if (nrhs > 3) {
		mexErrMsgTxt("Too many input arguments.");
	}
	/* first argument: param structure */
	if (!mxIsStruct(prhs[0])) {
		mexErrMsgTxt("First argument (grampc structure) must be a struct.");
	}
	/* Second argument: option name */
	if (!mxIsChar(prhs[1])) {
		mexErrMsgTxt("Second argument (option name) must be a string.");
	}
	/* Third argument: option value */
	if (!mxIsNumeric(prhs[2])) {
		mexErrMsgTxt("Third argument (option value) must be a number.");
	}
	/* check proper number of output arguments */
	if (nlhs > 1) {
		mexErrMsgTxt("Too much output arguments. Only out argument is grampc structure.");
	}

	/* Allocate grampc struct and copy mxData */
	mxuserparam = mxDuplicateArray(mxGetField(prhs[0], 0, "userparam"));
	mx2typeGRAMPC(&grampc, prhs[0], mxuserparam);
	paramname = mxArrayToString(prhs[1]);

	/* Set parameters */
	if (!strcmp(paramname, "Thor") || !strcmp(paramname, "Tmax") || !strcmp(paramname, "Tmin")
		|| !strcmp(paramname, "dt") || !strcmp(paramname, "t0")) {
		grampc_setparam_real(grampc, paramname, (typeRNum)mxGetScalar(prhs[2]));
	}
	else {
		/* Check datatype and set vector */
		if (mxGetClassID(prhs[2]) == mxtypeRNum_CLASS) {
			grampc_setparam_real_vector(grampc, paramname, (typeRNum*)mxGetPr(prhs[2]));
		}
		else if (mxGetClassID(prhs[2]) == mxDOUBLE_CLASS) {
			/* create temp array to cast the input values */
			CastDvec2Numvec(&realparamval, mxGetPr(prhs[2]), mxGetM(prhs[2])*mxGetN(prhs[2]));
			grampc_setparam_real_vector(grampc, paramname, realparamval);
			free(realparamval);
		}
		else {
			grampc_error_addstring(INVALID_PARAM_DATATYP, paramname);
		}
	}

	/* create mx structure */
	typeGRAMPC2mx(plhs, grampc, mxuserparam);

	/* free allocated memory */
	grampc_free(&grampc);
	mxFree(paramname);
}
