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
#include "grampc_setopt.h"
#include "grampc_mess.h"
#include "grampc_conversion_Cmex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mxArray *mxuserparam = NULL;
	typeGRAMPC *grampc = NULL;
	typeChar *optname, *optvalue;
	typeInt *intoptval;
	typeRNum *realoptval;

	/* check proper number of input arguments */
	if (nrhs < 3) {
		mexErrMsgTxt("Not enough input arguments.");
	}
	if (nrhs > 3) {
		mexErrMsgTxt("Too many input arguments.");
	}
	/* first argument: grampc structure  */
	if (!mxIsStruct(prhs[0])) {
		mexErrMsgTxt("First argument (grampc structure) must be a struct.");
	}
	/* second argument: MPC structure */
	if (!mxIsChar(prhs[1])) {
		mexErrMsgTxt("Second argument (option name) must be a string.");
	}
	/* third argument: option structure */
	if (!mxIsNumeric(prhs[2]) && !mxIsChar(prhs[2])) {
		mexErrMsgTxt("Third argument (option value) must be a string or a number.");
	}
	/* check proper number of output arguments */
	if (nlhs > 1) {
		mexErrMsgTxt("Too much output arguments. Only out argument is grampc structure.");
	}

	/* Allocate grampc struct and copy mxData */
	mxuserparam = mxDuplicateArray(mxGetField(prhs[0], 0, "userparam"));
	optname = mxArrayToString(prhs[1]);
	mx2typeGRAMPC(&grampc, prhs[0], mxuserparam);


	/* Set option using c-function */
	if (mxIsChar(prhs[2])) {
		optvalue = mxArrayToString(prhs[2]);
		grampc_setopt_string(grampc, optname, optvalue);
		mxFree(optvalue);
	}
	else if (!strcmp(optname, "MaxGradIter") || !strcmp(optname, "MaxMultIter") || !strcmp(optname, "Nhor") || !strcmp(optname, "IntegratorMaxSteps")) {
		grampc_setopt_int(grampc, optname, (typeInt)mxGetScalar(prhs[2]));
	}
	else if (!strcmp(optname, "FlagsRodas")) {
		if (mxGetClassID(prhs[2]) == mxtypeInt_CLASS) {
			grampc_setopt_int_vector(grampc, optname, (typeInt*)mxGetData(prhs[2]));
		}
		else if (mxGetClassID(prhs[2]) == mxDOUBLE_CLASS) {
			/* create temp array to cast the input values */
			CastDvec2Intvec(&intoptval, mxGetPr(prhs[2]), mxGetM(prhs[2])*mxGetN(prhs[2]));
			grampc_setopt_int_vector(grampc, optname, intoptval);
			free(intoptval);
		}
		else {
			grampc_error_addstring("Invalid datatyp of option ", optname);
		}
	}
	else if (!strcmp(optname, "xScale") || !strcmp(optname, "xOffset")
		|| !strcmp(optname, "uScale") || !strcmp(optname, "uOffset")
		|| !strcmp(optname, "pScale") || !strcmp(optname, "pOffset")
		|| !strcmp(optname, "cScale") || !strcmp(optname, "ConstraintsAbsTol")) {
		if (mxGetClassID(prhs[2]) == mxtypeRNum_CLASS) {
			grampc_setopt_real_vector(grampc, optname, (typeRNum*)mxGetPr(prhs[2]));
		}
		else if (mxGetClassID(prhs[2]) == mxDOUBLE_CLASS) {
			/* create temp array to cast the input values */
			CastDvec2Numvec(&realoptval, mxGetPr(prhs[2]), mxGetM(prhs[2])*mxGetN(prhs[2]));
			grampc_setopt_real_vector(grampc, optname, realoptval);
			free(realoptval);
		}
		else {
			grampc_error_addstring(INVALID_OPTION_DATATYP, optname);
		}
	}
	else {
		grampc_setopt_real(grampc, optname, (typeRNum)mxGetScalar(prhs[2]));
	}

	/* create mx structure */
	typeGRAMPC2mx(plhs, grampc, mxuserparam);

	/* free allocated memory */
	grampc_free(&grampc);
	mxFree(optname);
}
