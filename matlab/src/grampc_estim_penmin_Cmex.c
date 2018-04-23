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
	if (nlhs > 2) {
		mexErrMsgTxt("Too many output arguments.");
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
