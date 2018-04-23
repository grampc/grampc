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
