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
#include "probfct.h"
#include "grampc_init.h"
#include "grampc_conversion_Cmex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	typeRNum *t, *x, *u, *p;
	typeRNum *icost;
	double *dicost;
	typeUSERPARAM *userparam;
	typeBoolean cast, cpy_t, cpy_x, cpy_u, cpy_p;
	typeInt Nx, Nu, Np, Ng, Nh, NgT, NhT;
    typeGRAMPCparam *param;

	/* Userparam is required for ocp_dim */
	if (nrhs != 6)
		mexErrMsgTxt("Wrong number of input arguments. Input arguments are t,x,u,p,grampc.param,userparam");
	userparam = (typeUSERPARAM*)mxGetData(prhs[5]);

	ocp_dim(&Nx, &Nu, &Np, &Ng, &Nh, &NgT, &NhT, userparam);

	/* check input arguments */
	if (mxGetM(prhs[0]) > 1 || mxGetN(prhs[0]) > 1)
		mexErrMsgTxt("Time (input argument 1) must be a scalar");
	if (mxGetM(prhs[1]) != (unsigned int)Nx && mxGetN(prhs[1]) != (unsigned int)Nx)
		mexErrMsgTxt("States (input argument 2) must have minimum length Nx");
	if (mxGetM(prhs[2]) != (unsigned int)Nu && mxGetN(prhs[2]) != (unsigned int)Nu)
		mexErrMsgTxt("Inputs (input argument 3) must have minimum length Nu");
	if (mxGetM(prhs[3]) != (unsigned int)Np && mxGetN(prhs[3]) != (unsigned int)Np)
		mexErrMsgTxt("Parameters (input argument 4) must have minimum length Np");
    if (!mxIsStruct(prhs[4]))
        mexErrMsgTxt("grampc.param (input argument 5) must be a struct.");

	/* check proper number of output arguments */
	if (nlhs > 1) {
		mexErrMsgTxt("Too many output arguments.");
	}

	/* Assign pointers to input arguments */
    mx2typeGRAMPCparam(&param, prhs[4]);
	cpy_t = AssignRealmxInput(&t, prhs, 0);
	cpy_x = AssignRealmxInput(&x, prhs, 1);
	cpy_u = AssignRealmxInput(&u, prhs, 2);
	cpy_p = AssignRealmxInput(&p, prhs, 3);
	cast = (cpy_t == 1) || (cpy_x == 1) || (cpy_u == 1) || (cpy_p == 1);

	/* Assign pointer to output argument */
	AssignRealmxOutput(&icost, &dicost, plhs, cast, 1, 0);

	/* function call */
	lfct(icost, t[0], x, u, p, param, userparam);

	/* Cast output */
	if (cast) {
		dicost[0] = (double)icost[0];
		free(icost);
	}

	/* If allocated memory for type conversion free it */
	if (cpy_t) free(t);
	if (cpy_x) free(x);
	if (cpy_u) free(u);
	if (cpy_p) free(p);

	/* free allocated memory */
    grampc_free_param(&param);
}
