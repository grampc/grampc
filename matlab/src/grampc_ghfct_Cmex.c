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
#include "probfct.h"
#include "grampc_init.h"
#include "grampc_conversion_Cmex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	typeRNum *t, *x, *u, *p;
	typeRNum *cfct;
	double *dcfct;
	typeUSERPARAM *userparam;
	typeBoolean cast, cpy_t, cpy_x, cpy_u, cpy_p;
	typeInt i, Nx, Nu, Np, Ng, Nh, NgT, NhT;

	/* Userparam is required for ocp_dim */
	if (nrhs != 5)
		mexErrMsgTxt("Wrong number of input arguments. Input arguments are t,x,u,p,userparam");
	userparam = (typeUSERPARAM*)mxGetData(prhs[4]);

	ocp_dim(&Nx, &Nu, &Np, &Ng, &Nh, &NgT, &NhT, userparam);

	/* check input arguments */
	if (mxGetM(prhs[0]) > 1 || mxGetN(prhs[0]) > 1)
		mexErrMsgTxt("Time (input argument 1) must be a scalar");
	if (mxGetM(prhs[1]) != (unsigned int)Nx && mxGetN(prhs[1]) != (unsigned int)Nx)
		mexErrMsgTxt("States (input argument 2) must have minimum length Nx");
	if (mxGetM(prhs[2]) != (unsigned int)Nu && mxGetN(prhs[2]) != (unsigned int)Nu)
		mexErrMsgTxt("Inputs (input argument 3) must have minimum length Nu");
	if (mxGetM(prhs[3]) != (unsigned int)Np && mxGetN(prhs[3]) != (unsigned int)Np)
		mexErrMsgTxt("Parameters (input argument 4) must have minimum length p");

	/* check proper number of output arguments */
	if (nlhs > 1) {
		mexErrMsgTxt("Too many output arguments.");
	}

	/* Assign pointers to input arguments */
	cpy_t = AssignRealmxInput(&t, prhs, 0);
	cpy_x = AssignRealmxInput(&x, prhs, 1);
	cpy_u = AssignRealmxInput(&u, prhs, 2);
	cpy_p = AssignRealmxInput(&p, prhs, 3);
	cast = (cpy_t == 1) || (cpy_x == 1) || (cpy_u == 1) || (cpy_p == 1);

	/* Assign pointer to output argument */
	AssignRealmxOutput(&cfct, &dcfct, plhs, cast, Ng + Nh, 0);

	/* function call */
	gfct(cfct, t[0], x, u, p, userparam);
	hfct(cfct + Ng, t[0], x, u, p, userparam);

	/* cast result */
	if (cast) {
		for (i = 0; i < Ng + Nh; i++) {
			dcfct[i] = (double)cfct[i];
		}
		free(cfct);
	}

	/* If allocated memory for type conversion free it */
	if (cpy_t) free(t);
	if (cpy_x) free(x);
	if (cpy_u) free(u);
	if (cpy_p) free(p);
}
