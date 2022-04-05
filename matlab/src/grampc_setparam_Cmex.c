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
		|| !strcmp(paramname, "dt") || !strcmp(paramname, "t0"))
  {
    /* check dimension */
    const size_t loc_dimN = mxGetN(prhs[2]);
    const size_t loc_dimM = mxGetM(prhs[2]);
    if (!strcmp(paramname, "Thor"))
    {
      if (!((loc_dimN == 1) && (loc_dimM == 1)))
      {
        mexErrMsgTxt("Wrong dimension for parameter Thor");
      }
    }
    else if (!strcmp(paramname, "Tmax"))
    {
      if (!((loc_dimN == 1) && (loc_dimM == 1)))
      {
        mexErrMsgTxt("Wrong dimension for parameter Tmax");
      }
    }
    else if (!strcmp(paramname, "Tmin"))
    {
      if (!((loc_dimN == 1) && (loc_dimM == 1)))
      {
        mexErrMsgTxt("Wrong dimension for parameter Tmin");
      }
    }
    else if (!strcmp(paramname, "dt"))
    {
      if (!((loc_dimN == 1) && (loc_dimM == 1)))
      {
        mexErrMsgTxt("Wrong dimension for parameter dt");
      }
    }
    else if (!strcmp(paramname, "t0"))
    {
      if (!((loc_dimN == 1) && (loc_dimM == 1)))
      {
        mexErrMsgTxt("Wrong dimension for parameter t0");
      }
    }
		grampc_setparam_real(grampc, paramname, (typeRNum)mxGetScalar(prhs[2]));
	}
	else
  {
    /* check dimension */
    const size_t loc_dimN = mxGetN(prhs[2]);
    const size_t loc_dimM = mxGetM(prhs[2]);
    if (!strcmp(paramname, "x0"))
    {
      const int loc_dimVec = grampc->param->Nx;
      if (!((loc_dimN == loc_dimVec) && (loc_dimM == 1)) && !((loc_dimN == 1) && (loc_dimM == loc_dimVec)))
      {
        mexErrMsgTxt("Wrong dimension for parameter x0");
      }
    }
    else if (!strcmp(paramname, "xdes"))
    {
      const int loc_dimVec = grampc->param->Nx;
      if (!((loc_dimN == loc_dimVec) && (loc_dimM == 1)) && !((loc_dimN == 1) && (loc_dimM == loc_dimVec)))
      {
        mexErrMsgTxt("Wrong dimension for parameter xdes");
      }
    }
    else if (!strcmp(paramname, "u0"))
    {
      const int loc_dimVec = grampc->param->Nu;
      if (!((loc_dimN == loc_dimVec) && (loc_dimM == 1)) && !((loc_dimN == 1) && (loc_dimM == loc_dimVec)))
      {
        mexErrMsgTxt("Wrong dimension for parameter u0");
      }
    }
    else if (!strcmp(paramname, "udes"))
    {
      const int loc_dimVec = grampc->param->Nu;
      if (!((loc_dimN == loc_dimVec) && (loc_dimM == 1)) && !((loc_dimN == 1) && (loc_dimM == loc_dimVec)))
      {
        mexErrMsgTxt("Wrong dimension for parameter udes");
      }
    }
    else if (!strcmp(paramname, "umax"))
    {
      const int loc_dimVec = grampc->param->Nu;
      if (!((loc_dimN == loc_dimVec) && (loc_dimM == 1)) && !((loc_dimN == 1) && (loc_dimM == loc_dimVec)))
      {
        mexErrMsgTxt("Wrong dimension for parameter umax");
      }
    }
    else if (!strcmp(paramname, "umin"))
    {
      const int loc_dimVec = grampc->param->Nu;
      if (!((loc_dimN == loc_dimVec) && (loc_dimM == 1)) && !((loc_dimN == 1) && (loc_dimM == loc_dimVec)))
      {
        mexErrMsgTxt("Wrong dimension for parameter umin");
      }
    }
    else if (!strcmp(paramname, "p0"))
    {
      const int loc_dimVec = grampc->param->Np;
      if (!((loc_dimN == loc_dimVec) && (loc_dimM == 1)) && !((loc_dimN == 1) && (loc_dimM == loc_dimVec)))
      {
        mexErrMsgTxt("Wrong dimension for parameter p0");
      }
    }
    else if (!strcmp(paramname, "pmax"))
    {
      const int loc_dimVec = grampc->param->Np;
      if (!((loc_dimN == loc_dimVec) && (loc_dimM == 1)) && !((loc_dimN == 1) && (loc_dimM == loc_dimVec)))
      {
        mexErrMsgTxt("Wrong dimension for parameter pmax");
      }
    }
    else if (!strcmp(paramname, "pmin"))
    {
      const int loc_dimVec = grampc->param->Np;
      if (!((loc_dimN == loc_dimVec) && (loc_dimM == 1)) && !((loc_dimN == 1) && (loc_dimM == loc_dimVec)))
      {
        mexErrMsgTxt("Wrong dimension for parameter pmin");
      }
    }
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
