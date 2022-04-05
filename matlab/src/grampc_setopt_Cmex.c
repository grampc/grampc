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
		mexErrMsgTxt("Too much output arguments. Only output argument is grampc structure.");
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
	else if (!strcmp(optname, "MaxGradIter") || !strcmp(optname, "MaxMultIter") || !strcmp(optname, "Nhor") || !strcmp(optname, "IntegratorMaxSteps"))
  {
    /* check dimension of the input arguments, all should be scalar values */
    if ((1 == mxGetN(prhs[2])) && (1 == mxGetM(prhs[2])))
    {
      grampc_setopt_int(grampc, optname, (typeInt)mxGetScalar(prhs[2]));
    }
    else
    {
      if (!strcmp(optname, "MaxGradIter"))
      {
        mexErrMsgTxt("Option MaxGradIter must be a scalar value.");
      }
      else if (!strcmp(optname, "MaxMultIter"))
      {
        mexErrMsgTxt("Option MaxMultIter must be a scalar value.");
      }
      else if (!strcmp(optname, "Nhor"))
      {
        mexErrMsgTxt("Option Nhor must be a scalar value.");
      }
      else if (!strcmp(optname, "IntegratorMaxSteps"))
      {
        mexErrMsgTxt("Option IntegratorMaxSteps must be a scalar value.");
      }
    }
	}
  else if (!strcmp(optname, "FlagsRodas"))
  {
    /* check dimension */
    if (((8 == mxGetN(prhs[2])) && (1 == mxGetM(prhs[2]))) || ((1 == mxGetN(prhs[2])) && (8 == mxGetM(prhs[2]))))
    {
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
    else
    {
      mexErrMsgTxt("Option FlagsRodas must be a vector of dimension 8.");
    }
	}
	else if (!strcmp(optname, "xScale") || !strcmp(optname, "xOffset")
		|| !strcmp(optname, "uScale") || !strcmp(optname, "uOffset")
		|| !strcmp(optname, "pScale") || !strcmp(optname, "pOffset")
		|| !strcmp(optname, "cScale") || !strcmp(optname, "ConstraintsAbsTol"))
  {
    /* check dimension of input arguments */
    const size_t loc_dimN = mxGetN(prhs[2]);
    const size_t loc_dimM = mxGetM(prhs[2]);
    if (!strcmp(optname, "xScale"))
    {
      const int loc_dimVec = grampc->param->Nx;
      if (!((loc_dimN == loc_dimVec) && (loc_dimM == 1)) && !((loc_dimN == 1) && (loc_dimM == loc_dimVec)))
      {
        mexErrMsgTxt("Wrong dimension for option xScale");
      }
    }
    else if (!strcmp(optname, "xOffset"))
    {
      const int loc_dimVec = grampc->param->Nx;
      if (!((loc_dimN == loc_dimVec) && (loc_dimM == 1)) && !((loc_dimN == 1) && (loc_dimM == loc_dimVec)))
      {
        mexErrMsgTxt("Wrong dimension for option xOffset");
      }
    }
    else if (!strcmp(optname, "uScale"))
    {
      const int loc_dimVec = grampc->param->Nu;
      if (!((loc_dimN == loc_dimVec) && (loc_dimM == 1)) && !((loc_dimN == 1) && (loc_dimM == loc_dimVec)))
      {
        mexErrMsgTxt("Wrong dimension for option uScale");
      }
    }
    else if (!strcmp(optname, "uOffset"))
    {
      const int loc_dimVec = grampc->param->Nu;
      if (!((loc_dimN == loc_dimVec) && (loc_dimM == 1)) && !((loc_dimN == 1) && (loc_dimM == loc_dimVec)))
      {
        mexErrMsgTxt("Wrong dimension for option uOffset");
      }
    }
    else if (!strcmp(optname, "pScale"))
    {
      const int loc_dimVec = grampc->param->Np;
      if (!((loc_dimN == loc_dimVec) && (loc_dimM == 1)) && !((loc_dimN == 1) && (loc_dimM == loc_dimVec)))
      {
        mexErrMsgTxt("Wrong dimension for option pScale");
      }
    }
    else if (!strcmp(optname, "pOffset"))
    {
      const int loc_dimVec = grampc->param->Np;
      if (!((loc_dimN == loc_dimVec) && (loc_dimM == 1)) && !((loc_dimN == 1) && (loc_dimM == loc_dimVec)))
      {
        mexErrMsgTxt("Wrong dimension for option pOffset");
      }
    }
    else if (!strcmp(optname, "cScale"))
    {
      const int loc_dimVec = grampc->param->Nc;
      if (!((loc_dimN == loc_dimVec) && (loc_dimM == 1)) && !((loc_dimN == 1) && (loc_dimM == loc_dimVec)))
      {
        mexErrMsgTxt("Wrong dimension for option cScale");
      }
    }
    else if (!strcmp(optname, "ConstraintsAbsTol"))
    {
      const int loc_dimVec = grampc->param->Nc;
      if (!((loc_dimN == loc_dimVec) && (loc_dimM == 1)) && !((loc_dimN == 1) && (loc_dimM == loc_dimVec)))
      {
        mexErrMsgTxt("Wrong dimension for option ConstraintsAbsTol");
      }
    }
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
	else
  {
    if ((1 == mxGetN(prhs[2])) && (1 == mxGetM(prhs[2])))
    {
      grampc_setopt_real(grampc, optname, (typeRNum)mxGetScalar(prhs[2]));
    }
    else
    {
      if (!strcmp(optname, "IntegratorRelTol"))
      {
        mexErrMsgTxt("Option IntegratorRelTol must be a scalar value.");
      }
      else if (!strcmp(optname, "IntegratorAbsTol"))
      {
        mexErrMsgTxt("Option IntegratorAbsTol must be a scalar value.");
      }
      else if (!strcmp(optname, "IntegratorMinStepSize"))
      {
        mexErrMsgTxt("Option IntegratorMinStepSize must be a scalar value.");
      }
      else if (!strcmp(optname, "LineSearchMax"))
      {
        mexErrMsgTxt("Option LineSearchMax must be a scalar value.");
      }
      else if (!strcmp(optname, "LineSearchMin"))
      {
        mexErrMsgTxt("Option LineSearchMin must be a scalar value.");
      }
      else if (!strcmp(optname, "LineSearchInit"))
      {
        mexErrMsgTxt("Option LineSearchInit must be a scalar value.");
      }
			else if (!strcmp(optname, "LineSearchAdaptAbsTol"))
			{
				mexErrMsgTxt("Option LineSearchAdaptAbsTol must be a scalar value.");
			}
      else if (!strcmp(optname, "LineSearchAdaptFactor"))
      {
        mexErrMsgTxt("Option LineSearchAdaptFactor must be a scalar value.");
      }
      else if (!strcmp(optname, "LineSearchIntervalTol"))
      {
        mexErrMsgTxt("Option LineSearchIntervalTol must be a scalar value.");
      }
			else if (!strcmp(optname, "LineSearchIntervalFactor"))
			{
				mexErrMsgTxt("Option LineSearchIntervalFactor must be a scalar value.");
			}
      else if (!strcmp(optname, "OptimParamLineSearchFactor"))
      {
        mexErrMsgTxt("Option OptimParamLineSearchFactor must be a scalar value.");
      }
      else if (!strcmp(optname, "OptimTimeLineSearchFactor"))
      {
        mexErrMsgTxt("Option OptimTimeLineSearchFactor must be a scalar value.");
      }
      else if (!strcmp(optname, "TScale"))
      {
        mexErrMsgTxt("Option TScale must be a scalar value.");
      }
      else if (!strcmp(optname, "TOffset"))
      {
        mexErrMsgTxt("Option TOffset must be a scalar value.");
      }
      else if (!strcmp(optname, "JScale"))
      {
        mexErrMsgTxt("Option JScale must be a scalar value.");
      }
      else if (!strcmp(optname, "MultiplierMax"))
      {
        mexErrMsgTxt("Option MultiplierMax must be a scalar value.");
      }
      else if (!strcmp(optname, "MultiplierDampingFactor"))
      {
        mexErrMsgTxt("Option MultiplierDampingFactor must be a scalar value.");
      }
      else if (!strcmp(optname, "PenaltyMax"))
      {
        mexErrMsgTxt("Option PenaltyMax must be a scalar value.");
      }
      else if (!strcmp(optname, "PenaltyMin"))
      {
        mexErrMsgTxt("Option PenaltyMin must be a scalar value.");
      }
      else if (!strcmp(optname, "PenaltyIncreaseFactor"))
      {
        mexErrMsgTxt("Option PenaltyIncreaseFactor must be a scalar value.");
      }
      else if (!strcmp(optname, "PenaltyDecreaseFactor"))
      {
        mexErrMsgTxt("Option PenaltyDecreaseFactor must be a scalar value.");
      }
      else if (!strcmp(optname, "PenaltyIncreaseThreshold"))
      {
        mexErrMsgTxt("Option PenaltyIncreaseThreshold must be a scalar value.");
      }
      else if (!strcmp(optname, "AugLagUpdateGradientRelTol"))
      {
        mexErrMsgTxt("Option AugLagUpdateGradientRelTol must be a scalar value.");
      }
      else if (!strcmp(optname, "ConvergenceGradientRelTol"))
      {
        mexErrMsgTxt("Option ConvergenceGradientRelTol must be a scalar value.");
      }
      /* Undefined optName */
      else
      {
        grampc_error_addstring(INVALID_OPTION_NAME, optname);
      }
    }
	}

	/* create mx structure */
	typeGRAMPC2mx(plhs, grampc, mxuserparam);

	/* free allocated memory */
	grampc_free(&grampc);
	mxFree(optname);
}
