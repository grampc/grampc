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

#ifndef GRAMPC_CONVERSION_CMEX_H_
#define GRAMPC_CONVERSION_CMEX_H_

#include "grampc_util.h"
#include "grampc_alloc.h"
#include "grampc_mess.h"

#ifndef S_FUNCTION_NAME /* function not compatible with simulink coder */

void cAddNumScalar(typeRNum *cs, const char *fieldname, const mxArray *source) {
	*cs = (typeRNum)mxGetScalar(mxGetField(source, 0, fieldname));
}

void cAddIntScalar(typeInt *cs, const char *fieldname, const mxArray *source) {
	*cs = (typeInt)mxGetScalar(mxGetField(source, 0, fieldname));
}

typeInt cAddNumMatrix(typeRNum **cs, const char *fieldname, const mxArray *source) {
	mxArray *field = mxGetField(source, 0, fieldname);
	if (field != NULL) {
		size_t size = mxGetM(field)*mxGetN(field);
		if (size == 0) {
			*cs = NULL;
			return 1;
		}
		*cs = (typeRNum *)malloc(size * sizeof(typeRNum));
		if (*cs != NULL) {
			memcpy(*cs, (typeRNum*)mxGetPr(field), size * sizeof(typeRNum));
			return 1;
		}
		return -1;
	}
	else {
		grampc_error_addstring("struct from matlab does not have this field \n", fieldname);
		return -1;
	}
}

typeInt cAddIntMatrix(typeInt **cs, const char *fieldname, const mxArray *source) {
	mxArray *field = mxGetField(source, 0, fieldname);
	if (field != NULL) {
		size_t size = mxGetM(field)*mxGetN(field);
		if (size == 0) {
			*cs = NULL;
			return 1;
		}
		*cs = (typeInt *)malloc(size * sizeof(typeInt));
		if (*cs != NULL) {
			memcpy(*cs, (typeInt*)mxGetPr(field), size * sizeof(typeInt));
			return 1;
		}
		return -1;
	}
	else {
		grampc_error_addstring("struct from matlab does not have this field \n", fieldname);
		return -1;
	}
}


typeInt mxAddNumMatrix(mxArray *pm, const char *fieldname, ctypeRNum *value, ctypeInt rows, ctypeInt cols) {
	typeInt idx = mxAddField(pm, fieldname);
	if (idx != -1) {
		mxArray* field = mxCreateNumericMatrix(rows, cols, mxtypeRNum_CLASS, mxREAL);
		memcpy(mxGetPr(field), value, rows * cols * sizeof(typeRNum));
		mxSetField(pm, 0, fieldname, field);
	}
	return idx;
}

typeInt mxAddIntMatrix(mxArray *pm, const char *fieldname, ctypeInt*value, ctypeInt rows, ctypeInt cols) {
	typeInt idx = mxAddField(pm, fieldname);
	if (idx != -1) {
		mxArray* field = mxCreateNumericMatrix(rows, cols, mxtypeInt_CLASS, mxREAL);
		memcpy(mxGetPr(field), value, rows * cols * sizeof(typeInt));
		mxSetField(pm, 0, fieldname, field);
	}
	return idx;
}

typeInt mxAddStructField(mxArray *pm, const char *fieldname, mxArray*field) {
	typeInt idx = mxAddField(pm, fieldname);
	if (idx != -1) {
		mxSetField(pm, 0, fieldname, field);
	}
	return idx;
}


void mx2typeGRAMPC(typeGRAMPC**grampc, const mxArray *mxgrampc, const mxArray *mxuserparam) {

	typeGRAMPCparam* param;
	typeGRAMPCopt* opt;
	typeGRAMPCsol* sol;
	typeGRAMPCrws* rws;

	mxArray *mxparam = NULL;
	mxArray *mxopt = NULL;
	mxArray *mxsol = NULL;
	mxArray *mxrws = NULL;

	mxparam = mxGetField(mxgrampc, 0, "param");
	mxopt = mxGetField(mxgrampc, 0, "opt");
	mxsol = mxGetField(mxgrampc, 0, "sol");
	mxrws = mxGetField(mxgrampc, 0, "rws");

	/* STRUCTURE MEMORY ALLOCATION **********************************************/
	*grampc = (typeGRAMPC *)calloc(1, sizeof(**grampc));
	if (*grampc == NULL) {
		grampc_error(GRAMPC_ALLOC_FAILED);
	}
	(*grampc)->param = (typeGRAMPCparam *)calloc(1, sizeof(*(*grampc)->param));
	if ((*grampc)->param == NULL) {
		grampc_error(PARAM_ALLOC_FAILED);
	}
	(*grampc)->sol = (typeGRAMPCsol *)calloc(1, sizeof(*(*grampc)->sol));
	if ((*grampc)->sol == NULL) {
		grampc_error(SOL_ALLOC_FAILED);
	}
	(*grampc)->rws = (typeGRAMPCrws *)calloc(1, sizeof(*(*grampc)->rws));
	if ((*grampc)->rws == NULL) {
		grampc_error(RWS_ALLOC_FAILED);
	}
	(*grampc)->opt = (typeGRAMPCopt *)calloc(1, sizeof(*(*grampc)->opt));
	if ((*grampc)->opt == NULL) {
		grampc_error(OPT_ALLOC_FAILED);
	}

	param = ((*grampc)->param);
	opt = ((*grampc)->opt);
	sol = ((*grampc)->sol);
	rws = ((*grampc)->rws);

	/* PARAM STRUCTURE **********************************************************/
	cAddIntScalar(&param->Nx, "Nx", mxparam);
	cAddIntScalar(&param->Nu, "Nu", mxparam);
	cAddIntScalar(&param->Np, "Np", mxparam);
	cAddIntScalar(&param->Ng, "Ng", mxparam);
	cAddIntScalar(&param->Nh, "Nh", mxparam);
	cAddIntScalar(&param->NgT, "NgT", mxparam);
	cAddIntScalar(&param->NhT, "NhT", mxparam);
	cAddIntScalar(&param->Nc, "Nc", mxparam);

	cAddNumMatrix(&param->x0, "x0", mxparam);
	cAddNumMatrix(&param->xdes, "xdes", mxparam);

	cAddNumMatrix(&param->u0, "u0", mxparam);
	cAddNumMatrix(&param->udes, "udes", mxparam);
	cAddNumMatrix(&param->umax, "umax", mxparam);
	cAddNumMatrix(&param->umin, "umin", mxparam);

	cAddNumMatrix(&param->p0, "p0", mxparam);
	cAddNumMatrix(&param->pmax, "pmax", mxparam);
	cAddNumMatrix(&param->pmin, "pmin", mxparam);

	cAddNumScalar(&param->Thor, "Thor", mxparam);
	cAddNumScalar(&param->Tmax, "Tmax", mxparam);
	cAddNumScalar(&param->Tmin, "Tmin", mxparam);

	cAddNumScalar(&param->dt, "dt", mxparam);
	cAddNumScalar(&param->t0, "t0", mxparam);

	/* OPTIONS STRUCTURE ******************************************************/
	cAddIntScalar(&opt->Nhor, "Nhor", mxopt);
	cAddIntScalar(&opt->MaxGradIter, "MaxGradIter", mxopt);
	cAddIntScalar(&opt->MaxMultIter, "MaxMultIter", mxopt);
	cAddIntScalar(&opt->ShiftControl, "ShiftControl", mxopt);

	cAddIntScalar(&opt->TimeDiscretization, "TimeDiscretization", mxopt);

	cAddIntScalar(&opt->IntegralCost, "IntegralCost", mxopt);
	cAddIntScalar(&opt->TerminalCost, "TerminalCost", mxopt);
	cAddIntScalar(&opt->IntegratorCost, "IntegratorCost", mxopt);

	cAddIntScalar(&opt->Integrator, "Integrator", mxopt);
	cAddNumScalar(&opt->IntegratorRelTol, "IntegratorRelTol", mxopt);
	cAddNumScalar(&opt->IntegratorAbsTol, "IntegratorAbsTol", mxopt);
	cAddNumScalar(&opt->IntegratorMinStepSize, "IntegratorMinStepSize", mxopt);
	cAddIntScalar(&opt->IntegratorMaxSteps, "IntegratorMaxSteps", mxopt);
	/*cAddIntMatrix(&opt->FlagsRodas, "FlagsRodas", mxopt);*/

	cAddIntScalar(&opt->LineSearchType, "LineSearchType", mxopt);
	cAddIntScalar(&opt->LineSearchExpAutoFallback, "LineSearchExpAutoFallback", mxopt);
	cAddNumScalar(&opt->LineSearchMax, "LineSearchMax", mxopt);
	cAddNumScalar(&opt->LineSearchMin, "LineSearchMin", mxopt);
	cAddNumScalar(&opt->LineSearchInit, "LineSearchInit", mxopt);
	cAddNumScalar(&opt->LineSearchAdaptAbsTol, "LineSearchAdaptAbsTol", mxopt);
	cAddNumScalar(&opt->LineSearchAdaptFactor, "LineSearchAdaptFactor", mxopt);
	cAddNumScalar(&opt->LineSearchIntervalTol, "LineSearchIntervalTol", mxopt);
	cAddNumScalar(&opt->LineSearchIntervalFactor, "LineSearchIntervalFactor", mxopt);

	cAddIntScalar(&opt->OptimControl, "OptimControl", mxopt);
	cAddIntScalar(&opt->OptimParam, "OptimParam", mxopt);
	cAddNumScalar(&opt->OptimParamLineSearchFactor, "OptimParamLineSearchFactor", mxopt);
	cAddIntScalar(&opt->OptimTime, "OptimTime", mxopt);
	cAddNumScalar(&opt->OptimTimeLineSearchFactor, "OptimTimeLineSearchFactor", mxopt);

	cAddIntScalar(&opt->ScaleProblem, "ScaleProblem", mxopt);
	cAddNumMatrix(&opt->xScale, "xScale", mxopt);
	cAddNumMatrix(&opt->xOffset, "xOffset", mxopt);
	cAddNumMatrix(&opt->uScale, "uScale", mxopt);
	cAddNumMatrix(&opt->uOffset, "uOffset", mxopt);
	cAddNumMatrix(&opt->pScale, "pScale", mxopt);
	cAddNumMatrix(&opt->pOffset, "pOffset", mxopt);
	cAddNumScalar(&opt->TScale, "TScale", mxopt);
	cAddNumScalar(&opt->TOffset, "TOffset", mxopt);
	cAddNumScalar(&opt->JScale, "JScale", mxopt);
	cAddNumMatrix(&opt->cScale, "cScale", mxopt);

	cAddIntScalar(&opt->EqualityConstraints, "EqualityConstraints", mxopt);
	cAddIntScalar(&opt->InequalityConstraints, "InequalityConstraints", mxopt);
	cAddIntScalar(&opt->TerminalEqualityConstraints, "TerminalEqualityConstraints", mxopt);
	cAddIntScalar(&opt->TerminalInequalityConstraints, "TerminalInequalityConstraints", mxopt);
	cAddIntScalar(&opt->ConstraintsHandling, "ConstraintsHandling", mxopt);
	cAddNumMatrix(&opt->ConstraintsAbsTol, "ConstraintsAbsTol", mxopt);

	cAddNumScalar(&opt->MultiplierMax, "MultiplierMax", mxopt);
	cAddNumScalar(&opt->MultiplierDampingFactor, "MultiplierDampingFactor", mxopt);
	cAddNumScalar(&opt->PenaltyMax, "PenaltyMax", mxopt);
	cAddNumScalar(&opt->PenaltyMin, "PenaltyMin", mxopt);
	cAddNumScalar(&opt->PenaltyIncreaseFactor, "PenaltyIncreaseFactor", mxopt);
	cAddNumScalar(&opt->PenaltyDecreaseFactor, "PenaltyDecreaseFactor", mxopt);
	cAddNumScalar(&opt->PenaltyIncreaseThreshold, "PenaltyIncreaseThreshold", mxopt);
	cAddNumScalar(&opt->AugLagUpdateGradientRelTol, "AugLagUpdateGradientRelTol", mxopt);

	cAddIntScalar(&opt->ConvergenceCheck, "ConvergenceCheck", mxopt);
	cAddNumScalar(&opt->ConvergenceGradientRelTol, "ConvergenceGradientRelTol", mxopt);

	/* SOLUTION STRUCTURE ******************************************************/
	cAddNumMatrix(&sol->xnext, "xnext", mxsol);
	cAddNumMatrix(&sol->unext, "unext", mxsol);
	cAddNumMatrix(&sol->pnext, "pnext", mxsol);
	cAddNumScalar(&sol->Tnext, "Tnext", mxsol);
	/*cAddNumMatrix(&sol->J, "J", mxsol);*/
	cAddNumScalar(&sol->cfct, "cfct", mxsol);
	cAddNumScalar(&sol->pen, "pen", mxsol);
	cAddIntScalar(&sol->status, "status", mxsol);
	cAddIntMatrix(&sol->iter, "iter", mxsol);

	/* RWS STRUCTURE **********************************************************/
	cAddNumMatrix(&rws->t, "t", mxrws);
	cAddNumMatrix(&rws->tls, "tls", mxrws);

	cAddNumMatrix(&rws->x, "x", mxrws);
	cAddNumMatrix(&rws->adj, "adj", mxrws);
	cAddNumMatrix(&rws->dcdx, "dcdx", mxrws);

	cAddNumMatrix(&rws->u, "u", mxrws);
	cAddNumMatrix(&rws->uls, "uls", mxrws);
	cAddNumMatrix(&rws->uprev, "uprev", mxrws);
	cAddNumMatrix(&rws->gradu, "gradu", mxrws);
	cAddNumMatrix(&rws->graduprev, "graduprev", mxrws);
	cAddNumMatrix(&rws->dcdu, "dcdu", mxrws);

	cAddNumMatrix(&rws->p, "p", mxrws);
	cAddNumMatrix(&rws->pls, "pls", mxrws);
	cAddNumMatrix(&rws->pprev, "pprev", mxrws);
	cAddNumMatrix(&rws->gradp, "gradp", mxrws);
	cAddNumMatrix(&rws->gradpprev, "gradpprev", mxrws);
	cAddNumMatrix(&rws->dcdp, "dcdp", mxrws);

	cAddNumScalar(&rws->T, "T", mxrws);
	cAddNumScalar(&rws->Tprev, "Tprev", mxrws);
	cAddNumScalar(&rws->gradT, "gradT", mxrws);
	cAddNumScalar(&rws->gradTprev, "gradTprev", mxrws);
	cAddNumScalar(&rws->dcdt, "dcdt", mxrws);

	cAddNumMatrix(&rws->mult, "mult", mxrws);
	cAddNumMatrix(&rws->pen, "pen", mxrws);
	cAddNumMatrix(&rws->cfct, "cfct", mxrws);
	cAddNumMatrix(&rws->cfctprev, "cfctprev", mxrws);
	cAddNumMatrix(&rws->cfctAbsTol, "cfctAbsTol", mxrws);

	cAddNumMatrix(&rws->lsAdapt, "lsAdapt", mxrws);
	cAddNumMatrix(&rws->lsExplicit, "lsExplicit", mxrws);
	cAddNumMatrix(&rws->rwsScale, "rwsScale", mxrws);
	cAddIntScalar(&rws->lrwsGeneral, "lrwsGeneral", mxrws);
	cAddNumMatrix(&rws->rwsGeneral, "rwsGeneral", mxrws);

	cAddIntScalar(&rws->lworkRodas, "lworkRodas", mxrws);
	cAddIntScalar(&rws->liworkRodas, "liworkRodas", mxrws);
	cAddNumMatrix(&rws->rparRodas, "rparRodas", mxrws);
	cAddIntMatrix(&rws->iparRodas, "iparRodas", mxrws);
	cAddNumMatrix(&rws->workRodas, "workRodas", mxrws);
	cAddIntMatrix(&rws->iworkRodas, "iworkRodas", mxrws);

	/* USERPARAM **************************************************************/
	(*grampc)->userparam = (typeUSERPARAM*)mxGetData(mxuserparam);
}


void typeGRAMPC2mx(mxArray *plhs[], const typeGRAMPC*grampc, mxArray *mxuserparam) {

	/* generation of outputs */
	mxArray* mxgrampcstruct = mxCreateStructMatrix(1, 1, 0, NULL);
	mxArray* mxparamstruct = mxCreateStructMatrix(1, 1, 0, NULL);
	mxArray* mxoptstruct = mxCreateStructMatrix(1, 1, 0, NULL);
	mxArray* mxsolstruct = mxCreateStructMatrix(1, 1, 0, NULL);
	mxArray* mxrwsstruct = mxCreateStructMatrix(1, 1, 0, NULL);

	/* PARAM STRUCTURE **********************************************************/
	mxAddIntMatrix(mxparamstruct, "Nx", &grampc->param->Nx, 1, 1);
	mxAddIntMatrix(mxparamstruct, "Nu", &grampc->param->Nu, 1, 1);
	mxAddIntMatrix(mxparamstruct, "Np", &grampc->param->Np, 1, 1);
	mxAddIntMatrix(mxparamstruct, "Ng", &grampc->param->Ng, 1, 1);
	mxAddIntMatrix(mxparamstruct, "Nh", &grampc->param->Nh, 1, 1);
	mxAddIntMatrix(mxparamstruct, "NgT", &grampc->param->NgT, 1, 1);
	mxAddIntMatrix(mxparamstruct, "NhT", &grampc->param->NhT, 1, 1);
	mxAddIntMatrix(mxparamstruct, "Nc", &grampc->param->Nc, 1, 1);

	mxAddNumMatrix(mxparamstruct, "x0", grampc->param->x0, grampc->param->Nx, 1);
	mxAddNumMatrix(mxparamstruct, "xdes", grampc->param->xdes, grampc->param->Nx, 1);

	mxAddNumMatrix(mxparamstruct, "u0", grampc->param->u0, grampc->param->Nu, 1);
	mxAddNumMatrix(mxparamstruct, "udes", grampc->param->udes, grampc->param->Nu, 1);
	mxAddNumMatrix(mxparamstruct, "umax", grampc->param->umax, grampc->param->Nu, 1);
	mxAddNumMatrix(mxparamstruct, "umin", grampc->param->umin, grampc->param->Nu, 1);

	mxAddNumMatrix(mxparamstruct, "p0", grampc->param->p0, grampc->param->Np, 1);
	mxAddNumMatrix(mxparamstruct, "pmax", grampc->param->pmax, grampc->param->Np, 1);
	mxAddNumMatrix(mxparamstruct, "pmin", grampc->param->pmin, grampc->param->Np, 1);

	mxAddNumMatrix(mxparamstruct, "Thor", &grampc->param->Thor, 1, 1);
	mxAddNumMatrix(mxparamstruct, "Tmax", &grampc->param->Tmax, 1, 1);
	mxAddNumMatrix(mxparamstruct, "Tmin", &grampc->param->Tmin, 1, 1);

	mxAddNumMatrix(mxparamstruct, "dt", &grampc->param->dt, 1, 1);
	mxAddNumMatrix(mxparamstruct, "t0", &grampc->param->t0, 1, 1);

	mxAddStructField(mxgrampcstruct, "param", mxparamstruct);

	/* OPTIONS STRUCTURE ******************************************************/
	mxAddIntMatrix(mxoptstruct, "Nhor", &grampc->opt->Nhor, 1, 1);
	mxAddIntMatrix(mxoptstruct, "MaxGradIter", &grampc->opt->MaxGradIter, 1, 1);
	mxAddIntMatrix(mxoptstruct, "MaxMultIter", &grampc->opt->MaxMultIter, 1, 1);
	mxAddIntMatrix(mxoptstruct, "ShiftControl", &grampc->opt->ShiftControl, 1, 1);

	mxAddIntMatrix(mxoptstruct, "TimeDiscretization", &grampc->opt->TimeDiscretization, 1, 1);

	mxAddIntMatrix(mxoptstruct, "IntegralCost", &grampc->opt->IntegralCost, 1, 1);
	mxAddIntMatrix(mxoptstruct, "TerminalCost", &grampc->opt->TerminalCost, 1, 1);
	mxAddIntMatrix(mxoptstruct, "IntegratorCost", &grampc->opt->IntegratorCost, 1, 1);

	mxAddIntMatrix(mxoptstruct, "Integrator", &grampc->opt->Integrator, 1, 1);
	mxAddNumMatrix(mxoptstruct, "IntegratorRelTol", &grampc->opt->IntegratorRelTol, 1, 1);
	mxAddNumMatrix(mxoptstruct, "IntegratorAbsTol", &grampc->opt->IntegratorAbsTol, 1, 1);
	mxAddNumMatrix(mxoptstruct, "IntegratorMinStepSize", &grampc->opt->IntegratorMinStepSize, 1, 1);
	mxAddIntMatrix(mxoptstruct, "IntegratorMaxSteps", &grampc->opt->IntegratorMaxSteps, 1, 1);
	mxAddIntMatrix(mxoptstruct, "FlagsRodas", grampc->opt->FlagsRodas, 1, 8);

	mxAddIntMatrix(mxoptstruct, "LineSearchType", &grampc->opt->LineSearchType, 1, 1);
	mxAddIntMatrix(mxoptstruct, "LineSearchExpAutoFallback", &grampc->opt->LineSearchExpAutoFallback, 1, 1);
	mxAddNumMatrix(mxoptstruct, "LineSearchMax", &grampc->opt->LineSearchMax, 1, 1);
	mxAddNumMatrix(mxoptstruct, "LineSearchMin", &grampc->opt->LineSearchMin, 1, 1);
	mxAddNumMatrix(mxoptstruct, "LineSearchInit", &grampc->opt->LineSearchInit, 1, 1);
	mxAddNumMatrix(mxoptstruct, "LineSearchAdaptAbsTol", &grampc->opt->LineSearchAdaptAbsTol, 1, 1);
	mxAddNumMatrix(mxoptstruct, "LineSearchAdaptFactor", &grampc->opt->LineSearchAdaptFactor, 1, 1);
	mxAddNumMatrix(mxoptstruct, "LineSearchIntervalTol", &grampc->opt->LineSearchIntervalTol, 1, 1);
	mxAddNumMatrix(mxoptstruct, "LineSearchIntervalFactor", &grampc->opt->LineSearchIntervalFactor, 1, 1);

	mxAddIntMatrix(mxoptstruct, "OptimControl", &grampc->opt->OptimControl, 1, 1);
	mxAddIntMatrix(mxoptstruct, "OptimParam", &grampc->opt->OptimParam, 1, 1);
	mxAddNumMatrix(mxoptstruct, "OptimParamLineSearchFactor", &grampc->opt->OptimParamLineSearchFactor, 1, 1);
	mxAddIntMatrix(mxoptstruct, "OptimTime", &grampc->opt->OptimTime, 1, 1);
	mxAddNumMatrix(mxoptstruct, "OptimTimeLineSearchFactor", &grampc->opt->OptimTimeLineSearchFactor, 1, 1);

	mxAddIntMatrix(mxoptstruct, "ScaleProblem", &grampc->opt->ScaleProblem, 1, 1);
	mxAddNumMatrix(mxoptstruct, "xScale", grampc->opt->xScale, grampc->param->Nx, 1);
	mxAddNumMatrix(mxoptstruct, "xOffset", grampc->opt->xOffset, grampc->param->Nx, 1);
	mxAddNumMatrix(mxoptstruct, "uScale", grampc->opt->uScale, grampc->param->Nu, 1);
	mxAddNumMatrix(mxoptstruct, "uOffset", grampc->opt->uOffset, grampc->param->Nu, 1);
	mxAddNumMatrix(mxoptstruct, "pScale", grampc->opt->pScale, grampc->param->Np, 1);
	mxAddNumMatrix(mxoptstruct, "pOffset", grampc->opt->pOffset, grampc->param->Np, 1);
	mxAddNumMatrix(mxoptstruct, "TScale", &grampc->opt->TScale, 1, 1);
	mxAddNumMatrix(mxoptstruct, "TOffset", &grampc->opt->TOffset, 1, 1);
	mxAddNumMatrix(mxoptstruct, "JScale", &grampc->opt->JScale, 1, 1);
	mxAddNumMatrix(mxoptstruct, "cScale", grampc->opt->cScale, grampc->param->Nc, 1);

	mxAddIntMatrix(mxoptstruct, "EqualityConstraints", &grampc->opt->EqualityConstraints, 1, 1);
	mxAddIntMatrix(mxoptstruct, "InequalityConstraints", &grampc->opt->InequalityConstraints, 1, 1);
	mxAddIntMatrix(mxoptstruct, "TerminalEqualityConstraints", &grampc->opt->TerminalEqualityConstraints, 1, 1);
	mxAddIntMatrix(mxoptstruct, "TerminalInequalityConstraints", &grampc->opt->TerminalInequalityConstraints, 1, 1);
	mxAddIntMatrix(mxoptstruct, "ConstraintsHandling", &grampc->opt->ConstraintsHandling, 1, 1);
	mxAddNumMatrix(mxoptstruct, "ConstraintsAbsTol", grampc->opt->ConstraintsAbsTol, 1, grampc->param->Nc);

	mxAddNumMatrix(mxoptstruct, "MultiplierMax", &grampc->opt->MultiplierMax, 1, 1);
	mxAddNumMatrix(mxoptstruct, "MultiplierDampingFactor", &grampc->opt->MultiplierDampingFactor, 1, 1);
	mxAddNumMatrix(mxoptstruct, "PenaltyMax", &grampc->opt->PenaltyMax, 1, 1);
	mxAddNumMatrix(mxoptstruct, "PenaltyMin", &grampc->opt->PenaltyMin, 1, 1);
	mxAddNumMatrix(mxoptstruct, "PenaltyIncreaseFactor", &grampc->opt->PenaltyIncreaseFactor, 1, 1);
	mxAddNumMatrix(mxoptstruct, "PenaltyDecreaseFactor", &grampc->opt->PenaltyDecreaseFactor, 1, 1);
	mxAddNumMatrix(mxoptstruct, "PenaltyIncreaseThreshold", &grampc->opt->PenaltyIncreaseThreshold, 1, 1);
	mxAddNumMatrix(mxoptstruct, "AugLagUpdateGradientRelTol", &grampc->opt->AugLagUpdateGradientRelTol, 1, 1);

	mxAddIntMatrix(mxoptstruct, "ConvergenceCheck", &grampc->opt->ConvergenceCheck, 1, 1);
	mxAddNumMatrix(mxoptstruct, "ConvergenceGradientRelTol", &grampc->opt->ConvergenceGradientRelTol, 1, 1);
	mxAddStructField(mxgrampcstruct, "opt", mxoptstruct);

	/* SOLUTION STRUCTURE ******************************************************/
	mxAddNumMatrix(mxsolstruct, "xnext", grampc->sol->xnext, grampc->param->Nx, 1);
	mxAddNumMatrix(mxsolstruct, "unext", grampc->sol->unext, grampc->param->Nu, 1);
	mxAddNumMatrix(mxsolstruct, "pnext", grampc->sol->pnext, grampc->param->Np, 1);
	mxAddNumMatrix(mxsolstruct, "Tnext", &grampc->sol->Tnext, 1, 1);
	mxAddNumMatrix(mxsolstruct, "J", grampc->sol->J, 1, 2);
	mxAddNumMatrix(mxsolstruct, "cfct", &grampc->sol->cfct, 1, 1);
	mxAddNumMatrix(mxsolstruct, "pen", &grampc->sol->pen, 1, 1);
	mxAddIntMatrix(mxsolstruct, "iter", grampc->sol->iter, 1, grampc->opt->MaxMultIter);
	mxAddIntMatrix(mxsolstruct, "status", &grampc->sol->status, 1, 1);

	mxAddStructField(mxgrampcstruct, "sol", mxsolstruct);

	/* RWS STRUCTURE **********************************************************/
	mxAddNumMatrix(mxrwsstruct, "t", grampc->rws->t, 1, grampc->opt->Nhor);
	mxAddNumMatrix(mxrwsstruct, "tls", grampc->rws->tls, 1, grampc->opt->Nhor);

	mxAddNumMatrix(mxrwsstruct, "x", grampc->rws->x, grampc->param->Nx, grampc->opt->Nhor);
	mxAddNumMatrix(mxrwsstruct, "adj", grampc->rws->adj, grampc->param->Nx, grampc->opt->Nhor);
	mxAddNumMatrix(mxrwsstruct, "dcdx", grampc->rws->dcdx, grampc->param->Nx, grampc->opt->Nhor + 1);

	mxAddNumMatrix(mxrwsstruct, "u", grampc->rws->u, grampc->param->Nu, grampc->opt->Nhor);
	mxAddNumMatrix(mxrwsstruct, "uls", grampc->rws->uls, grampc->param->Nu, grampc->opt->Nhor);
	mxAddNumMatrix(mxrwsstruct, "uprev", grampc->rws->uprev, grampc->param->Nu, grampc->opt->Nhor);
	mxAddNumMatrix(mxrwsstruct, "gradu", grampc->rws->gradu, grampc->param->Nu, grampc->opt->Nhor);
	mxAddNumMatrix(mxrwsstruct, "graduprev", grampc->rws->graduprev, grampc->param->Nu, grampc->opt->Nhor);
	mxAddNumMatrix(mxrwsstruct, "dcdu", grampc->rws->dcdu, grampc->param->Nu, grampc->opt->Nhor);

	mxAddNumMatrix(mxrwsstruct, "p", grampc->rws->p, grampc->param->Np, 1);
	mxAddNumMatrix(mxrwsstruct, "pls", grampc->rws->pls, grampc->param->Np, 1);
	mxAddNumMatrix(mxrwsstruct, "pprev", grampc->rws->pprev, grampc->param->Np, 1);
	mxAddNumMatrix(mxrwsstruct, "gradp", grampc->rws->gradp, grampc->param->Np, 1);
	mxAddNumMatrix(mxrwsstruct, "gradpprev", grampc->rws->gradpprev, grampc->param->Np, 1);
	mxAddNumMatrix(mxrwsstruct, "dcdp", grampc->rws->dcdp, grampc->param->Np, grampc->opt->Nhor + 1);

	mxAddNumMatrix(mxrwsstruct, "T", &grampc->rws->T, 1, 1);
	mxAddNumMatrix(mxrwsstruct, "Tprev", &grampc->rws->Tprev, 1, 1);
	mxAddNumMatrix(mxrwsstruct, "gradT", &grampc->rws->gradT, 1, 1);
	mxAddNumMatrix(mxrwsstruct, "gradTprev", &grampc->rws->gradTprev, 1, 1);
	mxAddNumMatrix(mxrwsstruct, "dcdt", &grampc->rws->dcdt, 1, 1);

	mxAddNumMatrix(mxrwsstruct, "mult", grampc->rws->mult, grampc->param->Nc, grampc->opt->Nhor);
	mxAddNumMatrix(mxrwsstruct, "pen", grampc->rws->pen, grampc->param->Nc, grampc->opt->Nhor);
	mxAddNumMatrix(mxrwsstruct, "cfct", grampc->rws->cfct, grampc->param->Nc, grampc->opt->Nhor);
	mxAddNumMatrix(mxrwsstruct, "cfctprev", grampc->rws->cfctprev, grampc->param->Nc, grampc->opt->Nhor);
	mxAddNumMatrix(mxrwsstruct, "cfctAbsTol", grampc->rws->cfctAbsTol, 1, grampc->param->Nc);

	if (grampc->opt->LineSearchType == INT_ADAPTIVELS) {
		mxAddNumMatrix(mxrwsstruct, "lsAdapt", grampc->rws->lsAdapt, 1, 2 * (NALS + 1)*(1 + grampc->opt->MaxGradIter));
		mxAddNumMatrix(mxrwsstruct, "lsExplicit", grampc->rws->lsExplicit, 0, 0);
	}
	else {
		mxAddNumMatrix(mxrwsstruct, "lsAdapt", grampc->rws->lsAdapt, 0, 0);
		mxAddNumMatrix(mxrwsstruct, "lsExplicit", grampc->rws->lsExplicit, 1, NELS);
	}
	mxAddNumMatrix(mxrwsstruct, "rwsScale", grampc->rws->rwsScale, 1, 2 * (grampc->param->Nx + grampc->param->Nu + grampc->param->Np));
	mxAddIntMatrix(mxrwsstruct, "lrwsGeneral", &grampc->rws->lrwsGeneral, 1, 1);
	mxAddNumMatrix(mxrwsstruct, "rwsGeneral", grampc->rws->rwsGeneral, 1, grampc->rws->lrwsGeneral);

	mxAddIntMatrix(mxrwsstruct, "lworkRodas", &grampc->rws->lworkRodas, 1, 1);
	mxAddIntMatrix(mxrwsstruct, "liworkRodas", &grampc->rws->liworkRodas, 1, 1);
	if (grampc->opt->Integrator == INT_RODAS) {
		mxAddNumMatrix(mxrwsstruct, "rparRodas", grampc->rws->rparRodas, 1, grampc->param->Nx*grampc->opt->Nhor);
		mxAddIntMatrix(mxrwsstruct, "iparRodas", grampc->rws->iparRodas, 1, 20);
		mxAddNumMatrix(mxrwsstruct, "workRodas", grampc->rws->workRodas, 1, grampc->rws->lworkRodas);
		mxAddIntMatrix(mxrwsstruct, "iworkRodas", grampc->rws->iworkRodas, 1, grampc->rws->liworkRodas);
	}
	else {
		mxAddNumMatrix(mxrwsstruct, "rparRodas", grampc->rws->rparRodas, 0, 0);
		mxAddIntMatrix(mxrwsstruct, "iparRodas", grampc->rws->iparRodas, 0, 0);
		mxAddNumMatrix(mxrwsstruct, "workRodas", grampc->rws->workRodas, 0, 0);
		mxAddIntMatrix(mxrwsstruct, "iworkRodas", grampc->rws->iworkRodas, 0, 0);
	}

	mxAddStructField(mxgrampcstruct, "rws", mxrwsstruct);

	/* USERPARAM **************************************************************/
	mxAddStructField(mxgrampcstruct, "userparam", mxuserparam);

	plhs[0] = mxgrampcstruct;
	/**************************************************************************/
}
#endif

typeBoolean AssignRealmxInput(typeRNum **x, const mxArray *prhs[], ctypeInt inputNo) {

	if (mxGetClassID(prhs[inputNo]) == mxtypeRNum_CLASS) {
		*x = (typeRNum*)mxGetPr(prhs[inputNo]);
		return 0;
	}
	else if (mxGetClassID(prhs[inputNo]) == mxDOUBLE_CLASS) {
		CastDvec2Numvec(x, mxGetPr(prhs[inputNo]), mxGetM(prhs[inputNo])*mxGetN(prhs[inputNo]));
		return 1;
	}
	else {
		printf("Invalid data type of input argument %d", inputNo);
		printError("Invalid data type of input argument");
		return 0;
	}
}

void AssignRealmxOutput(typeRNum **output, double**doutput, mxArray *plhs[], ctypeInt cast, ctypeInt  outputSize, ctypeInt outputNo) {
	if (cast) {
		plhs[outputNo] = mxCreateNumericMatrix(outputSize, 1, mxDOUBLE_CLASS, mxREAL);
		createNumMatrix(output, outputSize);
		*doutput = (double*)mxGetPr(plhs[0]);
	}
	else {
		plhs[outputNo] = mxCreateNumericMatrix(outputSize, 1, mxtypeRNum_CLASS, mxREAL);
		*output = (typeRNum*)mxGetPr(plhs[0]);
		*doutput = NULL;
	}
}
#endif /* !GRAMPC_CONVERSION_CMEX_H_*/
