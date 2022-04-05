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


#include "grampc_alloc.h"
#include "grampc_init.h"
#include "grampc_mess.h"
#include "grampc_util.h"
#include "probfct.h"

void init_rws_time(const typeGRAMPC *grampc) {
	if (grampc->opt->ScaleProblem == INT_ON) {
		scale_time(&grampc->rws->T, grampc->param->Thor, grampc);
	}
	else {
		grampc->rws->T = grampc->param->Thor;
	}
	grampc->rws->Tprev = grampc->rws->T;
	discretize_time(grampc->rws->t, grampc->rws->T, grampc);
}

void init_rws_controls(const typeGRAMPC *grampc) {
	typeInt i, k;
	for (i = 0; i < grampc->opt->Nhor; i++) {
		k = i * grampc->param->Nu;
		if (grampc->opt->ScaleProblem == INT_ON) {
			scale_controls(grampc->rws->u + k, grampc->param->u0, grampc);
		}
		else {
			MatCopy(grampc->rws->u + k, grampc->param->u0, 1, grampc->param->Nu);
		}
	}
	MatCopy(grampc->rws->uprev, grampc->rws->u, grampc->opt->Nhor, grampc->param->Nu);
}

void init_rws_parameters(const typeGRAMPC *grampc) {
	if (grampc->opt->ScaleProblem == INT_ON) {
		scale_parameters(grampc->rws->p, grampc->param->p0, grampc);
	}
	else {
		MatCopy(grampc->rws->p, grampc->param->p0, 1, grampc->param->Np);
	}
	MatCopy(grampc->rws->pprev, grampc->rws->p, 1, grampc->param->Np);
}

void init_rws_linesearch(const typeGRAMPC *grampc) {
	typeInt i;
    if (grampc->opt->LineSearchType == INT_ADAPTIVELS) {
        /* init lsAdapt */
		for (i = 0; i < grampc->opt->MaxGradIter + 1; i++) {
			grampc->rws->lsAdapt[0 + i * 2 * (NALS + 1)] = grampc->opt->LineSearchInit*(1 - grampc->opt->LineSearchIntervalFactor);
			grampc->rws->lsAdapt[1 + i * 2 * (NALS + 1)] = grampc->opt->LineSearchInit;
			grampc->rws->lsAdapt[2 + i * 2 * (NALS + 1)] = grampc->opt->LineSearchInit*(1 + grampc->opt->LineSearchIntervalFactor);
			grampc->rws->lsAdapt[3 + i * 2 * (NALS + 1)] = grampc->opt->LineSearchInit;
		}
	}
    else {
		/* init lsExplicit */
		grampc->rws->lsExplicit[2] = grampc->opt->LineSearchInit;
		check_ControlLimits(grampc);
	}
}

void init_rws_multipliers(const typeGRAMPC *grampc) {
	MatSetScalar(grampc->rws->pen, grampc->opt->PenaltyMin, grampc->param->Nc, grampc->opt->Nhor);
}

void init_rws_constraints(const typeGRAMPC *grampc) {
	MatCopy(grampc->rws->cfctAbsTol, grampc->opt->ConstraintsAbsTol, 1, grampc->param->Nc);
	if (grampc->opt->ScaleProblem == INT_ON) {
		scale_constraints(grampc->rws->cfctAbsTol, grampc->opt->cScale, grampc->param->Nc);
	}
}

void grampc_init(typeGRAMPC **grampc, typeUSERPARAM *userparam)
{
    /* allocate memory */
    grampc_alloc_structs(grampc, userparam);

    /* initialize options that cannot be changed in fixed-size mode and that are relevant for allocation */
    #ifdef FIXEDSIZE
        (*grampc)->opt->Nhor = NHOR;
        (*grampc)->opt->MaxGradIter = MAXGRADITER;
        (*grampc)->opt->MaxMultIter = MAXMULTITER;
        (*grampc)->opt->Integrator = INTEGRATOR;
        (*grampc)->opt->IntegratorCost = INTEGRATORCOST;
        (*grampc)->opt->LineSearchType = LINESEARCHTYPE;
    #else
        (*grampc)->opt->Nhor = 30;
        (*grampc)->opt->MaxGradIter = 2;
        (*grampc)->opt->MaxMultIter = 1;
        (*grampc)->opt->Integrator = INT_HEUN;
        (*grampc)->opt->IntegratorCost = INT_TRAPZ;
        (*grampc)->opt->LineSearchType = INT_EXPLS2;

        /* alloc memory for line search-dependent fields */
        resize_rwsLinesearch(*grampc);
    #endif

    /* allocate memory */
    grampc_alloc_fields(grampc, userparam);

    /* set default parameters */
	MatSetScalar((*grampc)->param->umax, INF, 1, (*grampc)->param->Nu);
	MatSetScalar((*grampc)->param->umin, -INF, 1, (*grampc)->param->Nu);

	MatSetScalar((*grampc)->param->pmax, INF, 1, (*grampc)->param->Np);
	MatSetScalar((*grampc)->param->pmin, -INF, 1, (*grampc)->param->Np);

	(*grampc)->param->Thor = (typeRNum)-1.0;
	(*grampc)->param->Tmax = (typeRNum)1e8;
	(*grampc)->param->Tmin = (typeRNum)1e-8;

	(*grampc)->param->dt = (typeRNum)-1.0;
	(*grampc)->param->t0 = (typeRNum)0.0;

    /* initialize options */
	(*grampc)->opt->ShiftControl = INT_ON;

	(*grampc)->opt->TimeDiscretization = INT_UNIFORM;

	(*grampc)->opt->IntegralCost = INT_ON;
    (*grampc)->opt->TerminalCost = INT_ON;

	(*grampc)->opt->IntegratorRelTol = (typeRNum)1e-6;
	(*grampc)->opt->IntegratorAbsTol = (typeRNum)1e-8;
	(*grampc)->opt->IntegratorMinStepSize = EPS;
    (*grampc)->opt->IntegratorMaxSteps = (typeInt)1e8;

#if defined(FIXEDSIZE)
    /* rodas-specific options cannot be changed in fixed-size mode */
    (*grampc)->opt->FlagsRodas[0] = RODAS_IFCN;
    (*grampc)->opt->FlagsRodas[1] = RODAS_IDFX;
    (*grampc)->opt->FlagsRodas[2] = RODAS_IJAC;
    (*grampc)->opt->FlagsRodas[3] = RODAS_IMAS;
    (*grampc)->opt->FlagsRodas[4] = RODAS_MLJAC;
    (*grampc)->opt->FlagsRodas[5] = RODAS_MUJAC;
    (*grampc)->opt->FlagsRodas[6] = RODAS_MLMAS;
    (*grampc)->opt->FlagsRodas[7] = RODAS_MUMAS;

    (*grampc)->rws->lworkRodas = SIZE_RWORKRODAS;
    (*grampc)->rws->liworkRodas = SIZE_IWORKRODAS;
#else
	(*grampc)->opt->FlagsRodas[0] = 0;                     /* 0 --> right hand side independent of time t */
	(*grampc)->opt->FlagsRodas[1] = 0;                     /* 0 --> DF/DX is numerically computed */
	(*grampc)->opt->FlagsRodas[2] = 0;                     /* 1(0) -> analytical (numerical) jacobian (partial derivatives of right hand side w.r.t. state) */
	(*grampc)->opt->FlagsRodas[4] = (*grampc)->param->Nx;  /* no. of lower diagonals of jacobian */
	(*grampc)->opt->FlagsRodas[5] = (*grampc)->param->Nx;  /* no. of upper diagonals of jacobian */
	(*grampc)->opt->FlagsRodas[3] = 0;                     /* 1 --> mass matrix is supplied */
	(*grampc)->opt->FlagsRodas[6] = (*grampc)->param->Nx;  /* no. of lower diagonals of mass matrix */
	(*grampc)->opt->FlagsRodas[7] = (*grampc)->param->Nx;  /* no. of upper diagonals of mass matrix */

    (*grampc)->rws->lworkRodas = 2 * (*grampc)->param->Nx*(*grampc)->param->Nx + 14 * (*grampc)->param->Nx + 20;
    (*grampc)->rws->liworkRodas = (*grampc)->param->Nx + 20;

    /* alloc memory for rodas-dependent fields */
    resize_rwsRodas(*grampc);
#endif

	(*grampc)->opt->LineSearchExpAutoFallback = INT_ON;
	(*grampc)->opt->LineSearchMax = (typeRNum)0.75;
	(*grampc)->opt->LineSearchMin = (typeRNum)1e-10;
	(*grampc)->opt->LineSearchInit = (typeRNum)1e-4;
	(*grampc)->opt->LineSearchAdaptAbsTol = (typeRNum)1e-6;
	(*grampc)->opt->LineSearchAdaptFactor = (typeRNum)3.0 / 2.0;
	(*grampc)->opt->LineSearchIntervalTol = (typeRNum)1e-1;
	(*grampc)->opt->LineSearchIntervalFactor = (typeRNum)0.85;

	(*grampc)->opt->OptimControl = INT_ON;
	(*grampc)->opt->OptimParam = INT_OFF;
	(*grampc)->opt->OptimParamLineSearchFactor = (typeRNum)1.0;
	(*grampc)->opt->OptimTime = INT_OFF;
	(*grampc)->opt->OptimTimeLineSearchFactor = (typeRNum)1.0;

    (*grampc)->opt->ScaleProblem = INT_OFF;
    MatSetScalar((*grampc)->opt->xScale, 1.0, 1, (*grampc)->param->Nx);
    MatSetScalar((*grampc)->opt->uScale, 1.0, 1, (*grampc)->param->Nu);
	MatSetScalar((*grampc)->opt->pScale, 1.0, 1, (*grampc)->param->Np);
	(*grampc)->opt->TScale = (typeRNum)1.0;
	(*grampc)->opt->TOffset = (typeRNum)0.0;
	(*grampc)->opt->JScale = (typeRNum)1.0;
	MatSetScalar((*grampc)->opt->cScale, 1.0, 1, (*grampc)->param->Nc);

	(*grampc)->opt->EqualityConstraints = INT_ON;
	(*grampc)->opt->InequalityConstraints = INT_ON;
	(*grampc)->opt->TerminalEqualityConstraints = INT_ON;
	(*grampc)->opt->TerminalInequalityConstraints = INT_ON;
    (*grampc)->opt->ConstraintsHandling = INT_AUGLAG;
	MatSetScalar((*grampc)->opt->ConstraintsAbsTol, (typeRNum) 1e-4, 1, (*grampc)->param->Nc);

	(*grampc)->opt->MultiplierMax = (typeRNum)1e6;
	(*grampc)->opt->MultiplierDampingFactor = (typeRNum)0.0;
	(*grampc)->opt->PenaltyMax = (typeRNum)1e6;
	(*grampc)->opt->PenaltyMin = (typeRNum)1e0;
	(*grampc)->opt->PenaltyIncreaseFactor = (typeRNum)1.05;
	(*grampc)->opt->PenaltyDecreaseFactor = (typeRNum)0.95;
	(*grampc)->opt->PenaltyIncreaseThreshold = (typeRNum)1.0;
	(*grampc)->opt->AugLagUpdateGradientRelTol = (typeRNum)1e-2;

	(*grampc)->opt->ConvergenceCheck = INT_OFF;
	(*grampc)->opt->ConvergenceGradientRelTol = (typeRNum)1e-6;

    /* initialize sol */
    (*grampc)->sol->Tnext = (typeRNum)0.0;
	(*grampc)->sol->cfct = (typeRNum)0.0;
    (*grampc)->sol->pen = (typeRNum)0.0;

    /* initialize rws */
	init_rws_time(*grampc);
    init_rws_multipliers(*grampc);
    init_rws_linesearch(*grampc);

	(*grampc)->rws->gradT = 0;
	(*grampc)->rws->gradTprev = 0;
	(*grampc)->rws->dcdt = 0;

#ifdef FIXEDSIZE
    (*grampc)->rws->lrwsGeneral = SIZE_RWSGENERAL;
#else
    /* alloc memory for general workspace */
	resize_rwsGeneral(*grampc);
#endif
}
