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
		resizeNumMatrix(&grampc->rws->lsAdapt, 2 * (NALS + 1)*(1 + grampc->opt->MaxGradIter));
		resizeNumMatrix(&grampc->rws->lsExplicit, 0);
		for (i = 0; i < grampc->opt->MaxGradIter + 1; i++) {
			grampc->rws->lsAdapt[0 + i * 2 * (NALS + 1)] = grampc->opt->LineSearchInit*(1 - grampc->opt->LineSearchIntervalFactor);
			grampc->rws->lsAdapt[1 + i * 2 * (NALS + 1)] = grampc->opt->LineSearchInit;
			grampc->rws->lsAdapt[2 + i * 2 * (NALS + 1)] = grampc->opt->LineSearchInit*(1 + grampc->opt->LineSearchIntervalFactor);
			grampc->rws->lsAdapt[3 + i * 2 * (NALS + 1)] = grampc->opt->LineSearchInit;
		}
	}
	else {
		resizeNumMatrix(&grampc->rws->lsAdapt, 0);
		resizeNumMatrix(&grampc->rws->lsExplicit, NELS);
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


	/* USERPARAM ****************************************************************/
	(*grampc)->userparam = userparam;


	/* PARAM STRUCTURE **********************************************************/
	ocp_dim(&(*grampc)->param->Nx, &(*grampc)->param->Nu, &(*grampc)->param->Np, &(*grampc)->param->Ng, &(*grampc)->param->Nh, &(*grampc)->param->NgT, &(*grampc)->param->NhT, (*grampc)->userparam);
	/* check number of states and control */
	if ((*grampc)->param->Nx <= 0) {
		grampc_error(INVALID_NX);
	}
	if ((*grampc)->param->Nu < 0) {
		grampc_error(INVALID_NU);
	}
	if ((*grampc)->param->Np < 0) {
		grampc_error(INVALID_NP);
	}
	if ((*grampc)->param->Ng < 0) {
		grampc_error(INVALID_Ng);
	}
	if ((*grampc)->param->Nh < 0) {
		grampc_error(INVALID_Nh);
	}
	if ((*grampc)->param->NgT < 0) {
		grampc_error(INVALID_NgT);
	}
	if ((*grampc)->param->NhT < 0) {
		grampc_error(INVALID_NhT);
	}
	(*grampc)->param->Nc = (*grampc)->param->Ng + (*grampc)->param->Nh + (*grampc)->param->NgT + (*grampc)->param->NhT;

	createNumMatrix(&(*grampc)->param->x0, (*grampc)->param->Nx);
	createNumMatrix(&(*grampc)->param->xdes, (*grampc)->param->Nx);

	createNumMatrix(&(*grampc)->param->u0, (*grampc)->param->Nu);
	createNumMatrix(&(*grampc)->param->udes, (*grampc)->param->Nu);
	createNumMatrix(&(*grampc)->param->umax, (*grampc)->param->Nu);
	createNumMatrix(&(*grampc)->param->umin, (*grampc)->param->Nu);
	MatSetScalar((*grampc)->param->umax, INF, 1, (*grampc)->param->Nu);
	MatSetScalar((*grampc)->param->umin, -INF, 1, (*grampc)->param->Nu);

	createNumMatrix(&(*grampc)->param->p0, (*grampc)->param->Np);
	createNumMatrix(&(*grampc)->param->pmax, (*grampc)->param->Np);
	createNumMatrix(&(*grampc)->param->pmin, (*grampc)->param->Np);
	MatSetScalar((*grampc)->param->pmax, INF, 1, (*grampc)->param->Np);
	MatSetScalar((*grampc)->param->pmin, -INF, 1, (*grampc)->param->Np);

	(*grampc)->param->Thor = (typeRNum)-1.0;
	(*grampc)->param->Tmax = (typeRNum)1e8;
	(*grampc)->param->Tmin = (typeRNum)1e-8;

	(*grampc)->param->dt = (typeRNum)-1.0;
	(*grampc)->param->t0 = (typeRNum)0.0;


	/* OPTIONS STRUCTURE *******************************************************/
	(*grampc)->opt->Nhor = 30;
	(*grampc)->opt->MaxGradIter = 2;
	(*grampc)->opt->MaxMultIter = 1;
	(*grampc)->opt->ShiftControl = INT_ON;

	(*grampc)->opt->TimeDiscretization = INT_UNIFORM;

	(*grampc)->opt->IntegralCost = INT_ON;
	(*grampc)->opt->TerminalCost = INT_ON;
	(*grampc)->opt->IntegratorCost = INT_TRAPZ;

	(*grampc)->opt->Integrator = INT_HEUN;
	(*grampc)->opt->IntegratorRelTol = (typeRNum)1e-6;
	(*grampc)->opt->IntegratorAbsTol = (typeRNum)1e-8;
	(*grampc)->opt->IntegratorMinStepSize = EPS;
	(*grampc)->opt->IntegratorMaxSteps = (typeInt)1e8;
	createIntMatrix(&(*grampc)->opt->FlagsRodas, 8);
	(*grampc)->opt->FlagsRodas[0] = 0;                     /* 0 --> right hand side independent of time t */
	(*grampc)->opt->FlagsRodas[1] = 0;                     /* 0 --> DF/DX is numerically computed */
	(*grampc)->opt->FlagsRodas[2] = 0;                     /* 1(0) -> analytical (numerical) jacobian (partial derivatives of right hand side w.r.t. state) */
	(*grampc)->opt->FlagsRodas[4] = (*grampc)->param->Nx;  /* no. of lower diagonals of jacobian */
	(*grampc)->opt->FlagsRodas[5] = (*grampc)->param->Nx;  /* no. of upper diagonals of jacobian */
	(*grampc)->opt->FlagsRodas[3] = 0;                     /* 1 --> mass matrix is supplied */
	(*grampc)->opt->FlagsRodas[6] = (*grampc)->param->Nx;  /* no. of lower diagonals of mass matrix */
	(*grampc)->opt->FlagsRodas[7] = (*grampc)->param->Nx;  /* no. of upper diagonals of mass matrix */

	(*grampc)->opt->LineSearchType = INT_EXPLS2;
	(*grampc)->opt->LineSearchExpAutoFallback = INT_ON;
	(*grampc)->opt->LineSearchMax = (typeRNum)0.75;
	(*grampc)->opt->LineSearchMin = (typeRNum)1e-10;
	(*grampc)->opt->LineSearchInit = (typeRNum)1e-4;
	(*grampc)->opt->LineSearchIntervalFactor = (typeRNum)0.85;
	(*grampc)->opt->LineSearchAdaptFactor = (typeRNum)3.0 / 2.0;
	(*grampc)->opt->LineSearchIntervalTol = (typeRNum)1e-1;

	(*grampc)->opt->OptimControl = INT_ON;
	(*grampc)->opt->OptimParam = INT_OFF;
	(*grampc)->opt->OptimParamLineSearchFactor = (typeRNum)1.0;
	(*grampc)->opt->OptimTime = INT_OFF;
	(*grampc)->opt->OptimTimeLineSearchFactor = (typeRNum)1.0;

	(*grampc)->opt->ScaleProblem = INT_OFF;
	createNumMatrix(&(*grampc)->opt->xScale, (*grampc)->param->Nx);
	createNumMatrix(&(*grampc)->opt->xOffset, (*grampc)->param->Nx);
	MatSetScalar((*grampc)->opt->xScale, 1.0, 1, (*grampc)->param->Nx);

	createNumMatrix(&(*grampc)->opt->uScale, (*grampc)->param->Nu);
	createNumMatrix(&(*grampc)->opt->uOffset, (*grampc)->param->Nu);
	MatSetScalar((*grampc)->opt->uScale, 1.0, 1, (*grampc)->param->Nu);

	createNumMatrix(&(*grampc)->opt->pScale, (*grampc)->param->Np);
	createNumMatrix(&(*grampc)->opt->pOffset, (*grampc)->param->Np);
	MatSetScalar((*grampc)->opt->pScale, 1.0, 1, (*grampc)->param->Np);

	(*grampc)->opt->TScale = (typeRNum)1.0;
	(*grampc)->opt->TOffset = (typeRNum)0.0;

	(*grampc)->opt->JScale = (typeRNum)1.0;

	createNumMatrix(&(*grampc)->opt->cScale, (*grampc)->param->Nc);
	MatSetScalar((*grampc)->opt->cScale, 1.0, 1, (*grampc)->param->Nc);

	(*grampc)->opt->EqualityConstraints = INT_ON;
	(*grampc)->opt->InequalityConstraints = INT_ON;
	(*grampc)->opt->TerminalEqualityConstraints = INT_ON;
	(*grampc)->opt->TerminalInequalityConstraints = INT_ON;
	(*grampc)->opt->ConstraintsHandling = INT_AUGLAG;
	createNumMatrix(&(*grampc)->opt->ConstraintsAbsTol, (*grampc)->param->Nc);
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

	/* SOLUTION STRUCTURE ******************************************************/
	createNumMatrix(&(*grampc)->sol->xnext, (*grampc)->param->Nx);
	createNumMatrix(&(*grampc)->sol->unext, (*grampc)->param->Nu);
	createNumMatrix(&(*grampc)->sol->pnext, (*grampc)->param->Np);
	(*grampc)->sol->Tnext = (typeRNum)0.0;
	createNumMatrix(&(*grampc)->sol->J, 2);
	(*grampc)->sol->cfct = (typeRNum)0.0;
	(*grampc)->sol->pen = (typeRNum)0.0;
	createIntMatrix(&(*grampc)->sol->iter, (*grampc)->opt->MaxMultIter);


	/* RWS STRUCTURE **********************************************************/
	createNumMatrix(&(*grampc)->rws->t, (*grampc)->opt->Nhor);
	createNumMatrix(&(*grampc)->rws->tls, (*grampc)->opt->Nhor);
	init_rws_time(*grampc);

	createNumMatrix(&(*grampc)->rws->x, (*grampc)->param->Nx*(*grampc)->opt->Nhor);
	createNumMatrix(&(*grampc)->rws->adj, (*grampc)->param->Nx*(*grampc)->opt->Nhor);
	createNumMatrix(&(*grampc)->rws->dcdx, (*grampc)->param->Nx*((*grampc)->opt->Nhor + 1));

	createNumMatrix(&(*grampc)->rws->u, (*grampc)->param->Nu*(*grampc)->opt->Nhor);
	createNumMatrix(&(*grampc)->rws->uls, (*grampc)->param->Nu*(*grampc)->opt->Nhor);
	createNumMatrix(&(*grampc)->rws->uprev, (*grampc)->param->Nu*(*grampc)->opt->Nhor);
	createNumMatrix(&(*grampc)->rws->gradu, (*grampc)->param->Nu*(*grampc)->opt->Nhor);
	createNumMatrix(&(*grampc)->rws->graduprev, (*grampc)->param->Nu*(*grampc)->opt->Nhor);
	createNumMatrix(&(*grampc)->rws->dcdu, (*grampc)->param->Nu*(*grampc)->opt->Nhor);

	createNumMatrix(&(*grampc)->rws->p, (*grampc)->param->Np);
	createNumMatrix(&(*grampc)->rws->pls, (*grampc)->param->Np);
	createNumMatrix(&(*grampc)->rws->pprev, (*grampc)->param->Np);
	createNumMatrix(&(*grampc)->rws->gradp, (*grampc)->param->Np);
	createNumMatrix(&(*grampc)->rws->gradpprev, (*grampc)->param->Np);
	createNumMatrix(&(*grampc)->rws->dcdp, (*grampc)->param->Np * ((*grampc)->opt->Nhor + 1));

	(*grampc)->rws->gradT = 0;
	(*grampc)->rws->gradTprev = 0;
	(*grampc)->rws->dcdt = 0;

	createNumMatrix(&(*grampc)->rws->mult, (*grampc)->param->Nc*(*grampc)->opt->Nhor);
	createNumMatrix(&(*grampc)->rws->pen, (*grampc)->param->Nc*(*grampc)->opt->Nhor);
	createNumMatrix(&(*grampc)->rws->cfct, (*grampc)->param->Nc*(*grampc)->opt->Nhor);
	createNumMatrix(&(*grampc)->rws->cfctprev, (*grampc)->param->Nc*(*grampc)->opt->Nhor);
	createNumMatrix(&(*grampc)->rws->cfctAbsTol, (*grampc)->param->Nc);
	init_rws_multipliers(*grampc);

	init_rws_linesearch(*grampc);

	createNumMatrix(&(*grampc)->rws->rwsScale, 2 * ((*grampc)->param->Nx + (*grampc)->param->Nu + (*grampc)->param->Np));
	resize_rwsGeneral(*grampc);

	(*grampc)->rws->lworkRodas = 2 * (*grampc)->param->Nx*(*grampc)->param->Nx + 14 * (*grampc)->param->Nx + 20;
	(*grampc)->rws->liworkRodas = (*grampc)->param->Nx + 20;
}

void grampc_free(typeGRAMPC **grampc)
{
	free((*grampc)->param->x0);
	free((*grampc)->param->xdes);
	free((*grampc)->param->u0);
	free((*grampc)->param->udes);
	free((*grampc)->param->umax);
	free((*grampc)->param->umin);
	free((*grampc)->param->p0);
	free((*grampc)->param->pmax);
	free((*grampc)->param->pmin);

	free((*grampc)->opt->FlagsRodas);
	free((*grampc)->opt->xScale);
	free((*grampc)->opt->xOffset);
	free((*grampc)->opt->uScale);
	free((*grampc)->opt->uOffset);
	free((*grampc)->opt->pScale);
	free((*grampc)->opt->pOffset);
	free((*grampc)->opt->cScale);
	free((*grampc)->opt->ConstraintsAbsTol);

	free((*grampc)->sol->xnext);
	free((*grampc)->sol->unext);
	free((*grampc)->sol->pnext);
	free((*grampc)->sol->J);
	free((*grampc)->sol->iter);

	free((*grampc)->rws->t);
	free((*grampc)->rws->tls);
	free((*grampc)->rws->x);
	free((*grampc)->rws->adj);
	free((*grampc)->rws->dcdx);
	free((*grampc)->rws->u);
	free((*grampc)->rws->uls);
	free((*grampc)->rws->uprev);
	free((*grampc)->rws->gradu);
	free((*grampc)->rws->graduprev);
	free((*grampc)->rws->dcdu);
	free((*grampc)->rws->p);
	free((*grampc)->rws->pls);
	free((*grampc)->rws->pprev);
	free((*grampc)->rws->gradp);
	free((*grampc)->rws->gradpprev);
	free((*grampc)->rws->dcdp);
	free((*grampc)->rws->mult);
	free((*grampc)->rws->pen);
	free((*grampc)->rws->cfct);
	free((*grampc)->rws->cfctprev);
	free((*grampc)->rws->cfctAbsTol);
	free((*grampc)->rws->lsAdapt);
	free((*grampc)->rws->lsExplicit);
	free((*grampc)->rws->rwsScale);
	free((*grampc)->rws->rwsGeneral);
	free((*grampc)->rws->rparRodas);
	free((*grampc)->rws->iparRodas);
	free((*grampc)->rws->workRodas);
	free((*grampc)->rws->iworkRodas);

	free((*grampc)->param);
	free((*grampc)->opt);
	free((*grampc)->sol);
	free((*grampc)->rws);

	free(*grampc);
}
