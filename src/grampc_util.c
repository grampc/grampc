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


#include "grampc_util.h"
#include "grampc_alloc.h"


void discretize_time(typeRNum *tvec, typeRNum T, const typeGRAMPC *grampc)
{
	typeInt i;
	typeInt N = grampc->opt->Nhor;
	typeRNum m, c;
	typeRNum T_ = T;

	if (grampc->opt->ScaleProblem == INT_ON) {
		unscale_time(&T_, T, grampc);
	}

	/* Uniform discretization of interval [0, T_] */
	if (grampc->opt->TimeDiscretization == INT_UNIFORM || T_ <= (N - 1) * grampc->param->dt) {
		m = T_ / (N - 1);
		for (i = 0; i < grampc->opt->Nhor; i++) {
			tvec[i] = m * i;
		}
	}
	/* Nonuniform discretization of interval [0, T_] */
	else if (grampc->opt->TimeDiscretization == INT_NONUNIFORM) {
		m = (T_ / (N - 1) - grampc->param->dt) / (N - 2);
		c = grampc->param->dt - m;
		for (i = 0; i < grampc->opt->Nhor; i++) {
			tvec[i] = m * i * i + c * i;
		}
	}
}

void check_ControlLimits(const typeGRAMPC* grampc) {
	typeInt i;
	/* if any control limit is INF or -INF disable auto fallback strategy of the explicit linesearch */
	if (grampc->opt->LineSearchType == INT_EXPLS1 || grampc->opt->LineSearchType == INT_EXPLS2) {
		grampc->rws->lsExplicit[3] = 1;

		for (i = 0; i < grampc->param->Nu; i++) {
			if (grampc->param->umax[i] >= INF || grampc->param->umin[i] <= -INF) {
				grampc->rws->lsExplicit[3] = 0;
			}
		}
	}
}

char* IntegratorInt2Str(ctypeInt INT_Integrator) {
	switch (INT_Integrator) {
	case INT_EULER:    return "euler";
	case INT_MODEULER: return "modeuler";
	case INT_HEUN:     return "heun";
	case INT_RODAS:    return "rodas";
	case INT_RUKU45:   return "ruku45";
	}
	return "";
}

char* LineSearchTypeInt2Str(ctypeInt INT_LineSearchType) {
	switch (INT_LineSearchType) {
	case INT_ADAPTIVELS: return "adaptive";
	case INT_EXPLS1:     return "explicit1";
	case INT_EXPLS2:     return "explicit2";
	}
	return "";
}


void lsearch_fit(typeRNum *kfit, typeRNum *Jfit, ctypeRNum *k, ctypeRNum *J)
{
	ctypeRNum aux = (2 * J[1] - J[0] - J[2]);
	ctypeRNum aux2 = (J[2] - 4 * J[1]);

	/* exists a minimum? (positive curvature) */
	if (aux <= -aEPS) {
		kfit[0] = k[1] + (k[1] - k[0]) / 2 * (J[2] - J[0]) / aux;
		Jfit[0] = (J[0] * J[0] + aux2 * aux2 - 2 * J[0] * (4 * J[1] + J[2])) / (8 * aux);
	}
	/* smallest J */
	if (aux > -aEPS  || kfit[0] < k[0] || kfit[0] > k[2]) {
		if (J[0] <= J[1] && J[0] <= J[2]) {
			kfit[0] = k[0];
			Jfit[0] = J[0];
		}
		else if (J[2] <= J[0]  && J[2] <= J[1]) {
			kfit[0] = k[2];
			Jfit[0] = J[2];
		}
		else {
			kfit[0] = k[1];
			Jfit[0] = J[1];
		}
	}
}

void interplin(typeRNum *varint, ctypeRNum *tvec, ctypeRNum *varvec, ctypeRNum tint,
	ctypeInt Nvar, ctypeInt Nvec, ctypeInt searchdir)
{
	/* option to determine position ioff in time vector such that tvec[ioff] <= tint <= tvec[ioff+1]
	 * searchdir= 1: search in forward direction starting at ioff=0									 *
	 * searchdir=-1: search in backward direction starting at ioff=Nvec-1
	 */
	typeInt i;
	typeInt ioff;
	typeRNum dtratio;
	ctypeRNum *var0, *var1;

	if (tint <= tvec[0]) {
		for (i = 0; i < Nvar; i++) {
			varint[i] = varvec[i];
		}
	}
	else if (tint >= tvec[Nvec - 1]) {
		var0 = varvec + (Nvec - 1)*Nvar;
		for (i = 0; i < Nvar; i++)
			varint[i] = var0[i];
	}
	else {
		if (searchdir == 1) {
			ioff = 0;
			while (tvec[ioff] < tint) {
				ioff += 1;
			}
			ioff -= 1;
		}
		else {
			ioff = Nvec - 2;
			while (tvec[ioff] > tint) {
				ioff -= 1;
			}
		}
		dtratio = (tint - tvec[ioff]) / (tvec[ioff + 1] - tvec[ioff]);
		var0 = varvec + ioff * Nvar;
		var1 = var0 + Nvar;
		for (i = 0; i < Nvar; i++) {
			varint[i] = var0[i] + dtratio * (var1[i] - var0[i]);
		}
	}
}

void unscale_states(typeRNum *out, ctypeRNum *x, const typeGRAMPC *grampc)
{
	typeInt i;
	for (i = 0; i < grampc->param->Nx; i++) {
		out[i] = x[i] * grampc->opt->xScale[i] + grampc->opt->xOffset[i];
	}
}

void unscale_adjoints(typeRNum *out, ctypeRNum *adj, const typeGRAMPC *grampc)
{
	typeInt i;
	for (i = 0; i < grampc->param->Nx; i++) {
		out[i] = adj[i] / grampc->opt->xScale[i];
	}
}

void unscale_controls(typeRNum *out, ctypeRNum *u, const typeGRAMPC *grampc)
{
	typeInt i;
	for (i = 0; i < grampc->param->Nu; i++) {
		out[i] = u[i] * grampc->opt->uScale[i] + grampc->opt->uOffset[i];
	}
}

void unscale_parameters(typeRNum *out, ctypeRNum *p, const typeGRAMPC *grampc)
{
	typeInt i;
	for (i = 0; i < grampc->param->Np; i++) {
		out[i] = p[i] * grampc->opt->pScale[i] + grampc->opt->pOffset[i];
	}
}

void unscale_time(typeRNum *out, ctypeRNum T, const typeGRAMPC *grampc)
{
	out[0] = T * grampc->opt->TScale + grampc->opt->TOffset;
}


void scale_states(typeRNum *out, ctypeRNum *x, const typeGRAMPC *grampc)
{
	typeInt i;
	for (i = 0; i < grampc->param->Nx; i++) {
		out[i] = (x[i] - grampc->opt->xOffset[i]) / grampc->opt->xScale[i];
	}
}

void scale_adjoints(typeRNum *out, ctypeRNum *adj, const typeGRAMPC *grampc)
{
	typeInt i;
	for (i = 0; i < grampc->param->Nx; i++) {
		out[i] = adj[i] * grampc->opt->xScale[i];
	}
}

void scale_controls(typeRNum *out, ctypeRNum *u, const typeGRAMPC *grampc)
{
	typeInt i;
	for (i = 0; i < grampc->param->Nu; i++) {
		out[i] = (u[i] - grampc->opt->uOffset[i]) / grampc->opt->uScale[i];
	}
}

void scale_parameters(typeRNum *out, ctypeRNum *p, const typeGRAMPC *grampc)
{
	typeInt i;
	for (i = 0; i < grampc->param->Np; i++) {
		out[i] = (p[i] - grampc->opt->pOffset[i]) / grampc->opt->pScale[i];
	}
}

void scale_time(typeRNum *out, typeRNum T, const typeGRAMPC *grampc)
{
	out[0] = (T - grampc->opt->TOffset) / grampc->opt->TScale;
}

void scale_constraints(typeRNum *c, ctypeRNum *cScale, ctypeInt Ncon)
{
	typeInt i;
	for (i = 0; i < Ncon; i++) {
		c[i] = c[i] / cScale[i];
	}
}

void scale_cost(typeRNum *J, ctypeRNum JScale, ctypeInt Ncost)
{
	typeInt i;
	for (i = 0; i < Ncost; i++) {
		J[i] = J[i] / JScale;
	}
}

void MatSetScalar(typeRNum *C, ctypeRNum a, ctypeInt n1, ctypeInt n2)
{
	/* matrix set C = a *
	 * C: (n1 x n2)     */
	typeInt i;
	for (i = 0; i < n1*n2; i++) {
		C[i] = a;
	}
}

void MatCopy(typeRNum *C, ctypeRNum *A, ctypeInt n1, ctypeInt n2)
{
	/* matrix copy C = A *
	 * A,C: (n1 x n2)    */
	memcpy(C, A, n1*n2 * sizeof(*C));
}

void MatAdd(typeRNum *C, ctypeRNum *A, ctypeRNum *B, ctypeInt n1, ctypeInt n2)
{
	/* matrix summation C = A+B *
	 * A,B,C: (n1 x n2)         */
	typeInt i;
	for (i = 0; i < n1*n2; i++) {
		C[i] = A[i] + B[i];
	}
}

void MatMult(typeRNum *C, ctypeRNum *A, ctypeRNum *B, ctypeInt n1, ctypeInt n2, ctypeInt n3)
{
	/* matrix multiplication C = A*B			   *
	 * A: (n1 x n2)    B: (n2 x n3)   C: (n1 x n3) */
	typeInt i, j, k;
	typeRNum sigma;
	for (i = 0; i < n1; i++) {
		for (j = 0; j < n3; j++) {
			sigma = 0;
			for (k = 0; k < n2; k++) {
				sigma += A[i*n2 + k] * B[k*n3 + j];
			}
			C[i*n3 + j] = sigma;
		}
	}
}

void MatNorm(typeRNum *norm, ctypeRNum *A, ctypeInt n1, ctypeInt n2)
{
	/* A: (n1 x n2) */
	typeInt i, j;
	norm[0] = 0;
	for (i = 0; i < n1; i++) {
		for (j = 0; j < n2; j++) {
			norm[0] = norm[0] + A[i*n2 + j] * A[i*n2 + j];
		}
	}
	norm[0] = SQRT(norm[0]);
}

void MatDiffNorm(typeRNum *norm, ctypeRNum *A, ctypeRNum *B, ctypeInt n1, ctypeInt n2)
{
	/* A,B: (n1 x n2) */
	typeInt i, j;
	norm[0] = 0;
	for (i = 0; i < n1; i++) {
		for (j = 0; j < n2; j++) {
			norm[0] = norm[0] + (A[i*n2 + j] - B[i*n2 + j]) * (A[i*n2 + j] - B[i*n2 + j]);
		}
	}
	norm[0] = SQRT(norm[0]);
}

typeInt grampc_estim_penmin(typeGRAMPC *grampc, ctypeInt rungrampc) {
    typeRNum *xT_, *x_, *u_, *p_, h;
	typeRNum PenaltyMin_constr = grampc->opt->PenaltyMin;
	typeRNum PenaltyMin_tol = grampc->opt->PenaltyMin;
	typeRNum PenaltyMin = grampc->opt->PenaltyMin;
	typeRNum NormTol = 0;
	typeRNum NormConstraints = 0;
	typeInt i, j, MaxGradIter, MaxMultIter, Status;

#if defined(FIXEDSIZE) && (NC > 0)
    typeRNum Constraints[NC] = {};
#else
    typeRNum *Constraints = NULL;
#endif

	/* run grampc if flag is set to determine J and the trajectories */
	if (rungrampc) {
		/* limit the maximum number of iterations (important if e.g. the system is instable) */
		MaxGradIter = grampc->opt->MaxGradIter;
		MaxMultIter = grampc->opt->MaxMultIter;
		if (MaxGradIter > 20) {
			grampc->opt->MaxGradIter = 20;
		}
		if (MaxMultIter > 20) {
			grampc->opt->MaxMultIter = 20;
		}
		grampc_run(grampc);
	}

#ifndef FIXEDSIZE
    /* allocate memory for constraint evaluation */
    createNumMatrix(&Constraints, grampc->param->Nc);
#endif

	/* assign parameters and final states */
	if (grampc->opt->ScaleProblem == INT_ON) {
		ASSIGN_P(p_, grampc->rws->p, grampc);
		ASSIGN_X(xT_, grampc->rws->x + grampc->param->Nx*(grampc->opt->Nhor - 1), grampc);
	}
	else {
		p_ = grampc->rws->p;
		xT_ = grampc->rws->x + grampc->param->Nx*(grampc->opt->Nhor - 1);
	}

	/* evaluate and integrate the general constraints */
	for (i = 0; i < grampc->opt->Nhor; i++) {
		if (grampc->opt->ScaleProblem == INT_ON) {
			ASSIGN_X(x_, grampc->rws->x + i * grampc->param->Nx, grampc);
			ASSIGN_U(u_, grampc->rws->u + i * grampc->param->Nu, grampc);
		}
		else {
			x_ = grampc->rws->x + i * grampc->param->Nx;
			u_ = grampc->rws->u + i * grampc->param->Nu;
		}
		gfct(Constraints, grampc->rws->t[i], x_, u_, p_, grampc->userparam);
		hfct(Constraints + grampc->param->Ng, grampc->rws->t[i], x_, u_, p_, grampc->userparam);

		/* scale constraints */
		if (grampc->opt->ScaleProblem == INT_ON) {
			scale_constraints(Constraints, grampc->opt->cScale, grampc->param->Ng + grampc->param->Nh);
		}

		/* integrate the squared constraints via trapezoidal rule */
		if (i == 0) {
			h = (grampc->rws->t[i + 1] - grampc->rws->t[i]) / 2;
		}
		else if (i <= grampc->opt->Nhor - 2) {
			h = (grampc->rws->t[i + 1] - grampc->rws->t[i - 1]) / 2;
		}
		else {
			h = (grampc->rws->t[i] - grampc->rws->t[i - 1]) / 2;
		}
		for (j = 0; j < grampc->param->Ng + grampc->param->Nh; j++) {
			NormConstraints += Constraints[j] * Constraints[j] * h;
		}
	}

	/* evaluate the terminal constraints */
	gTfct(Constraints + grampc->param->Ng + grampc->param->Nh, grampc->param->Thor, xT_, p_, grampc->userparam);
	hTfct(Constraints + grampc->param->Ng + grampc->param->Nh + grampc->param->NgT, grampc->param->Thor, xT_, p_, grampc->userparam);

	/* scale constraints */
	if (grampc->opt->ScaleProblem == INT_ON) {
		scale_constraints(Constraints + grampc->param->Ng + grampc->param->Nh, grampc->opt->cScale + grampc->param->Ng + grampc->param->Nh, grampc->param->NgT + grampc->param->NhT);
	}

	/* sum the squared terminal constarints */
	for (i = grampc->param->Ng + grampc->param->Nh; i < grampc->param->Nc; i++) {
		NormConstraints += (Constraints[i] * Constraints[i]);
	}

	/* set the constraints in relation to the costs */
	if (NormConstraints > 0) {
		PenaltyMin_constr = 2 * ABS(grampc->sol->J[0]) / NormConstraints;
	}

#ifndef FIXEDSIZE
    /* free memory for constraint evaluation */
	free(Constraints);
#endif

	/* set the constraint tolerances in relation to the costs */
	for (i = 0; i < grampc->param->Nc; i++) {
		if (i < grampc->param->Ng + grampc->param->Nh) {
			NormTol = NormTol + grampc->param->Thor*grampc->opt->ConstraintsAbsTol[i] * grampc->opt->ConstraintsAbsTol[i];
		}
		else {
			NormTol = NormTol + grampc->opt->ConstraintsAbsTol[i] * grampc->opt->ConstraintsAbsTol[i];
		}
	}
	if (NormTol > 0) {
		PenaltyMin_tol = 2 * (typeRNum)1e-6 * ABS(grampc->sol->J[0]) / NormTol;
	}

	/* determine the resulting parameter and set the status */
	if (PenaltyMin_constr > PenaltyMin_tol) {
		PenaltyMin = PenaltyMin_constr;
		Status = 1;
	}
	else {
		PenaltyMin = PenaltyMin_tol;
		Status = 2;
	}
	if (PenaltyMin > grampc->opt->PenaltyMax / 500) {
		PenaltyMin = grampc->opt->PenaltyMax / 500;
		Status = 0;
	}
	grampc_setopt_real(grampc, "PenaltyMin", PenaltyMin);

	/* reset grampc if grampc_run was evaluated */
	if (rungrampc) {
		/* reset iteration limits */
		grampc->opt->MaxGradIter = MaxGradIter;
		grampc->opt->MaxMultIter = MaxMultIter;

		/* reset sol structure */
		MatSetScalar(grampc->sol->xnext, 0, grampc->param->Nx, 1);
		MatSetScalar(grampc->sol->unext, 0, grampc->param->Nu, 1);
		MatSetScalar(grampc->sol->pnext, 0, grampc->param->Np, 1);
		grampc->sol->Tnext = 0;
		MatSetScalar(grampc->sol->J, 0, 2, 1);
		grampc->sol->cfct = 0;
        grampc->sol->pen = 0;
        memset(grampc->sol->iter, 0, grampc->opt->MaxMultIter);
		grampc->sol->status = 0;

		/* reset rws structure */
		MatSetScalar(grampc->rws->tls, 0, grampc->opt->Nhor, 1);

		MatSetScalar(grampc->rws->x, 0, grampc->opt->Nhor, grampc->param->Nx);
		MatSetScalar(grampc->rws->adj, 0, grampc->opt->Nhor, grampc->param->Nx);
		MatSetScalar(grampc->rws->dcdx, 0, (grampc->opt->Nhor + 1), grampc->param->Nx);

		MatSetScalar(grampc->rws->uls, 0, grampc->opt->Nhor, grampc->param->Nu);
		MatSetScalar(grampc->rws->gradu, 0, grampc->opt->Nhor, grampc->param->Nu);
		MatSetScalar(grampc->rws->graduprev, 0, grampc->opt->Nhor, grampc->param->Nu);
		MatSetScalar(grampc->rws->dcdu, 0, grampc->opt->Nhor, grampc->param->Nu);

		MatSetScalar(grampc->rws->pls, 0, grampc->param->Np, 1);
		MatSetScalar(grampc->rws->gradp, 0, grampc->param->Np, 1);
		MatSetScalar(grampc->rws->gradpprev, 0, grampc->param->Np, 1);
		MatSetScalar(grampc->rws->dcdp, 0, (grampc->opt->Nhor + 1), grampc->param->Np);

		grampc->rws->gradT = 0;
		grampc->rws->gradTprev = 0;
		grampc->rws->dcdt = 0;

		MatSetScalar(grampc->rws->mult, 0, grampc->opt->Nhor, grampc->param->Nc);
		MatSetScalar(grampc->rws->pen, 0, grampc->opt->Nhor, grampc->param->Nc);
		MatSetScalar(grampc->rws->cfct, 0, grampc->opt->Nhor, grampc->param->Nc);
		MatSetScalar(grampc->rws->cfctprev, 0, grampc->opt->Nhor, grampc->param->Nc);

        MatSetScalar(grampc->rws->rwsScale, 0, 2 * (grampc->param->Nx + grampc->param->Nu + grampc->param->Np), 1);

		/* Initialization */
		init_rws_time(grampc);
		init_rws_controls(grampc);
		init_rws_parameters(grampc);
		init_rws_multipliers(grampc);
        init_rws_linesearch(grampc);
	}
	return Status;
}
