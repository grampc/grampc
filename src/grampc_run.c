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


#include "grampc_run.h"


void grampc_run(const typeGRAMPC *grampc)
{
	typeInt i, j, k, imult, igrad;
	typeBoolean sysintegrated = 0;
	typeRNum alpha;
	typeRNum cfct_norm;
	typeRNum pen_norm;

	typeRNum *t = grampc->rws->t;
	typeRNum *u = grampc->rws->u;
	typeRNum *gradu = grampc->rws->gradu;
	typeRNum *p = grampc->rws->p;
	typeRNum *gradp = grampc->rws->gradp;
	typeRNum *T = &grampc->rws->T;
	typeRNum *gradT = &grampc->rws->gradT;

	/* Reset solution structure */
	typeBoolean converged_grad = 0;
	typeBoolean converged_const = 0;
	for (imult = 0; imult < grampc->opt->MaxMultIter; imult++) {
		grampc->sol->iter[imult] = 0;
	}
	grampc->sol->status = STATUS_NONE;

	/* Validate necessary parameters */
	if (grampc->param->dt <= 0.0) {
		grampc_error(DT_NOT_VALID);
	}
	if (grampc->param->Thor < grampc->param->dt) {
		grampc_error(THOR_NOT_VALID);
	}

	/* Initial condition */
	if (grampc->opt->ScaleProblem == INT_ON) {
		scale_states(grampc->rws->x, grampc->param->x0, grampc);
	}
	else {
		MatCopy(grampc->rws->x, grampc->param->x0, 1, grampc->param->Nx);
	}

	/* Shift of input, Lagrange multiplier, and penalty parameter trajectories */
	/* Attention shiftTrajecotry and shortenTrajectory require a uniform grid size over the horizon */
	if (grampc->opt->ShiftControl == INT_ON) {
		if (grampc->opt->OptimTime == INT_OFF) {
			shiftTrajectory(u, grampc->opt->Nhor, grampc->param->Nu, grampc->param->Nu, grampc->param->dt, t);
			shiftTrajectory(grampc->rws->uprev, grampc->opt->Nhor, grampc->param->Nu, grampc->param->Nu, grampc->param->dt, t);
			shiftTrajectory(grampc->rws->graduprev, grampc->opt->Nhor, grampc->param->Nu, grampc->param->Nu, grampc->param->dt, t);
			shiftTrajectory(grampc->rws->cfctprev, grampc->opt->Nhor, grampc->param->Nc, (grampc->param->Ng + grampc->param->Nh), grampc->param->dt, t);
			shiftTrajectory(grampc->rws->mult, grampc->opt->Nhor, grampc->param->Nc, (grampc->param->Ng + grampc->param->Nh), grampc->param->dt, t);
			shiftTrajectory(grampc->rws->pen, grampc->opt->Nhor, grampc->param->Nc, (grampc->param->Ng + grampc->param->Nh), grampc->param->dt, t);
		}
		else {
			T[0] = T[0] < grampc->param->Tmin + grampc->param->dt ? grampc->param->Tmin : T[0] - grampc->param->dt;
			shortenTrajectory(u, grampc->opt->Nhor, grampc->param->Nu, grampc->param->Nu, grampc->param->dt, t);
			shortenTrajectory(grampc->rws->uprev, grampc->opt->Nhor, grampc->param->Nu, grampc->param->Nu, grampc->param->dt, t);
			shortenTrajectory(grampc->rws->graduprev, grampc->opt->Nhor, grampc->param->Nu, grampc->param->Nu, grampc->param->dt, t);
			shortenTrajectory(grampc->rws->cfctprev, grampc->opt->Nhor, grampc->param->Nc, (grampc->param->Ng + grampc->param->Nh), grampc->param->dt, t);
			shortenTrajectory(grampc->rws->mult, grampc->opt->Nhor, grampc->param->Nc, (grampc->param->Ng + grampc->param->Nh), grampc->param->dt, t);
			shortenTrajectory(grampc->rws->pen, grampc->opt->Nhor, grampc->param->Nc, (grampc->param->Ng + grampc->param->Nh), grampc->param->dt, t);
			discretize_time(t, T[0], grampc);
		}
	}

	/* LOOP OVER NO. OF MULTIPLIER STEPS *******************************************/
	for (imult = 0; imult < grampc->opt->MaxMultIter; imult++) {

		/* LOOP OVER NO. OF GRADIENT STEPS *********************************************/
		for (igrad = 0; igrad < grampc->opt->MaxGradIter; igrad++) {

			/* Forward integration of system and evaluation of constraints */
			if (!sysintegrated) {
				evaluate_sys(t, u, p, grampc);
				evaluate_constraints(t, u, p, 1, 0, grampc);
				sysintegrated = 1;
			}

			/* Backward integration of adjoint system */
			evaluate_adjsys(t, u, p, grampc);

			if (grampc->opt->OptimControl == INT_ON) {
				/* Gradient w.r.t. u */
				evaluate_gradu(grampc);
			}
			if (grampc->opt->OptimParam == INT_ON) {
				/* Gradient w.r.t. p */
				evaluate_gradp(grampc);
			}
			if (grampc->opt->OptimTime == INT_ON) {
				/* Gradient w.r.t. T */
				evaluate_gradT(grampc);
			}

			/* Determine step size alpha by linesearch */
			alpha = 0;
			if (grampc->opt->LineSearchType == INT_ADAPTIVELS) {
				/* Adaptive line search */
				linesearch_adaptive(&alpha, igrad, grampc);
			}
			else {
				/* Explicit line search */
				linesearch_explicit(&alpha, grampc);
			}

			if (grampc->opt->OptimControl == INT_ON) {
				/* Save previous values */
				MatCopy(grampc->rws->uprev, u, grampc->opt->Nhor, grampc->param->Nu);
				MatCopy(grampc->rws->graduprev, gradu, grampc->opt->Nhor, grampc->param->Nu);
				/* Update control */
				for (i = 0; i < grampc->opt->Nhor; i++) {
					for (j = 0; j < grampc->param->Nu; j++) {
						u[i*grampc->param->Nu + j] = u[i*grampc->param->Nu + j] - alpha * gradu[i*grampc->param->Nu + j];
					}
				}
				inputproj(u, grampc);
			}
			if (grampc->opt->OptimParam == INT_ON) {
				/* Save previous values */
				MatCopy(grampc->rws->pprev, p, 1, grampc->param->Np);
				MatCopy(grampc->rws->gradpprev, gradp, 1, grampc->param->Np);
				/* Update parameters */
				for (j = 0; j < grampc->param->Np; j++) {
					p[j] = p[j] - grampc->opt->OptimParamLineSearchFactor * alpha * gradp[j];
				}
				paramproj(p, grampc);
			}
			if (grampc->opt->OptimTime == INT_ON) {
				/* Save previous values */
				grampc->rws->Tprev = T[0];
				grampc->rws->gradTprev = gradT[0];
				/* Update time */
				T[0] = T[0] - grampc->opt->OptimTimeLineSearchFactor * alpha * gradT[0];
				timeproj(T, grampc);
				discretize_time(t, T[0], grampc);
			}

			/* Trajectories updated, integration necessary */
			sysintegrated = 0;

			/* Convergence test for gradient */
			if (grampc->opt->ConvergenceCheck == INT_ON) {
				converged_grad = convergence_test_gradient(grampc->opt->ConvergenceGradientRelTol, grampc);
				if (converged_grad) {
					grampc->sol->status |= STATUS_GRADIENT_CONVERGED;
					igrad++;
					break;
				}
			}
		}
		/* END GRADIENT LOOP ***********************************************************/
		grampc->sol->iter[imult] = igrad;

		/* Forward integration of system and evaluation of constraints */
		evaluate_sys(t, u, p, grampc);
		evaluate_constraints(t, u, p, (imult + 1 < grampc->opt->MaxMultIter), 1, grampc);
		sysintegrated = 1;

		/* Convergence test for constraints */
		if (grampc->opt->ConvergenceCheck == INT_ON && converged_grad) {
			converged_const = convergence_test_constraints(grampc->rws->cfctAbsTol, grampc);
			if (converged_const) {
				grampc->sol->status |= STATUS_CONSTRAINTS_CONVERGED;
				break;
			}
		}
	}
	/* END MULTIPLIER LOOP ***********************************************************/

	/* Calculation of xnext */
	interplin(grampc->sol->xnext, t, grampc->rws->x, grampc->param->dt, grampc->param->Nx, grampc->opt->Nhor, 1);

	/* Unscaling */
	if (grampc->opt->ScaleProblem == INT_ON) {
		unscale_states(grampc->sol->xnext, grampc->sol->xnext, grampc);
		unscale_controls(grampc->sol->unext, u, grampc);
		unscale_parameters(grampc->sol->pnext, p, grampc);
		unscale_time(&grampc->sol->Tnext, T[0], grampc);
	}
	else {
		MatCopy(grampc->sol->unext, u, 1, grampc->param->Nu);
		MatCopy(grampc->sol->pnext, p, 1, grampc->param->Np);
		grampc->sol->Tnext = T[0];
	}

	/* Calculation of cost J */
	evaluate_cost(grampc->sol->J, t, u, p, grampc);
	/* Calculation of constraint norm and penalty norm*/
	cfct_norm = 0;
	for (i = 0; i < grampc->opt->Nhor; i++) {
		k = i * grampc->param->Nc;
		for (j = 0; j < grampc->param->Ng; j++) {
			cfct_norm = cfct_norm + grampc->rws->cfct[k + j] * grampc->rws->cfct[k + j];
		}
		for (; j < grampc->param->Ng + grampc->param->Nh; j++) {
			cfct_norm = cfct_norm + MAX(grampc->rws->cfct[k + j], 0) * MAX(grampc->rws->cfct[k + j], 0);
		}
		for (; j < grampc->param->Ng + grampc->param->Nh + grampc->param->NgT; j++) {
			cfct_norm = cfct_norm + grampc->rws->cfct[k + j] * grampc->rws->cfct[k + j];
		}
		for (; j < grampc->param->Ng + grampc->param->Nh + grampc->param->NgT + grampc->param->NhT; j++) {
			cfct_norm = cfct_norm + MAX(grampc->rws->cfct[k + j], 0) * MAX(grampc->rws->cfct[k + j], 0);
		}
	}
	cfct_norm = SQRT(cfct_norm);
	MatNorm(&pen_norm, grampc->rws->pen, grampc->opt->Nhor, grampc->param->Nc);
	if (cfct_norm > grampc->sol->cfct && pen_norm >= grampc->sol->pen && !converged_const) {
		grampc->sol->status |= STATUS_INFEASIBLE;
	}
	grampc->sol->cfct = cfct_norm;
	grampc->sol->pen = pen_norm;
}


void evaluate_constraints(ctypeRNum *t, ctypeRNum *u, ctypeRNum *p, const typeBoolean evaljac, const typeBoolean updatemultiplier, const typeGRAMPC *grampc)
{
	typeInt i;
    typeInt N;

	ctypeRNum *x_ = NULL;
	ctypeRNum *u_ = NULL;
	ctypeRNum *p_ = p;
	typeRNum *mult = NULL;
	typeRNum *pen = NULL;
	typeRNum *cfct = NULL;
	typeRNum *cfctprev = NULL;
	typeRNum *dcdx = NULL;
	typeRNum *dcdu = NULL;
	typeRNum *dcdp = NULL;
	typeRNum *dcdt = &grampc->rws->dcdt;
	typeRNum *thresholds = NULL;
	typeRNum *cScale = NULL;
	typeBoolean converged_grad = 0;

	typeRNum *c = grampc->rws->rwsGeneral; /* size:  Nc+ 2*(Nx+Nu+Np) */
	typeRNum *dgdxvec = c + grampc->param->Nc;
	typeRNum *dhdxvec = dgdxvec + grampc->param->Nx;
	typeRNum *dgdpvec = dhdxvec + grampc->param->Nx;
	typeRNum *dhdpvec = dgdpvec + grampc->param->Np;
	typeRNum *dgduvec = dhdpvec + grampc->param->Np;
	typeRNum *dhduvec = dgduvec + grampc->param->Nu;
	typeRNum dgTdT, dhTdT;

	/* return if no constraints are defined */
	if (grampc->param->Nc == 0) {
		return;
	}

	/* Init constraint jacobians */
	MatSetScalar(grampc->rws->dcdx, 0, grampc->opt->Nhor + 1, grampc->param->Nx);
	MatSetScalar(grampc->rws->dcdu, 0, grampc->opt->Nhor, grampc->param->Nu);
	MatSetScalar(grampc->rws->dcdp, 0, grampc->opt->Nhor + 1, grampc->param->Np);
	MatSetScalar(dgdxvec, 0, 1, 2 * (grampc->param->Nx + grampc->param->Np + grampc->param->Nu));

	/* Multiplier update depends on the convergence of the subproblem */
	if (updatemultiplier) {
		converged_grad = convergence_test_gradient(grampc->opt->AugLagUpdateGradientRelTol, grampc);
		if (converged_grad) {
			grampc->sol->status |= STATUS_MULTIPLIER_UPDATE;
		}
	}

	/* Unscaling */
	if (grampc->opt->ScaleProblem == INT_ON) {
		ASSIGN_P(p_, p, grampc);
	}

	/* loop over the prediction horizon */
	if ((grampc->param->Ng + grampc->param->Nh > 0) && ((grampc->opt->EqualityConstraints == INT_ON) || (grampc->opt->InequalityConstraints == INT_ON))) {
        /* if there are terminal constraints, the integral constraints are not evaluated for the last point */
        if (grampc->param->NgT + grampc->param->NhT > 0) {
            N = grampc->opt->Nhor - 1;
        }
        /* if there are no terminal constraints, the integral constraints are evaluated for all points */
        else {
            N = grampc->opt->Nhor;
        }
        for (i = 0; i < N; i++)
		{
			mult = grampc->rws->mult + i * grampc->param->Nc;
			pen = grampc->rws->pen + i * grampc->param->Nc;
			cfct = grampc->rws->cfct + i * grampc->param->Nc;
			cfctprev = grampc->rws->cfctprev + i * grampc->param->Nc;
			dcdx = grampc->rws->dcdx + i * grampc->param->Nx;
			dcdu = grampc->rws->dcdu + i * grampc->param->Nu;
			dcdp = grampc->rws->dcdp + i * grampc->param->Np;
			thresholds = grampc->rws->cfctAbsTol;
			cScale = grampc->opt->cScale;

			/* Unscaling */
			if (grampc->opt->ScaleProblem == INT_ON) {
				ASSIGN_X(x_, grampc->rws->x + i * grampc->param->Nx, grampc);
				ASSIGN_U(u_, u + i * grampc->param->Nu, grampc);
			}
			else {
				x_ = grampc->rws->x + i * grampc->param->Nx;
				u_ = u + i * grampc->param->Nu;
			}

			/* EQUALITY CONSTRAINTS **************************************************/
			if (grampc->param->Ng > 0 && grampc->opt->EqualityConstraints == INT_ON) {

				/* Evaluate constraints */
				gfct(cfct, t[i], x_, u_, p_, grampc->userparam);
				if (grampc->opt->ScaleProblem == INT_ON) {
					scale_constraints(cfct, cScale, grampc->param->Ng);
				}

				/* Update multipliers */
				if (updatemultiplier) {
					update_multiplier_eqc(mult, pen, cfct, cfctprev, thresholds, grampc->param->Ng, converged_grad, grampc);
				}

				/* Evaluate jacobians */
				if (evaljac) {
					compute_jacobian_multiplier(c, mult, pen, cfct, grampc->param->Ng);
					if (grampc->opt->ScaleProblem == INT_ON) {
						scale_constraints(c, cScale, grampc->param->Ng);
					}

					dgdx_vec(dgdxvec, t[i], x_, u_, p_, c, grampc->userparam);
					MatAdd(dcdx, dcdx, dgdxvec, 1, grampc->param->Nx);

					if (grampc->opt->OptimControl == INT_ON) {
						dgdu_vec(dgduvec, t[i], x_, u_, p_, c, grampc->userparam);
						MatAdd(dcdu, dcdu, dgduvec, 1, grampc->param->Nu);
					}

					if (grampc->opt->OptimParam == INT_ON) {
						dgdp_vec(dgdpvec, t[i], x_, u_, p_, c, grampc->userparam);
						MatAdd(dcdp, dcdp, dgdpvec, 1, grampc->param->Np);
					}
				}
			}

			mult = mult + grampc->param->Ng;
			pen = pen + grampc->param->Ng;
			cfct = cfct + grampc->param->Ng;
			cfctprev = cfctprev + grampc->param->Ng;
			thresholds = thresholds + grampc->param->Ng;
			cScale = cScale + grampc->param->Ng;

			/* INEQUALITY CONSTRAINTS ************************************************/
			if (grampc->param->Nh > 0 && grampc->opt->InequalityConstraints == INT_ON) {

				/* Evaluate constraints */
				hfct(cfct, t[i], x_, u_, p_, grampc->userparam);
				if (grampc->opt->ScaleProblem == INT_ON) {
					scale_constraints(cfct, cScale, grampc->param->Nh);
				}
				update_cfct_for_ieqc(mult, pen, cfct, grampc->param->Nh);

				/* Update multipliers */
				if (updatemultiplier) {
					update_multiplier_ieqc(mult, pen, cfct, cfctprev, thresholds, grampc->param->Nh, converged_grad, grampc);
				}

				/* Evaluate jacobians */
				if (evaljac) {
					compute_jacobian_multiplier(c, mult, pen, cfct, grampc->param->Nh);
					if (grampc->opt->ScaleProblem == INT_ON) {
						scale_constraints(c, cScale, grampc->param->Nh);
					}

					dhdx_vec(dhdxvec, t[i], x_, u_, p_, c, grampc->userparam);
					MatAdd(dcdx, dcdx, dhdxvec, 1, grampc->param->Nx);

					if (grampc->opt->OptimControl == INT_ON) {
						dhdu_vec(dhduvec, t[i], x_, u_, p_, c, grampc->userparam);
						MatAdd(dcdu, dcdu, dhduvec, 1, grampc->param->Nu);
					}

					if (grampc->opt->OptimParam == INT_ON) {
						dhdp_vec(dhdpvec, t[i], x_, u_, p_, c, grampc->userparam);
						MatAdd(dcdp, dcdp, dhdpvec, 1, grampc->param->Np);
					}
				}
			}
		}
	}

	/* Terminal Constraints */
	if ((grampc->param->NgT + grampc->param->NhT > 0) && ((grampc->opt->TerminalEqualityConstraints == INT_ON) || (grampc->opt->TerminalInequalityConstraints == INT_ON)) ){

		i = grampc->opt->Nhor - 1;
		mult = grampc->rws->mult + i * grampc->param->Nc + grampc->param->Ng + grampc->param->Nh;
		pen = grampc->rws->pen + i * grampc->param->Nc + grampc->param->Ng + grampc->param->Nh;
		cfct = grampc->rws->cfct + i * grampc->param->Nc + grampc->param->Ng + grampc->param->Nh;
		cfctprev = grampc->rws->cfctprev + i * grampc->param->Nc + grampc->param->Ng + grampc->param->Nh;
		/* dcdx for Terminal constraint must be different from dcdx for integral constraints */
		dcdx = grampc->rws->dcdx + grampc->opt->Nhor * grampc->param->Nx;
		/* dcdp for Terminal constraint must be different from dcdp for integral constraints */
		dcdp = grampc->rws->dcdp + grampc->opt->Nhor * grampc->param->Np;
		/* dcdt for Terminal constraint must be different from dcdt for integral constraints */
		*dcdt = 0;
		thresholds = grampc->rws->cfctAbsTol + grampc->param->Ng + grampc->param->Nh;
		cScale = grampc->opt->cScale + grampc->param->Ng + grampc->param->Nh;

		/* Init constraint jacobians */;
		MatSetScalar(dgdxvec, 0, 1, 2 * (grampc->param->Nx + grampc->param->Np));
		dgTdT = 0;
		dhTdT = 0;

		/* Unscaling */
		if (grampc->opt->ScaleProblem == INT_ON) {
			ASSIGN_X(x_, grampc->rws->x + i * grampc->param->Nx, grampc);
		}
		else {
			x_ = grampc->rws->x + i * grampc->param->Nx;
		}

		/* TERMINAL EQUALITY CONSTRAINTS *******************************************/
		if (grampc->param->NgT > 0 && grampc->opt->TerminalEqualityConstraints == INT_ON) {

			/* Evaluate constraints */
			gTfct(cfct, t[i], x_, p_, grampc->userparam);
			if (grampc->opt->ScaleProblem == INT_ON) {
				scale_constraints(cfct, cScale, grampc->param->NgT);
			}

			/* Update multipliers */
			if (updatemultiplier) {
				update_multiplier_eqc(mult, pen, cfct, cfctprev, thresholds, grampc->param->NgT, converged_grad, grampc);
			}

			/* Evaluate jacobians */
			if (evaljac) {
				compute_jacobian_multiplier(c, mult, pen, cfct, grampc->param->NgT);
				if (grampc->opt->ScaleProblem == INT_ON) {
					scale_constraints(c, cScale, grampc->param->NgT);
				}

				dgTdx_vec(dgdxvec, t[i], x_, p_, c, grampc->userparam);
				MatAdd(dcdx, dcdx, dgdxvec, 1, grampc->param->Nx);

				if (grampc->opt->OptimParam == INT_ON) {
					dgTdp_vec(dgdpvec, t[i], x_, p_, c, grampc->userparam);
					MatAdd(dcdp, dcdp, dgdpvec, 1, grampc->param->Np);
				}

				if (grampc->opt->OptimTime == INT_ON) {
					dgTdT_vec(&dgTdT, t[i], x_, p_, c, grampc->userparam);
					*dcdt = *dcdt + dgTdT;
				}
			}
		}

		mult = mult + grampc->param->NgT;
		pen = pen + grampc->param->NgT;
		cfct = cfct + grampc->param->NgT;
		cfctprev = cfctprev + grampc->param->NgT;
		thresholds = thresholds + grampc->param->NgT;
		cScale = cScale + grampc->param->NgT;

		/* TERMINAL INEQUALITY CONSTRAINTS *****************************************/
		if (grampc->param->NhT > 0 && grampc->opt->TerminalInequalityConstraints == INT_ON) {

			/* Evaluate constraints */
			hTfct(cfct, t[i], x_, p_, grampc->userparam);
			if (grampc->opt->ScaleProblem == INT_ON) {
				scale_constraints(cfct, cScale, grampc->param->NhT);
			}
			update_cfct_for_ieqc(mult, pen, cfct, grampc->param->NhT);

			/* Update multipliers */
			if (updatemultiplier) {
				update_multiplier_ieqc(mult, pen, cfct, cfctprev, thresholds, grampc->param->NhT, converged_grad, grampc);
			}

			/* Evaluate jacobians */
			if (evaljac) {
				compute_jacobian_multiplier(c, mult, pen, cfct, grampc->param->NhT);
				if (grampc->opt->ScaleProblem == INT_ON) {
					scale_constraints(c, cScale, grampc->param->NhT);
				}

				dhTdx_vec(dhdxvec, t[i], x_, p_, c, grampc->userparam);
				MatAdd(dcdx, dcdx, dhdxvec, 1, grampc->param->Nx);

				if (grampc->opt->OptimParam == INT_ON) {
					dhTdp_vec(dhdpvec, t[i], x_, p_, c, grampc->userparam);
					MatAdd(dcdp, dcdp, dhdpvec, 1, grampc->param->Np);
				}

				if (grampc->opt->OptimTime == INT_ON) {
					dhTdT_vec(&dhTdT, t[i], x_, p_, c, grampc->userparam);
					*dcdt = *dcdt + dhTdT;
				}
			}
		}
	}
}

void update_multiplier_eqc(typeRNum *mult, typeRNum *pen, ctypeRNum *cfct, typeRNum *cfctprev,
	ctypeRNum *thresholds, ctypeInt Ncon, typeBoolean converged_grad, const typeGRAMPC *grampc)
{
	typeInt i;
	for (i = 0; i < Ncon; i++) {
		/* Increase multipliers for violated constraints if minimization converged */
		if (ABS(cfct[i]) > thresholds[i] && converged_grad) {
			if (grampc->opt->ConstraintsHandling == INT_AUGLAG) {
				mult[i] = mult[i] + (1 - grampc->opt->MultiplierDampingFactor) * pen[i] * cfct[i];
				if (mult[i] > grampc->opt->MultiplierMax) {
					mult[i] = grampc->opt->MultiplierMax;
					grampc->sol->status |= STATUS_MULTIPLIER_MAX;
				}
				else if (mult[i] < -grampc->opt->MultiplierMax) {
					mult[i] = -grampc->opt->MultiplierMax;
					grampc->sol->status |= STATUS_MULTIPLIER_MAX;
				}
			}
			if (ABS(cfct[i]) > grampc->opt->PenaltyIncreaseThreshold * ABS(cfctprev[i])) {
				pen[i] = pen[i] * grampc->opt->PenaltyIncreaseFactor;
				if (pen[i] > grampc->opt->PenaltyMax) {
					pen[i] = grampc->opt->PenaltyMax;
					grampc->sol->status |= STATUS_PENALTY_MAX;
				}
			}
			/* Save cfct as cfctprev for next multiplier update */
			cfctprev[i] = cfct[i];
		}
		/* Decrease multipliers for satisfied constraints */
		if (ABS(cfct[i]) < thresholds[i] / 10.0) {
			pen[i] = MAX(pen[i] * grampc->opt->PenaltyDecreaseFactor, grampc->opt->PenaltyMin);
		}
	}
}

void update_multiplier_ieqc(typeRNum *mult, typeRNum *pen, ctypeRNum *cfct, typeRNum *cfctprev,
	ctypeRNum *thresholds, ctypeInt Ncon, typeBoolean converged_grad, const typeGRAMPC *grampc)
{
	typeInt i;
	for (i = 0; i < Ncon; i++) {
		/* Increase multipliers for violated constraints if minimization converged */
		if (cfct[i] > thresholds[i] && converged_grad) {
			if (grampc->opt->ConstraintsHandling == INT_AUGLAG) {
				mult[i] = mult[i] + (1 - grampc->opt->MultiplierDampingFactor) * pen[i] * cfct[i];
				if (mult[i] > grampc->opt->MultiplierMax) {
					mult[i] = grampc->opt->MultiplierMax;
					grampc->sol->status |= STATUS_MULTIPLIER_MAX;
				}
			}
			if (cfct[i] > grampc->opt->PenaltyIncreaseThreshold * cfctprev[i]) {
				pen[i] = pen[i] * grampc->opt->PenaltyIncreaseFactor;
				if (pen[i] > grampc->opt->PenaltyMax) {
					pen[i] = grampc->opt->PenaltyMax;
					grampc->sol->status |= STATUS_PENALTY_MAX;
				}
			}
			/* Save cfct as cfctprev for next multiplier update */
			cfctprev[i] = cfct[i];
		}
		/* Decrease multipliers for satisfied constraints */
		if (cfct[i] < thresholds[i] / 10.0) {
			if (grampc->opt->ConstraintsHandling == INT_AUGLAG && cfct[i] < 0) {
				mult[i] = mult[i] + (1 - grampc->opt->MultiplierDampingFactor) * pen[i] * cfct[i];
			}
			pen[i] = MAX(pen[i] * grampc->opt->PenaltyDecreaseFactor, grampc->opt->PenaltyMin);
		}
	}
}

void update_cfct_for_ieqc(ctypeRNum *mult, ctypeRNum *pen, typeRNum *cfct, ctypeInt Ncon)
{
	typeInt i;
	for (i = 0; i < Ncon; i++) {
		cfct[i] = MAX(cfct[i], -mult[i] / pen[i]);
	}
}

void compute_jacobian_multiplier(typeRNum *c, ctypeRNum *mult, ctypeRNum *pen, ctypeRNum *cfct, ctypeInt Ncon)
{
	typeInt i;
	for (i = 0; i < Ncon; i++) {
		c[i] = mult[i] + pen[i] * cfct[i];
	}
}


void evaluate_sys(ctypeRNum *t, ctypeRNum *u, ctypeRNum *p, const typeGRAMPC *grampc)
{
	typeIntffctPtr pIntSys;
	ctypeRNum *p_ = p;

	/* integrator for system integration */
	if (grampc->opt->Integrator == INT_EULER) {
		pIntSys = &intsysEuler;
	}
	else if (grampc->opt->Integrator == INT_MODEULER) {
		pIntSys = &intsysModEuler;
	}
	else if (grampc->opt->Integrator == INT_HEUN) {
		pIntSys = &intsysHeun;
	}
	else if (grampc->opt->Integrator == INT_RODAS) {
		pIntSys = &intsysRodas;
	}
	else {
		pIntSys = &intsysRuKu45;
	}

	/* Unscaling */
	if (grampc->opt->ScaleProblem == INT_ON) {
		ASSIGN_P(p_, p, grampc);
	}

	(*pIntSys)(grampc->rws->x, FWINT, grampc->opt->Nhor, t, grampc->rws->x, u, p_, grampc, &Wsys);
}

void evaluate_adjsys(ctypeRNum *t, ctypeRNum *u, ctypeRNum *p, const typeGRAMPC *grampc)
{
	typeInt i = grampc->opt->Nhor - 1;
	ctypeRNum *x_ = grampc->rws->x + i * grampc->param->Nx;
	ctypeRNum *p_ = p;
	typeRNum *adj_ = grampc->rws->adj + i * grampc->param->Nx;

	typeIntffctPtr pIntSys;

	/* integrator for adjoint system integration */
	if (grampc->opt->Integrator == INT_EULER) {
		pIntSys = &intsysEuler;
	}
	else if (grampc->opt->Integrator == INT_MODEULER) {
		pIntSys = &intsysModEuler;
	}
	else if (grampc->opt->Integrator == INT_HEUN) {
		pIntSys = &intsysHeun;
	}
	else if (grampc->opt->Integrator == INT_RODAS) {
		pIntSys = &intsysRodas;
	}
	else {
		pIntSys = &intsysRuKu45;
	}

	/* Unscaling */
	if (grampc->opt->ScaleProblem == INT_ON) {
		ASSIGN_X(x_, grampc->rws->x + i * grampc->param->Nx, grampc);
		ASSIGN_P(p_, p, grampc);
	}

	/* Terminal condition for adjoint states */
	MatSetScalar(adj_, 0, 1, grampc->param->Nx);

	/* Jacobian of Terminal cost */
	if (grampc->opt->TerminalCost == INT_ON) {
		dVdx(adj_, t[i], x_, p_, grampc->param->xdes, grampc->userparam);
		if (grampc->opt->ScaleProblem == INT_ON) {
			scale_cost(adj_, grampc->opt->JScale, grampc->param->Nx);
		}
	}

	/* Jacobian of Terminal equality and inequality constraints */
	if ((grampc->param->NgT + grampc->param->NhT > 0) && (grampc->opt->TerminalEqualityConstraints == INT_ON
		|| grampc->opt->TerminalInequalityConstraints == INT_ON)) {
		MatAdd(adj_, adj_, grampc->rws->dcdx + grampc->opt->Nhor * grampc->param->Nx, 1, grampc->param->Nx);
	}

	/* Scaling */
	if (grampc->opt->ScaleProblem == INT_ON) {
		scale_adjoints(adj_, adj_, grampc);
	}

	/* integration of adjoint system in reverse time */
	(*pIntSys)(adj_, BWINT, grampc->opt->Nhor, t + i, grampc->rws->x + i * grampc->param->Nx, u + i * grampc->param->Nu, p_, grampc, &Wadjsys);
}

void Wsys(typeRNum *s, ctypeRNum *x, ctypeRNum *t, ctypeRNum *dummy,
	ctypeRNum *u, ctypeRNum *p_, ctypeRNum *dcdx, const typeGRAMPC *grampc)
{
	ctypeRNum *x_ = x;
	ctypeRNum *u_ = u;
	typeInt i;

	/* remove warning unreferenced formal parameter */
	(void)(dummy);
	(void)(dcdx);

	/* Unscaling */
	if (grampc->opt->ScaleProblem == INT_ON) {
		ASSIGN_X(x_, x, grampc);
		ASSIGN_U(u_, u, grampc);
	}

	ffct(s, t[0] + grampc->param->t0, x_, u_, p_, grampc->userparam);

	/* Scaling */
	if (grampc->opt->ScaleProblem == INT_ON) {
		for (i = 0; i < grampc->param->Nx; i++) {
			s[i] = s[i] / grampc->opt->xScale[i];
		}
	}
}

void Wadjsys(typeRNum *s, ctypeRNum *adj, ctypeRNum *t, ctypeRNum *x,
	ctypeRNum *u, ctypeRNum *p_, ctypeRNum *dcdx, const typeGRAMPC *grampc)
{
	ctypeRNum *x_ = x;
	ctypeRNum *adj_ = adj;
	ctypeRNum *u_ = u;
	typeRNum *dLdx = grampc->rws->rwsGeneral; /* size:  Nx */

	typeInt i;

	/* Unscaling */
	if (grampc->opt->ScaleProblem == INT_ON) {
		ASSIGN_X(x_, x, grampc);
		ASSIGN_ADJ(adj_, adj, grampc);
		ASSIGN_U(u_, u, grampc);
	}

	MatSetScalar(dLdx, 0, 1, grampc->param->Nx);

	if (grampc->opt->IntegralCost == INT_ON) {
		dldx(dLdx, t[0], x_, u_, p_, grampc->param->xdes, grampc->param->udes, grampc->userparam);
		if (grampc->opt->ScaleProblem == INT_ON) {
			scale_cost(dLdx, grampc->opt->JScale, grampc->param->Nx);
		}
	}

	dfdx_vec(s, t[0] + grampc->param->t0, x_, adj_, u_, p_, grampc->userparam);

	/* Scaling */
	if (grampc->opt->ScaleProblem == INT_ON) {
		for (i = 0; i < grampc->param->Nx; i++) {
			s[i] = (-dLdx[i] - s[i] - dcdx[i]) * grampc->opt->xScale[i];
		}
	}
	else {
		for (i = 0; i < grampc->param->Nx; i++) {
			s[i] = -dLdx[i] - s[i] - dcdx[i];
		}
	}
}

void evaluate_gradu(const typeGRAMPC *grampc)
{
	typeInt i, j;

	ctypeRNum *t = grampc->rws->t;
	ctypeRNum *x_ = NULL;
	ctypeRNum *u_ = NULL;
	ctypeRNum *p_ = grampc->rws->p;
	ctypeRNum *adj_ = NULL;
	ctypeRNum *dcdu = NULL;

	typeRNum *dLdu = grampc->rws->rwsGeneral;
	typeRNum *s = dLdu + grampc->param->Nu;
	/*                      dLdu  s  */
	/* sizeof(rwsGeneral) = Nu + Nu */

	MatSetScalar(dLdu, 0, 1, grampc->param->Nu);

	/* Unscaling */
	if (grampc->opt->ScaleProblem == INT_ON) {
		ASSIGN_P(p_, grampc->rws->p, grampc);
	}

	for (i = 0; i < grampc->opt->Nhor; i++) {

		dcdu = grampc->rws->dcdu + i * grampc->param->Nu;

		/* Unscaling */
		if (grampc->opt->ScaleProblem == INT_ON) {
			ASSIGN_X(x_, grampc->rws->x + i * grampc->param->Nx, grampc);
			ASSIGN_ADJ(adj_, grampc->rws->adj + i * grampc->param->Nx, grampc);
			ASSIGN_U(u_, grampc->rws->u + i * grampc->param->Nu, grampc);
		}
		else {
			x_ = grampc->rws->x + i * grampc->param->Nx;
			adj_ = grampc->rws->adj + i * grampc->param->Nx;
			u_ = grampc->rws->u + i * grampc->param->Nu;
		}

		if (grampc->opt->IntegralCost == INT_ON) {
			dldu(dLdu, t[i], x_, u_, p_, grampc->param->xdes, grampc->param->udes, grampc->userparam);
			if (grampc->opt->ScaleProblem == INT_ON) {
				scale_cost(dLdu, grampc->opt->JScale, grampc->param->Nu);
			}
		}

		dfdu_vec(s, t[i] + grampc->param->t0, x_, adj_, u_, p_, grampc->userparam);

		/* Scaling */
		if (grampc->opt->ScaleProblem == INT_ON) {
			for (j = 0; j < grampc->param->Nu; j++) {
				grampc->rws->gradu[i*grampc->param->Nu + j] = (dLdu[j] + s[j] + dcdu[j]) * grampc->opt->uScale[j];
			}
		}
		else {
			for (j = 0; j < grampc->param->Nu; j++) {
				grampc->rws->gradu[i*grampc->param->Nu + j] = dLdu[j] + s[j] + dcdu[j];
			}
		}
	}
}

void inputproj(typeRNum *u, const typeGRAMPC *grampc)
{
	typeInt i, j;
	typeRNum *umin_ = grampc->param->umin;
	typeRNum *umax_ = grampc->param->umax;

	if (grampc->opt->ScaleProblem == INT_ON) {
		umin_ = grampc->rws->rwsScale + 2 * grampc->param->Nx;
		umax_ = umin_ + grampc->param->Nu;
		scale_controls(umin_, grampc->param->umin, grampc);
		scale_controls(umax_, grampc->param->umax, grampc);
	}

	for (i = 0; i < grampc->opt->Nhor; i++) {
		for (j = 0; j < grampc->param->Nu; j++) {
			/* lower bound */
			if (u[j + i * grampc->param->Nu] < umin_[j]) {
				u[j + i * grampc->param->Nu] = umin_[j];
			}
			/* upper bound */
			else if (u[j + i * grampc->param->Nu] > umax_[j]) {
				u[j + i * grampc->param->Nu] = umax_[j];
			}
		}
	}
}


void evaluate_gradp(const typeGRAMPC *grampc)
{
	typeInt i, j;
	typeRNum h;

	ctypeRNum *t = grampc->rws->t;
	typeRNum *x = NULL;
	typeRNum *adj = NULL;
	typeRNum *u = NULL;
	typeRNum *p_ = grampc->rws->p;
	ctypeRNum *dcdp = NULL;

	typeRNum *gradp = grampc->rws->gradp;
	typeRNum *s = grampc->rws->rwsGeneral;

	/* Initialize to zero */
	MatSetScalar(gradp, 0.0, 1, grampc->param->Np);

	/* Unscaling */
	if (grampc->opt->ScaleProblem == INT_ON) {
		ASSIGN_P(p_, grampc->rws->p, grampc);
	}

	/* Compute integral over gradp using trapezoidal rule */
	for (i = 0; i < grampc->opt->Nhor; i++)
	{
		/* Integration step size */
		if (i == 0) {
			h = (t[i + 1] - t[i]) / 2;
		}
		else if (i <= grampc->opt->Nhor - 2) {
			h = (t[i + 1] - t[i - 1]) / 2;
		}
		else {
			h = (t[i] - t[i - 1]) / 2;
		}

		x = grampc->rws->x + i * grampc->param->Nx;
		adj = grampc->rws->adj + i * grampc->param->Nx;
		u = grampc->rws->u + i * grampc->param->Nu;
		dcdp = grampc->rws->dcdp + i * grampc->param->Np;

		WintParam(s, t[i], x, adj, u, p_, dcdp, grampc);
		for (j = 0; j < grampc->param->Np; j++) {
			gradp[j] = gradp[j] + h * s[j];
		}
	}

	x = grampc->rws->x + (grampc->opt->Nhor - 1) * grampc->param->Nx;
	dcdp = grampc->rws->dcdp + grampc->opt->Nhor * grampc->param->Np;

	WtermParam(s, t[i], x, p_, dcdp, grampc);
	MatAdd(gradp, gradp, s, 1, grampc->param->Np);

	/* Scaling */
	if (grampc->opt->ScaleProblem == INT_ON) {
		for (j = 0; j < grampc->param->Np; j++) {
			gradp[j] = gradp[j] * grampc->opt->pScale[j];
		}
	}
}

void WintParam(typeRNum *s, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj,
	ctypeRNum *u, ctypeRNum *p_, ctypeRNum *dcdp, const typeGRAMPC *grampc)
{
	typeInt i;
	ctypeRNum *x_ = x;
	ctypeRNum *adj_ = adj;
	ctypeRNum *u_ = u;

	typeRNum *dldp_val = grampc->rws->rwsGeneral + grampc->param->Np;
	typeRNum *dfdp_val = dldp_val + grampc->param->Np;
	MatSetScalar(dldp_val, 0, 1, grampc->param->Np);

	/* Unscaling */
	if (grampc->opt->ScaleProblem == INT_ON) {
		ASSIGN_X(x_, x, grampc);
		ASSIGN_ADJ(adj_, adj, grampc);
		ASSIGN_U(u_, u, grampc);
	}

	/* dl(x,u,p)/dp + df(x,u,p)/dp' * adj + dc(x,u,p)/dp' * mult */
	if (grampc->opt->IntegralCost == INT_ON) {
		dldp(dldp_val, t, x_, u_, p_, grampc->param->xdes, grampc->param->udes, grampc->userparam);
		if (grampc->opt->ScaleProblem == INT_ON) {
			scale_cost(dldp_val, grampc->opt->JScale, grampc->param->Np);
		}
	}

	dfdp_vec(dfdp_val, t + grampc->param->t0, x_, adj_, u_, p_, grampc->userparam);

	for (i = 0; i < grampc->param->Np; i++) {
		s[i] = dldp_val[i] + dfdp_val[i] + dcdp[i];
	}
}

void WtermParam(typeRNum *s, ctypeRNum t, ctypeRNum *x, ctypeRNum *p_, ctypeRNum *dcdp, const typeGRAMPC *grampc)
{
	typeInt i;
	ctypeRNum *x_ = x;

	typeRNum *dldp = grampc->rws->rwsGeneral + grampc->param->Np;
	MatSetScalar(dldp, 0, 1, grampc->param->Np);

	/* Unscaling */
	if (grampc->opt->ScaleProblem == INT_ON) {
		ASSIGN_X(x_, x, grampc);
	}

	/* dV(x(T),p)/dp + dc(x(T),p)/dp' * mult(T) */
	if (grampc->opt->TerminalCost == INT_ON) {
		dVdp(dldp, t, x_, p_, grampc->param->xdes, grampc->userparam);
		if (grampc->opt->ScaleProblem == INT_ON) {
			scale_cost(dldp, grampc->opt->JScale, grampc->param->Np);
		}
	}

	for (i = 0; i < grampc->param->Np; i++) {
		s[i] = dldp[i] + dcdp[i];
	}
}

void paramproj(typeRNum *p, const typeGRAMPC *grampc)
{
	typeInt i;
	typeRNum *pmin_ = grampc->param->pmin;
	typeRNum *pmax_ = grampc->param->pmax;

	if (grampc->opt->ScaleProblem == INT_ON) {
		pmin_ = grampc->rws->rwsScale;
		pmax_ = pmin_ + grampc->param->Np;
		scale_parameters(pmin_, grampc->param->pmin, grampc);
		scale_parameters(pmax_, grampc->param->pmax, grampc);
	}

	for (i = 0; i < grampc->param->Np; i++) {
		/* lower bound */
		if (p[i] < pmin_[i]) {
			p[i] = pmin_[i];
		}
		/* upper bound */
		else if (p[i] > pmax_[i]) {
			p[i] = pmax_[i];
		}
	}
}

void evaluate_gradT(const typeGRAMPC *grampc)
{
	typeInt i, k;

	ctypeRNum *t = grampc->rws->t;
	ctypeRNum *x_ = grampc->rws->x + (grampc->opt->Nhor - 1) * grampc->param->Nx;
	ctypeRNum *adj_ = grampc->rws->adj + (grampc->opt->Nhor - 1) * grampc->param->Nx;
	ctypeRNum *u_ = grampc->rws->u + (grampc->opt->Nhor - 1) * grampc->param->Nu;
	typeRNum *p_ = grampc->rws->p;

	ctypeRNum *mult = grampc->rws->mult + (grampc->opt->Nhor - 1) * grampc->param->Nc;
	ctypeRNum *pen = grampc->rws->pen + (grampc->opt->Nhor - 1) * grampc->param->Nc;
	ctypeRNum *cfct = grampc->rws->cfct + (grampc->opt->Nhor - 1) * grampc->param->Nc;

	typeRNum *s1 = grampc->rws->rwsGeneral;
	typeRNum l = 0.0;
	typeRNum f = 0.0;
	typeRNum c = 0.0;
	typeRNum dVdt = 0.0;

	/* Unscaling */
	if (grampc->opt->ScaleProblem == INT_ON) {
		ASSIGN_X(x_, x_, grampc);
		ASSIGN_ADJ(adj_, adj_, grampc);
		ASSIGN_U(u_, u_, grampc);
		ASSIGN_P(p_, p_, grampc);
	}

	/* H(x(T),u(T),p,adj(T),T) = */
	/* l(x,u,p,t)  */
	if (grampc->opt->IntegralCost == INT_ON) {
		lfct(&l, t[grampc->opt->Nhor - 1], x_, u_, p_, grampc->param->xdes, grampc->param->udes, grampc->userparam);
		if (grampc->opt->ScaleProblem == INT_ON) {
			scale_cost(&l, grampc->opt->JScale, 1);
		}
	}

	/* + adj' * f(x,u,p,t)   */
	ffct(s1, t[grampc->opt->Nhor - 1], x_, u_, p_, grampc->userparam);
	MatMult(&f, adj_, s1, 1, grampc->param->Nx, 1);

	/* + c(x,u,p,t)' * (mult + c(x,u,p,t)' * diag(pen) / 2); */
	if (grampc->opt->EqualityConstraints == INT_ON) {
		for (i = 0; i < grampc->param->Ng; i++) {
			k = i;
			c = c + cfct[k] * (mult[k] + cfct[k] * pen[k] / 2);
		}
	}
	if (grampc->opt->InequalityConstraints == INT_ON) {
		for (i = 0; i < grampc->param->Nh; i++) {
			k = i + grampc->param->Ng;
			c = c + cfct[k] * (mult[k] + cfct[k] * pen[k] / 2);
		}
	}

	/* dV(x,p,T)/dT */
	if (grampc->opt->TerminalCost == INT_ON) {
		dVdT(&dVdt, t[grampc->opt->Nhor - 1], x_, p_, grampc->param->xdes, grampc->userparam);
		if (grampc->opt->ScaleProblem == INT_ON) {
			scale_cost(&dVdt, grampc->opt->JScale, 1);
		}
	}

	/* + dc(x,p,T)/dT' * mult(T) */
	grampc->rws->gradT = l + f + c + dVdt + grampc->rws->dcdt;

	/* Scaling */
	if (grampc->opt->ScaleProblem == INT_ON) {
		grampc->rws->gradT = grampc->rws->gradT * grampc->opt->TScale;
	}
}

void timeproj(typeRNum *T, const typeGRAMPC *grampc)
{
	typeRNum Tmin_ = grampc->param->Tmin;
	typeRNum Tmax_ = grampc->param->Tmax;

	if (grampc->opt->ScaleProblem == INT_ON) {
		scale_time(&Tmin_, Tmin_, grampc);
		scale_time(&Tmax_, Tmax_, grampc);
	}

	/* lower bound */
	if (T[0] < Tmin_) {
		T[0] = Tmin_;
	}
	/* upper bound */
	else if (T[0] > Tmax_) {
		T[0] = Tmax_;
	}
}

void linesearch_adaptive(typeRNum *alpha, ctypeInt igrad, const typeGRAMPC *grampc)
{
	typeInt ils, j;
	ctypeInt NLS_adapt = 2 * (NALS + 1);
	typeRNum Jls[2] = { 0 , 0 };

	typeRNum *t = grampc->rws->t;
	typeRNum *tls = grampc->rws->tls;
	typeRNum *u = grampc->rws->u;
	typeRNum *uls = grampc->rws->uls;
	ctypeRNum *gradu = grampc->rws->gradu;
	typeRNum *p = grampc->rws->p;
	typeRNum *pls = grampc->rws->pls;
	ctypeRNum *gradp = grampc->rws->gradp;
	ctypeRNum *T = &grampc->rws->T;
	typeRNum Tls;
	ctypeRNum *gradT = &grampc->rws->gradT;
	typeRNum *lsAdapt = grampc->rws->lsAdapt + grampc->opt->MaxGradIter*NLS_adapt;

	/* Adapt interval */
	if (ABS(lsAdapt[NALS + 1] - lsAdapt[NALS + 3]) > grampc->opt->LineSearchAdaptAbsTol) {
		if (lsAdapt[NALS] >= lsAdapt[0] + (1 - grampc->opt->LineSearchIntervalTol)*(lsAdapt[NALS - 1] - lsAdapt[0])) {
			if (lsAdapt[NALS - 1] <= grampc->opt->LineSearchMax) {
				for (ils = 0; ils < NALS; ils++) {
					lsAdapt[ils] = lsAdapt[ils] * grampc->opt->LineSearchAdaptFactor;
				}
			}
			else {
				grampc->sol->status |= STATUS_LINESEARCH_MAX;
			}
		}
		else if (lsAdapt[NALS] <= lsAdapt[0] + grampc->opt->LineSearchIntervalTol*(lsAdapt[NALS - 1] - lsAdapt[0])) {
			if (lsAdapt[0] >= grampc->opt->LineSearchMin) {
				for (ils = 0; ils < NALS; ils++) {
					lsAdapt[ils] = lsAdapt[ils] / grampc->opt->LineSearchAdaptFactor;
				}
			}
			else {
				grampc->sol->status |= STATUS_LINESEARCH_MIN;
			}
		}
	}

	/* Evaluate cost for NLS step sizes */
	for (ils = 0; ils < NALS; ils++) {

		if (grampc->opt->OptimControl == INT_ON) {
			for (j = 0; j < grampc->opt->Nhor*grampc->param->Nu; j++) {
				uls[j] = u[j] - lsAdapt[ils] * gradu[j];
			}
			inputproj(uls, grampc);
		}
		else {
			uls = u;
		}

		if (grampc->opt->OptimParam == INT_ON) {
			for (j = 0; j < grampc->param->Np; j++) {
				pls[j] = p[j] - grampc->opt->OptimParamLineSearchFactor * lsAdapt[ils] * gradp[j];
			}
			paramproj(pls, grampc);
		}
		else {
			pls = p;
		}

		if (grampc->opt->OptimTime == INT_ON) {
			Tls = T[0] - grampc->opt->OptimTimeLineSearchFactor * lsAdapt[ils] * gradT[0];
			timeproj(&Tls, grampc);
			discretize_time(tls, Tls, grampc);
		}
		else {
			tls = t;
		}

		evaluate_sys(tls, uls, pls, grampc);
		evaluate_constraints(tls, uls, pls, 0, 0, grampc);
		evaluate_cost(Jls, tls, uls, pls, grampc);
		lsAdapt[ils + NALS + 1] = Jls[1];
	}

	/* Determine optimal step size by curve fitting */
	lsearch_fit(lsAdapt + NALS, lsAdapt + 2 * NALS + 1, lsAdapt, lsAdapt + NALS + 1);
	*alpha = lsAdapt[NALS];

	/* Save history */
	for (ils = 0; ils < 2 * (NALS + 1); ils++) {
		grampc->rws->lsAdapt[ils + igrad * NLS_adapt] = lsAdapt[ils];
	}

}

void linesearch_explicit(typeRNum *alpha, const typeGRAMPC *grampc)
{
	typeInt i, j;
	typeRNum *lsExplicit = grampc->rws->lsExplicit;
	typeRNum graduMax, graduAbs;
	typeRNum NomDenum[2];

	if (grampc->opt->OptimControl == INT_ON) {
		update_lsExplicit(lsExplicit, grampc->rws->u, grampc->rws->uprev, grampc->rws->gradu, grampc->rws->graduprev, grampc->opt->Nhor*grampc->param->Nu, grampc);
	}
	if (grampc->opt->OptimParam == INT_ON) {
		update_lsExplicit(NomDenum, grampc->rws->p, grampc->rws->pprev, grampc->rws->gradp, grampc->rws->gradpprev, grampc->param->Np, grampc);
		lsExplicit[0] += NomDenum[0] * grampc->opt->OptimParamLineSearchFactor;
		lsExplicit[1] += NomDenum[1] * grampc->opt->OptimParamLineSearchFactor * grampc->opt->OptimParamLineSearchFactor;
	}
	if (grampc->opt->OptimTime == INT_ON) {
		update_lsExplicit(NomDenum, &grampc->rws->T, &grampc->rws->Tprev, &grampc->rws->gradT, &grampc->rws->gradTprev, 1, grampc);
		lsExplicit[0] += NomDenum[0] * grampc->opt->OptimTimeLineSearchFactor;
		lsExplicit[1] += NomDenum[1] * grampc->opt->OptimTimeLineSearchFactor * grampc->opt->OptimTimeLineSearchFactor;
	}

	if (lsExplicit[0] > 0 && lsExplicit[1] > 0) {
		lsExplicit[2] = lsExplicit[0] / lsExplicit[1];
	}
	else {
		/* AutoFallback Method if the option is set OptimControl is on and limits for every input are set */
		if (grampc->opt->LineSearchExpAutoFallback == INT_ON && grampc->opt->OptimControl == INT_ON && lsExplicit[3] == 1) {
			/* Limit the stepsize to 10% of the maximum value */
			lsExplicit[2] = grampc->opt->LineSearchMax / 10;

			/* Set the stepsize, that the control update is not more than 1.0% of the control range */
			for (j = 0; j < grampc->param->Nu; j++) {
				graduMax = 0;
				for (i = 0; i < grampc->opt->Nhor; i++) {
					graduAbs = ABS(grampc->rws->gradu[j + i * grampc->param->Nu]);
					if (graduAbs > graduMax) {
						graduMax = graduAbs;
					}
				}
				if (grampc->opt->ScaleProblem == INT_ON) {
					graduMax = graduMax / grampc->opt->uScale[j];
				}
				lsExplicit[2] = MIN(lsExplicit[2], (grampc->param->umax[j] - grampc->param->umin[j]) / (graduMax * 100));
			}
		}
		else {
			lsExplicit[2] = grampc->opt->LineSearchInit;
			grampc->sol->status |= STATUS_LINESEARCH_INIT;
		}
	}

	if (lsExplicit[2] > grampc->opt->LineSearchMax) {
		lsExplicit[2] = grampc->opt->LineSearchMax;
		grampc->sol->status |= STATUS_LINESEARCH_MAX;
	}
	else if (lsExplicit[2] < grampc->opt->LineSearchMin) {
		lsExplicit[2] = grampc->opt->LineSearchMin;
		grampc->sol->status |= STATUS_LINESEARCH_MIN;
	}

	*alpha = lsExplicit[2];
}

void update_lsExplicit(typeRNum* NomDenum, ctypeRNum* a, ctypeRNum* aprev, ctypeRNum* dHda, ctypeRNum* dHdaprev, ctypeInt length, const typeGRAMPC *grampc) {
	typeRNum diffa, diffdHda;
	typeInt j;

	NomDenum[0] = 0;
	NomDenum[1] = 0;

	if (grampc->opt->LineSearchType == INT_EXPLS1) {
		/* Formula #1 */
		for (j = 0; j < length; j++) {
			diffa = a[j] - aprev[j];
			diffdHda = dHda[j] - dHdaprev[j];
			NomDenum[0] = NomDenum[0] + diffa * diffa;
			NomDenum[1] = NomDenum[1] + diffa * diffdHda;
		}
	}
	else {
		/* Formula #2 */
		for (j = 0; j < length; j++) {
			diffa = a[j] - aprev[j];
			diffdHda = dHda[j] - dHdaprev[j];
			NomDenum[0] = NomDenum[0] + diffa * diffdHda;
			NomDenum[1] = NomDenum[1] + diffdHda * diffdHda;
		}
	}
}

void evaluate_cost(typeRNum *s, ctypeRNum *t, ctypeRNum *u, ctypeRNum *p, const typeGRAMPC *grampc)
{
	typeInVfctPtr pIntCost;
	typeRNum Jint[2] = { 0, 0 };
	typeRNum Jterm[2] = { 0, 0 };
	ctypeRNum *p_ = p;

	/* integrator for integral cost function */
	if (grampc->opt->IntegratorCost == INT_TRAPZ) {
		pIntCost = &trapezodial;
	}
	else {
		pIntCost = &simpson;
	}

	/* Unscaling */
	if (grampc->opt->ScaleProblem == INT_ON) {
		ASSIGN_P(p_, p, grampc);
	}

	/* integrate cost */
	(*pIntCost)(Jint, t, grampc->rws->x, u, p_, grampc);

	/* terminal cost */
	WtermCost(Jterm, t[grampc->opt->Nhor - 1],
		grampc->rws->x + (grampc->opt->Nhor - 1) * grampc->param->Nx, p_,
		grampc->rws->mult + (grampc->opt->Nhor - 1) * grampc->param->Nc,
		grampc->rws->pen + (grampc->opt->Nhor - 1) * grampc->param->Nc,
		grampc->rws->cfct + (grampc->opt->Nhor - 1) * grampc->param->Nc, grampc);

	s[0] = Jint[0] + Jterm[0];
	s[1] = Jint[1] + Jterm[1];
}

void WintCost(typeRNum *s, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p_,
	ctypeRNum *mult, ctypeRNum *pen, ctypeRNum *cfct, const typeGRAMPC *grampc)
{
	typeInt i, k;
	ctypeRNum *x_ = x;
	ctypeRNum *u_ = u;

	/* Unscaling */
	if (grampc->opt->ScaleProblem == INT_ON) {
		ASSIGN_X(x_, x, grampc);
		ASSIGN_U(u_, u, grampc);
	}

	s[0] = 0;
	s[1] = 0;

	/* Integral Cost */
	if (grampc->opt->IntegralCost == INT_ON) {
		lfct(s, t, x_, u_, p_, grampc->param->xdes, grampc->param->udes, grampc->userparam);
		if (grampc->opt->ScaleProblem == INT_ON) {
			scale_cost(s, grampc->opt->JScale, 1);
		}
	}

	/* Equality Constraints */
	if (grampc->opt->EqualityConstraints == INT_ON) {
		for (i = 0; i < grampc->param->Ng; i++) {
			k = i;
            s[1] = s[1] + cfct[k] * (mult[k] + cfct[k] * pen[k] / 2);
		}
	}

	/* Inequality Constraints */
	if (grampc->opt->InequalityConstraints == INT_ON) {
		for (i = 0; i < grampc->param->Nh; i++) {
			k = grampc->param->Ng + i;
			s[1] = s[1] + cfct[k] * (mult[k] + cfct[k] * pen[k] / 2);
		}
	}
	s[1] = s[1] + s[0];
}

void WtermCost(typeRNum *s, ctypeRNum t, ctypeRNum *x, ctypeRNum *p_,
	ctypeRNum *mult, ctypeRNum *pen, ctypeRNum *cfct, const typeGRAMPC *grampc)
{
	typeInt i, k;
	ctypeRNum *x_ = x;

	/* Unscaling */
	if (grampc->opt->ScaleProblem == INT_ON) {
		ASSIGN_X(x_, x, grampc);
	}

	s[0] = 0;
	s[1] = 0;

	/* Terminal Cost */
	if (grampc->opt->TerminalCost == INT_ON) {
		Vfct(s, t, x_, p_, grampc->param->xdes, grampc->userparam);
		if (grampc->opt->ScaleProblem == INT_ON) {
			scale_cost(s, grampc->opt->JScale, 1);
		}
	}

	/* Terminal Equality Constraints */
	if (grampc->opt->TerminalEqualityConstraints == INT_ON) {
		for (i = 0; i < grampc->param->NgT; i++) {
			k = grampc->param->Ng + grampc->param->Nh + i;
            s[1] = s[1] + cfct[k] * (mult[k] + cfct[k] * pen[k] / 2);
		}
	}

	/* Terminal Inequality Constraints */
	if (grampc->opt->TerminalInequalityConstraints == INT_ON) {
		for (i = 0; i < grampc->param->NhT; i++) {
			k = grampc->param->Ng + grampc->param->Nh + grampc->param->NgT + i;
			s[1] = s[1] + cfct[k] * (mult[k] + cfct[k] * pen[k] / 2);
		}
	}
	s[1] = s[1] + s[0];
}

typeBoolean convergence_test_gradient(typeRNum rel_tol, const typeGRAMPC *grampc)
{
	typeRNum u_diff_norm = 0;
	typeRNum u_norm = 0;
	typeRNum u_rel_norm = 0;
	typeRNum p_diff_norm = 0;
	typeRNum p_norm = 0;
	typeRNum p_rel_norm = 0;
	typeRNum T_diff_norm = 0;
	typeRNum T_norm = 0;
	typeRNum T_rel_norm = 0;

	if (grampc->opt->OptimControl == INT_ON) {
		/* ||u - uprev|| */
		MatDiffNorm(&u_diff_norm, grampc->rws->u, grampc->rws->uprev, grampc->opt->Nhor, grampc->param->Nu);
		/* ||u|| */
		MatNorm(&u_norm, grampc->rws->u, grampc->opt->Nhor, grampc->param->Nu);
		u_rel_norm = u_norm > 0 ? u_diff_norm / u_norm : 0;
	}
	if (grampc->opt->OptimParam == INT_ON) {
		/* ||p - pprev|| */
		MatDiffNorm(&p_diff_norm, grampc->rws->p, grampc->rws->pprev, 1, grampc->param->Np);
		/* ||p|| */
		MatNorm(&p_norm, grampc->rws->p, 1, grampc->param->Np);
		p_rel_norm = p_norm > 0 ? p_diff_norm / p_norm : 0;
	}
	if (grampc->opt->OptimTime == INT_ON) {
		/* ||T - Tprev|| */
		T_diff_norm = (grampc->rws->T - grampc->rws->Tprev) * (grampc->rws->T - grampc->rws->Tprev);
		/* ||T|| */
		T_norm = grampc->rws->T * grampc->rws->T;
		T_rel_norm = T_norm > 0 ? SQRT(T_diff_norm / T_norm) : 0;
	}
	return MAX(MAX(u_rel_norm, p_rel_norm), T_rel_norm) < rel_tol;
}

typeBoolean convergence_test_constraints(ctypeRNum *abs_tol, const typeGRAMPC *grampc)
{
	typeInt i, j, k;
	ctypeRNum *cfct = grampc->rws->cfct;

	if (grampc->opt->EqualityConstraints == INT_ON) {
		for (j = 0; j < grampc->param->Ng; j++) {
			k = j;
			for (i = 1; i < grampc->opt->Nhor; i++) {
				if (ABS(cfct[i * grampc->param->Nc + k]) > abs_tol[k]) {
					return 0;
				}
			}
		}
	}
	if (grampc->opt->InequalityConstraints == INT_ON) {
		for (j = 0; j < grampc->param->Nh; j++) {
			k = grampc->param->Ng + j;
			for (i = 1; i < grampc->opt->Nhor; i++) {
				if (cfct[i * grampc->param->Nc + k] > abs_tol[k]) {
					return 0;
				}
			}
		}
	}

	i = grampc->opt->Nhor - 1;
	if (grampc->opt->TerminalEqualityConstraints == INT_ON) {
		for (j = 0; j < grampc->param->NgT; j++) {
			k = grampc->param->Ng + grampc->param->Nh + j;
			if (ABS(cfct[i * grampc->param->Nc + k]) > abs_tol[k]) {
				return 0;
			}
		}
	}
	if (grampc->opt->TerminalInequalityConstraints == INT_ON) {
		for (j = 0; j < grampc->param->NhT; j++) {
			k = grampc->param->Ng + grampc->param->Nh + grampc->param->NgT + j;
			if (cfct[i * grampc->param->Nc + k] > abs_tol[k]) {
				return 0;
			}
		}
	}
	return 1;
}

void shiftTrajectory(typeRNum *trajectory, ctypeInt Nhor, ctypeInt Nrows, ctypeInt Nshiftrows, ctypeRNum dt, ctypeRNum *t)
{
	/* function assumes uniform grid size over prediction horizon */
	typeInt i, j;

	/* compute how far the sampling points must be shifted */
	ctypeInt shift = (typeInt)(dt / (t[1] - t[0]));
	ctypeInt shiftvar = shift * Nrows;
	ctypeRNum interpfact = (dt / (t[1] - t[0])) - shift;

	if (shift >= Nhor) {
		grampc_error("Horizon too short for the current sampling time.");
	}

	/* Interpolation between the grid points */
	for (i = 0; i < Nhor - 1 - shift; i++) {
		for (j = 0; j < Nshiftrows; j++) {
			trajectory[j] = trajectory[shiftvar + j] + (trajectory[shiftvar + j + Nrows] - trajectory[shiftvar + j]) * interpfact;
		}
		trajectory += Nrows;
	}

	/* next element: extrapolation */
 /*
	if (shift == Nhor - 1) {
		for (j = 0; j < Nrows; j++) {
			trajectory[j] = trajectory[shiftvar + j] + (trajectory[shiftvar + j] - trajectory[shiftvar + j - Nrows]) * interpfact;
		}
	}
	else {
		for (j = 0; j < Nrows; j++) {
			trajectory[j] = trajectory[shiftvar + j] + (trajectory[shiftvar + j] - trajectory[j - Nrows]) / (1 / interpfact - 1);
		}
	}
	trajectory += Nrows;
	i++;
	*/

	/* if there are elements left hold last value */
	for (; i < Nhor; i++) {
		for (j = 0; j < Nshiftrows; j++) {
			trajectory[j] = trajectory[j - Nrows];
		}
		trajectory += Nrows;
	}
}

void shortenTrajectory(typeRNum *trajectory, ctypeInt Nhor, ctypeInt Nrows, ctypeInt Nshortenrows, ctypeRNum dt, ctypeRNum *t)
{
	/* function assumes uniform grid size over prediction horizon */
	typeInt i, j, ind;
	typeRNum takt;
	typeRNum dtgn = (t[Nhor - 1] - dt) / (Nhor - 1); /* compute new grid step size */

	/* Interpolation between the grid points from the beginning to the end */
	for (i = 0; i < Nhor - 1; i++) {
		takt = dt + dtgn * i;
		ind = (typeInt)(takt / (t[1] - t[0]));

		for (j = 0; j < Nshortenrows; j++) {
			trajectory[j + i * Nrows] = trajectory[ind*Nrows + j] + (trajectory[(ind + 1)* Nrows + j] - trajectory[ind*Nrows + j]) *(takt - t[ind]) / (t[1] - t[0]);
		}
	}
}
