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


#include "grampc_setopt.h"

#define COMPARE_INCLUSIVE 1
#define COMPARE_EXCLUSIVE 0

void setNumOpt(typeRNum *cs, const typeChar *optName, ctypeRNum optValue, ctypeRNum lb, ctypeRNum ub, ctypeInt cmp) {
	if ((cmp == COMPARE_INCLUSIVE && lb <= optValue && optValue <= ub) ||
		(cmp == COMPARE_EXCLUSIVE && lb < optValue && optValue < ub)) {
		*cs = optValue;
	}
	else {
		grampc_error_addstring(INVALID_OPTION_VALUE, optName);
	}
}

void setIntOpt(typeInt *cs, const typeChar *optName, ctypeInt optValue) {
	if (optValue > 0) {
		*cs = optValue;
	}
	else {
		grampc_error_addstring(INVALID_OPTION_VALUE, optName);
	}
}

void setOnOffOpt(typeInt *cs, const typeChar *optName, const typeChar *optValue) {
	if (!strcmp(optValue, "on")) {
		*cs = INT_ON;
	}
	else if (!strcmp(optValue, "off")) {
		*cs = INT_OFF;
	}
	else {
		grampc_error_addstring(INVALID_OPTION_VALUE, optName);
	}
}

void grampc_setopt_real(const typeGRAMPC *grampc, const typeChar *optName, ctypeRNum optValue)
{
	/* Integrator relative tolerance */
	if (!strcmp(optName, "IntegratorRelTol")) {
		setNumOpt(&grampc->opt->IntegratorRelTol, optName, optValue, 0.0, INF, COMPARE_EXCLUSIVE);
	}
	/* Integrator absolute tolerance */
	else if (!strcmp(optName, "IntegratorAbsTol")) {
		setNumOpt(&grampc->opt->IntegratorAbsTol, optName, optValue, 0.0, INF, COMPARE_EXCLUSIVE);
	}
	/* minimal step size of integrator ruku45 */
	else if (!strcmp(optName, "IntegratorMinStepSize")) {
		setNumOpt(&grampc->opt->IntegratorMinStepSize, optName, optValue, 0.0, INF, COMPARE_EXCLUSIVE);
	}
	/* Max value linesearch */
	else if (!strcmp(optName, "LineSearchMax")) {
		setNumOpt(&grampc->opt->LineSearchMax, optName, optValue, 0.0, INF, COMPARE_EXCLUSIVE);
	}
	/* Min value linesearch */
	else if (!strcmp(optName, "LineSearchMin")) {
		setNumOpt(&grampc->opt->LineSearchMin, optName, optValue, 0.0, INF, COMPARE_EXCLUSIVE);
	}
	/* Init value linesearch */
	else if (!strcmp(optName, "LineSearchInit")) {
		setNumOpt(&grampc->opt->LineSearchInit, optName, optValue, 0.0, INF, COMPARE_EXCLUSIVE);
		init_rws_linesearch(grampc);
	}
	/* LineSearchIntervalFactor */
	else if (!strcmp(optName, "LineSearchIntervalFactor")) {
		setNumOpt(&grampc->opt->LineSearchIntervalFactor, optName, optValue, 0.0, 1.0, COMPARE_EXCLUSIVE);
		init_rws_linesearch(grampc);
	}
	/* LineSearchAdaptFactor */
	else if (!strcmp(optName, "LineSearchAdaptFactor")) {
		setNumOpt(&grampc->opt->LineSearchAdaptFactor, optName, optValue, 1.0, INF, COMPARE_EXCLUSIVE);
	}
	/* LineSearchIntervalTol */
	else if (!strcmp(optName, "LineSearchIntervalTol")) {
		setNumOpt(&grampc->opt->LineSearchIntervalTol, optName, optValue, 0.0, 1.0, COMPARE_EXCLUSIVE);
	}
	/* OptimParamLineSearchFactor */
	else if (!strcmp(optName, "OptimParamLineSearchFactor")) {
		setNumOpt(&grampc->opt->OptimParamLineSearchFactor, optName, optValue, 0.0, INF, COMPARE_EXCLUSIVE);
	}
	/* OptimTimeLineSearchFactor */
	else if (!strcmp(optName, "OptimTimeLineSearchFactor")) {
		setNumOpt(&grampc->opt->OptimTimeLineSearchFactor, optName, optValue, 0.0, INF, COMPARE_EXCLUSIVE);
	}
	/* TScale */
	else if (!strcmp(optName, "TScale")) {
		setNumOpt(&grampc->opt->TScale, optName, optValue, 0.0, INF, COMPARE_EXCLUSIVE);
		init_rws_time(grampc);
	}
	/* TOffset */
	else if (!strcmp(optName, "TOffset")) {
		setNumOpt(&grampc->opt->TOffset, optName, optValue, -INF, INF, COMPARE_EXCLUSIVE);
		init_rws_time(grampc);
	}
	/* JScale */
	else if (!strcmp(optName, "JScale")) {
		setNumOpt(&grampc->opt->JScale, optName, optValue, 0.0, INF, COMPARE_EXCLUSIVE);
	}
	/* MultiplierMax */
	else if (!strcmp(optName, "MultiplierMax")) {
		setNumOpt(&grampc->opt->MultiplierMax, optName, optValue, 0.0, INF, COMPARE_EXCLUSIVE);
	}
	/* MultiplierDampingFactor */
	else if (!strcmp(optName, "MultiplierDampingFactor")) {
		setNumOpt(&grampc->opt->MultiplierDampingFactor, optName, optValue, 0.0, 1.0, COMPARE_INCLUSIVE);
	}
	/* PenaltyMax */
	else if (!strcmp(optName, "PenaltyMax")) {
		setNumOpt(&grampc->opt->PenaltyMax, optName, optValue, 0.0, INF, COMPARE_EXCLUSIVE);
	}
	/* PenaltyMin */
	else if (!strcmp(optName, "PenaltyMin")) {
		setNumOpt(&grampc->opt->PenaltyMin, optName, optValue, 0.0, INF, COMPARE_EXCLUSIVE);
		init_rws_multipliers(grampc);
	}
	/* PenaltyIncreaseFactor */
	else if (!strcmp(optName, "PenaltyIncreaseFactor")) {
		setNumOpt(&grampc->opt->PenaltyIncreaseFactor, optName, optValue, 1.0, INF, COMPARE_INCLUSIVE);
	}
	/* PenaltyDecreaseFactor */
	else if (!strcmp(optName, "PenaltyDecreaseFactor")) {
		setNumOpt(&grampc->opt->PenaltyDecreaseFactor, optName, optValue, 0.0, 1.0, COMPARE_INCLUSIVE);
	}
	/* PenaltyIncreaseThreshold */
	else if (!strcmp(optName, "PenaltyIncreaseThreshold")) {
		setNumOpt(&grampc->opt->PenaltyIncreaseThreshold, optName, optValue, 0.0, INF, COMPARE_INCLUSIVE);
	}
	/* AugLagUpdateGradientRelTol */
	else if (!strcmp(optName, "AugLagUpdateGradientRelTol")) {
		setNumOpt(&grampc->opt->AugLagUpdateGradientRelTol, optName, optValue, 0.0, 1, COMPARE_INCLUSIVE);
	}
	/* ConvergenceGradientRelTol */
	else if (!strcmp(optName, "ConvergenceGradientRelTol")) {
		setNumOpt(&grampc->opt->ConvergenceGradientRelTol, optName, optValue, 0.0, 1, COMPARE_INCLUSIVE);
	}
	/* Undefined optName */
	else {
		grampc_error_addstring(INVALID_OPTION_NAME, optName);
	}
}

void grampc_setopt_int(const typeGRAMPC *grampc, const typeChar *optName, ctypeInt optValue)
{
	/* MaxGradIter */
	if (!strcmp(optName, "MaxGradIter")) {
		setIntOpt(&grampc->opt->MaxGradIter, optName, optValue);
		/* Reallocaton of lsAdapt */
		resizeNumMatrix(&grampc->rws->lsAdapt, 2 * (NALS + 1)*(1 + grampc->opt->MaxGradIter));
		init_rws_linesearch(grampc);
	}
	/* MaxMultIter */
	else if (!strcmp(optName, "MaxMultIter")) {
		setIntOpt(&grampc->opt->MaxMultIter, optName, optValue);
		/* Reallocation of iter */
		resizeIntMatrix(&grampc->sol->iter, grampc->opt->MaxMultIter);
	}
	else if (!strcmp(optName, "Nhor")) {
		if (optValue > 1) {
			grampc->opt->Nhor = optValue;
		}
		else {
			grampc_error_addstring(INVALID_OPTION_VALUE, optName);
		}

		/* Reallocation of all fields that depend on Nhor */
		resizeNumMatrix(&grampc->rws->t, grampc->opt->Nhor);
		resizeNumMatrix(&grampc->rws->tls, grampc->opt->Nhor);

		resizeNumMatrix(&grampc->rws->x, grampc->opt->Nhor * grampc->param->Nx);
		resizeNumMatrix(&grampc->rws->adj, grampc->opt->Nhor * grampc->param->Nx);
		resizeNumMatrix(&grampc->rws->dcdx, (grampc->opt->Nhor + 1) * grampc->param->Nx);

		resizeNumMatrix(&grampc->rws->u, grampc->opt->Nhor * grampc->param->Nu);
		resizeNumMatrix(&grampc->rws->uls, grampc->opt->Nhor * grampc->param->Nu);
		resizeNumMatrix(&grampc->rws->uprev, grampc->opt->Nhor * grampc->param->Nu);
		resizeNumMatrix(&grampc->rws->gradu, grampc->opt->Nhor * grampc->param->Nu);
		resizeNumMatrix(&grampc->rws->graduprev, grampc->opt->Nhor * grampc->param->Nu);
		resizeNumMatrix(&grampc->rws->dcdu, grampc->opt->Nhor * grampc->param->Nu);

		resizeNumMatrix(&grampc->rws->dcdp, (grampc->opt->Nhor + 1) * grampc->param->Np);

		resizeNumMatrix(&grampc->rws->mult, grampc->opt->Nhor * grampc->param->Nc);
		resizeNumMatrix(&grampc->rws->pen, grampc->opt->Nhor * grampc->param->Nc);
		resizeNumMatrix(&grampc->rws->cfct, grampc->opt->Nhor * grampc->param->Nc);
		resizeNumMatrix(&grampc->rws->cfctprev, grampc->opt->Nhor * grampc->param->Nc);

		/* Initialization over Nhor */
		init_rws_time(grampc);
		init_rws_controls(grampc);
		init_rws_multipliers(grampc);
		resize_rwsRodas(grampc);
	}
	else if (!strcmp(optName, "IntegratorMaxSteps")) {
		setIntOpt(&grampc->opt->IntegratorMaxSteps, optName, optValue);
	}
	/* Undefined optName */
	else {
		grampc_error_addstring(INVALID_OPTION_NAME, optName);
	}
}

void grampc_setopt_string(const typeGRAMPC *grampc, const typeChar *optName, const typeChar *optValue)
{
	/* ControlShift */
	if (!strcmp(optName, "ShiftControl")) {
		setOnOffOpt(&grampc->opt->ShiftControl, optName, optValue);
	}
	/* ProblemScale */
	else if (!strcmp(optName, "ScaleProblem")) {
		setOnOffOpt(&grampc->opt->ScaleProblem, optName, optValue);
		init_rws_controls(grampc);
		init_rws_parameters(grampc);
		init_rws_time(grampc);
		init_rws_constraints(grampc);
	}
	/* TimeDiscretization */
	else if (!strcmp(optName, "TimeDiscretization")) {
		if (!strcmp(optValue, "uniform")) {
			grampc->opt->TimeDiscretization = INT_UNIFORM;
		}
		else if (!strcmp(optValue, "nonuniform")) {
			grampc->opt->TimeDiscretization = INT_NONUNIFORM;
		}
		else {
			grampc_error_addstring(INVALID_OPTION_VALUE, optName);
		}
		init_rws_time(grampc);
	}
	/* IntegratorCost */
	else if (!strcmp(optName, "IntegratorCost")) {
		if (!strcmp(optValue, "trapezodial")) {
			grampc->opt->IntegratorCost = INT_TRAPZ;
		}
		else if (!strcmp(optValue, "simpson")) {
			grampc->opt->IntegratorCost = INT_SIMPSON;
		}
		else {
			grampc_error_addstring(INVALID_OPTION_VALUE, optName);
		}
		resize_rwsGeneral(grampc);
	}
	/* Integrator type */
	else if (!strcmp(optName, "Integrator")) {
		if (!strcmp(optValue, "euler")) {
			grampc->opt->Integrator = INT_EULER;
		}
		else if (!strcmp(optValue, "modeuler")) {
			grampc->opt->Integrator = INT_MODEULER;
		}
		else if (!strcmp(optValue, "heun")) {
			grampc->opt->Integrator = INT_HEUN;
		}
		else if (!strcmp(optValue, "rodas")) {
			grampc->opt->Integrator = INT_RODAS;
		}
		else if (!strcmp(optValue, "ruku45")) {
			grampc->opt->Integrator = INT_RUKU45;
		}
		else {
			grampc_error_addstring(INVALID_OPTION_VALUE, optName);
		}
		resize_rwsGeneral(grampc);
		resize_rwsRodas(grampc);
	}
	/* LineSearchType */
	else if (!strcmp(optName, "LineSearchType")) {
		if (!strcmp(optValue, "adaptive")) {
			grampc->opt->LineSearchType = INT_ADAPTIVELS;
		}
		else if (!strcmp(optValue, "explicit1")) {
			grampc->opt->LineSearchType = INT_EXPLS1;
		}
		else if (!strcmp(optValue, "explicit2")) {
			grampc->opt->LineSearchType = INT_EXPLS2;
		}
		else {
			grampc_error_addstring(INVALID_OPTION_VALUE, optName);
		}
		init_rws_linesearch(grampc);
	}
	/* LineSearchExpAutoFallback */
	else if (!strcmp(optName, "LineSearchExpAutoFallback")) {
		setOnOffOpt(&grampc->opt->LineSearchExpAutoFallback, optName, optValue);
	}
	/* activate/deactivate control Optimization */
	else if (!strcmp(optName, "OptimControl")) {
		setOnOffOpt(&grampc->opt->OptimControl, optName, optValue);
	}
	/* activate/deactivate parameter Optimization */
	else if (!strcmp(optName, "OptimParam")) {
		setOnOffOpt(&grampc->opt->OptimParam, optName, optValue);
	}
	/* activate/deactivate time Optimization */
	else if (!strcmp(optName, "OptimTime")) {
		setOnOffOpt(&grampc->opt->OptimTime, optName, optValue);
	}
	/* activate/deactivate IntegralCost */
	else if (!strcmp(optName, "IntegralCost")) {
		setOnOffOpt(&grampc->opt->IntegralCost, optName, optValue);
	}
	/* activate/deactivate TerminalCost */
	else if (!strcmp(optName, "TerminalCost")) {
		setOnOffOpt(&grampc->opt->TerminalCost, optName, optValue);
	}
	/* EqualityConstraints */
	else if (!strcmp(optName, "EqualityConstraints")) {
		setOnOffOpt(&grampc->opt->EqualityConstraints, optName, optValue);
	}
	/* InequalityConstraints */
	else if (!strcmp(optName, "InequalityConstraints")) {
		setOnOffOpt(&grampc->opt->InequalityConstraints, optName, optValue);
	}
	/* TerminalEqualityConstraints */
	else if (!strcmp(optName, "TerminalEqualityConstraints")) {
		setOnOffOpt(&grampc->opt->TerminalEqualityConstraints, optName, optValue);
	}
	/* TerminalInequalityConstraints */
	else if (!strcmp(optName, "TerminalInequalityConstraints")) {
		setOnOffOpt(&grampc->opt->TerminalInequalityConstraints, optName, optValue);
	}
	/* choose between external penalty and augmented lagrangian */
	else if (!strcmp(optName, "ConstraintsHandling")) {
		if (!strcmp(optValue, "extpen")) {
			grampc->opt->ConstraintsHandling = INT_EXTPEN;
		}
		else if (!strcmp(optValue, "auglag")) {
			grampc->opt->ConstraintsHandling = INT_AUGLAG;
		}
		else {
			grampc_error_addstring(INVALID_OPTION_VALUE, optName);
		}
		init_rws_multipliers(grampc);
	}
	/* activate/deactivate convergence test */
	else if (!strcmp(optName, "ConvergenceCheck")) {
		setOnOffOpt(&grampc->opt->ConvergenceCheck, optName, optValue);
	}
	/* Undefined optName */
	else {
		grampc_error_addstring(INVALID_OPTION_NAME, optName);
	}
}

void grampc_setopt_real_vector(const typeGRAMPC *grampc, const typeChar *optName, ctypeRNum *optValue)
{
	/* xScale */
	if (!strcmp(optName, "xScale")) {
		MatCopy(grampc->opt->xScale, optValue, 1, grampc->param->Nx);
	}
	/* xOffset */
	else if (!strcmp(optName, "xOffset")) {
		MatCopy(grampc->opt->xOffset, optValue, 1, grampc->param->Nx);
	}
	/* uScale */
	else if (!strcmp(optName, "uScale")) {
		MatCopy(grampc->opt->uScale, optValue, 1, grampc->param->Nu);
		init_rws_controls(grampc);
	}
	/* uOffset */
	else if (!strcmp(optName, "uOffset")) {
		MatCopy(grampc->opt->uOffset, optValue, 1, grampc->param->Nu);
		init_rws_controls(grampc);
	}
	/* pScale */
	else if (!strcmp(optName, "pScale")) {
		MatCopy(grampc->opt->pScale, optValue, 1, grampc->param->Np);
		init_rws_parameters(grampc);
	}
	/* pOffset */
	else if (!strcmp(optName, "pOffset")) {
		MatCopy(grampc->opt->pOffset, optValue, 1, grampc->param->Np);
		init_rws_parameters(grampc);
	}
	/* cScale */
	else if (!strcmp(optName, "cScale")) {
		MatCopy(grampc->opt->cScale, optValue, 1, grampc->param->Nc);
		init_rws_constraints(grampc);
	}
	/* ConstraintsAbsTol */
	else if (!strcmp(optName, "ConstraintsAbsTol")) {
		MatCopy(grampc->opt->ConstraintsAbsTol, optValue, 1, grampc->param->Nc);
		init_rws_constraints(grampc);
	}
	/* Undefined optName */
	else {
		grampc_error_addstring(INVALID_OPTION_NAME, optName);
	}
}

void grampc_setopt_int_vector(const typeGRAMPC *grampc, const typeChar *optName, ctypeInt *optValue)
{
	/* FlagsRodas */
	if (!strcmp(optName, "FlagsRodas")) {
		memcpy(grampc->opt->FlagsRodas, optValue, 8 * sizeof(*grampc->opt->FlagsRodas));
		setLWorkRodas(grampc);
	}
	/* Undefined optName */
	else {
		grampc_error_addstring(INVALID_OPTION_NAME, optName);
	}
}

void grampc_printopt(const typeGRAMPC *grampc)
{
	myPrint("%s", "-- MPC OPTIONS --\n");
	myPrint("                         Nhor: %d\n", grampc->opt->Nhor);
	myPrint("                  MaxGradIter: %d\n", grampc->opt->MaxGradIter);
	myPrint("                  MaxMultIter: %d\n", grampc->opt->MaxMultIter);
	myPrint("                 ShiftControl: %s\n", grampc->opt->ShiftControl == INT_ON ? "on" : "off");

	myPrint("           TimeDiscretization: %s\n", grampc->opt->TimeDiscretization == INT_UNIFORM ? "uniform" : "nonuniform");

	myPrint("                 IntegralCost: %s\n", grampc->opt->IntegralCost == INT_ON ? "on" : "off");
	myPrint("                 TerminalCost: %s\n", grampc->opt->TerminalCost == INT_ON ? "on" : "off");
	myPrint("               IntegratorCost: %s\n", grampc->opt->IntegratorCost == INT_TRAPZ ? "trapezodial" : "simpson");

	myPrint("                   Integrator: %s\n", IntegratorInt2Str(grampc->opt->Integrator));
	myPrint("             IntegratorRelTol: %.2e\n", grampc->opt->IntegratorRelTol);
	myPrint("             IntegratorAbsTol: %.2e\n", grampc->opt->IntegratorAbsTol);
	myPrint("        IntegratorMinStepSize: %.2e\n", grampc->opt->IntegratorMinStepSize);
	myPrint("           IntegratorMaxSteps: %.2e\n", (typeRNum)grampc->opt->IntegratorMaxSteps);

	if (grampc->opt->Integrator == INT_RODAS) {
		myPrint("                         IFCN: %d\n", grampc->opt->FlagsRodas[0]);
		myPrint("                         IDFX: %d\n", grampc->opt->FlagsRodas[1]);
		myPrint("                         IJAC: %d\n", grampc->opt->FlagsRodas[2]);
		myPrint("                        MLJAC: %d\n", grampc->opt->FlagsRodas[4]);
		myPrint("                        MUJAC: %d\n", grampc->opt->FlagsRodas[5]);
		myPrint("                         IMAS: %d\n", grampc->opt->FlagsRodas[3]);
		myPrint("                        MLMAS: %d\n", grampc->opt->FlagsRodas[6]);
		myPrint("                        MUMAS: %d\n", grampc->opt->FlagsRodas[7]);
	}

	myPrint("               LineSearchType: %s\n", LineSearchTypeInt2Str(grampc->opt->LineSearchType));
	myPrint("    LineSearchExpAutoFallback: %s\n", grampc->opt->LineSearchExpAutoFallback == INT_ON ? "on" : "off");
	myPrint("                LineSearchMax: %.2e\n", grampc->opt->LineSearchMax);
	myPrint("                LineSearchMin: %.2e\n", grampc->opt->LineSearchMin);
	myPrint("               LineSearchInit: %.2e\n", grampc->opt->LineSearchInit);
	myPrint("     LineSearchIntervalFactor: %.3f\n", grampc->opt->LineSearchIntervalFactor);
	myPrint("        LineSearchAdaptFactor: %.3f\n", grampc->opt->LineSearchAdaptFactor);
	myPrint("        LineSearchIntervalTol: %.3f\n", grampc->opt->LineSearchIntervalTol);

	myPrint("                 OptimControl: %s\n", grampc->opt->OptimControl == INT_ON ? "on" : "off");
	myPrint("                   OptimParam: %s\n", grampc->opt->OptimParam == INT_ON ? "on" : "off");
	myPrint("   OptimParamLineSearchFactor: %.3f\n", grampc->opt->OptimParamLineSearchFactor);
	myPrint("                    OptimTime: %s\n", grampc->opt->OptimTime == INT_ON ? "on" : "off");
	myPrint("    OptimTimeLineSearchFactor: %.3f\n", grampc->opt->OptimTimeLineSearchFactor);

	myPrint("                 ScaleProblem: %s\n", grampc->opt->ScaleProblem == INT_ON ? "on" : "off");
	print_vector("                       xScale: ", grampc->opt->xScale, grampc->param->Nx);
	print_vector("                      xOffset: ", grampc->opt->xOffset, grampc->param->Nx);
	print_vector("                       uScale: ", grampc->opt->uScale, grampc->param->Nu);
	print_vector("                      uOffset: ", grampc->opt->uOffset, grampc->param->Nu);
	print_vector("                       pScale: ", grampc->opt->pScale, grampc->param->Np);
	print_vector("                      pOffset: ", grampc->opt->pOffset, grampc->param->Np);
	myPrint("                       TScale: %.3f\n", grampc->opt->TScale);
	myPrint("                      TOffset: %.3f\n", grampc->opt->TOffset);
	myPrint("                       JScale: %.3f\n", grampc->opt->JScale);
	print_vector("                       cScale: ", grampc->opt->cScale, grampc->param->Nc);

	myPrint("          EqualityConstraints: %s\n", grampc->opt->EqualityConstraints == INT_ON ? "on" : "off");
	myPrint("        InequalityConstraints: %s\n", grampc->opt->InequalityConstraints == INT_ON ? "on" : "off");
	myPrint("  TerminalEqualityConstraints: %s\n", grampc->opt->TerminalEqualityConstraints == INT_ON ? "on" : "off");
	myPrint("TerminalInequalityConstraints: %s\n", grampc->opt->TerminalInequalityConstraints == INT_ON ? "on" : "off");
	myPrint("          ConstraintsHandling: %s\n", grampc->opt->ConstraintsHandling == INT_EXTPEN ? "extpen" : "auglag");
	print_vector("             ConstraintAbsTol: ", grampc->opt->ConstraintsAbsTol, grampc->param->Nc);

	myPrint("                MultiplierMax: %.2e\n", grampc->opt->MultiplierMax);
	myPrint("      MultiplierDampingFactor: %.3f\n", grampc->opt->MultiplierDampingFactor);
	myPrint("                   PenaltyMax: %.2e\n", grampc->opt->PenaltyMax);
	myPrint("                   PenaltyMin: %.2e\n", grampc->opt->PenaltyMin);
	myPrint("        PenaltyIncreaseFactor: %.3f\n", grampc->opt->PenaltyIncreaseFactor);
	myPrint("        PenaltyDecreaseFactor: %.3f\n", grampc->opt->PenaltyDecreaseFactor);
	myPrint("     PenaltyIncreaseThreshold: %.3f\n", grampc->opt->PenaltyIncreaseThreshold);
	myPrint("   AugLagUpdateGradientRelTol: %.2e\n", grampc->opt->AugLagUpdateGradientRelTol);

	myPrint("             ConvergenceCheck: %s\n", grampc->opt->ConvergenceCheck == INT_ON ? "on" : "off");
	myPrint("    ConvergenceGradientRelTol: %.2e\n", grampc->opt->ConvergenceGradientRelTol);

}
