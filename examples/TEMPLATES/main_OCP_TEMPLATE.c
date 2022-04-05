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


#include "grampc.h"
#include <time.h>

#define NX	
#define NU	
#define NC  
#define NP  

#define PRINTRES

#ifdef PRINTRES
void openFile(FILE **file, const char *name) {
	*file = fopen(name, "w");
	if (*file == NULL) {
		myPrint("Open file %s", name);
		printError(" failed!");
	}
}

void printNumVector2File(FILE *file, ctypeRNum *const val, ctypeInt size) {
	typeInt i;
	for (i = 0; i < size - 1; i++) {
		fprintf(file, "%.5f ,", val[i]);
	}
	fprintf(file, "%.5f;\n", val[size - 1]); /* new line */
}
void printIntVector2File(FILE *file, ctypeInt *const val, ctypeInt size) {
	typeInt i;
	for (i = 0; i < size - 1; i++) {
		fprintf(file, "%d ,", val[i]);
	}
	fprintf(file, "%d;\n", val[size - 1]); /* new line */
}
#endif

int main()
{
	typeGRAMPC *grampc;
	typeInt iOCP, i;
	typeInt MaxMultSimIter = 400;
	typeBoolean converged_grad = 0;
	typeBoolean converged_const = 0;

#ifdef PRINTRES
	FILE *file_x, *file_u, *file_cfct, *file_pen, *file_lag, *file_p, *file_T, *file_J, *file_Ncfct, *file_Npen,
		*file_iter, *file_status, *file_t, *file_veciter;
#endif

	clock_t tic, toc;
	typeRNum *CPUtimeVec;
	typeRNum CPUtime = 0;

	/********* Parameter definition *********/
	/* Initial values and setpoints of the states, inputs, parameters, penalties and Lagrangian mmultipliers, setpoints for the states and inputs */
	ctypeRNum x0[NX] = {};
	ctypeRNum xdes[NX] = {};

	/* Initial values, setpoints and limits of the inputs */
	ctypeRNum u0[NU] = {};
	ctypeRNum udes[NU] = {};
	ctypeRNum umax[NU] = {};
	ctypeRNum umin[NU] = {};

	/* Initial values and limits of the parameters */
	ctypeRNum p0[NP] = {};
	ctypeRNum pmax[NP] = {};
	ctypeRNum pmin[NP] = {};

	/* Time variables */
	ctypeRNum Thor = -10;  /* Prediction horizon */
	ctypeRNum Tmax = (typeRNum)1e8;
	ctypeRNum Tmin = (typeRNum)1e-8;

	ctypeRNum dt = (typeRNum)-1; /* Sampling time */
	typeRNum t = (typeRNum)0.0;  /* Initial time */

	/********* Option definition *********/
	/* Basic algorithmic options */
	ctypeInt Nhor = (typeInt)30;        /* Number of steps for the System Integration */
	ctypeInt MaxGradIter = (typeInt)2;  /* Maximum number of gradient Iterations */
	ctypeInt MaxMultIter = (typeInt)1;  /* Maximum number of augmented Lagrangain Interations */
	const char* ShiftControl = "on";

	/* Cost integration */
	const char* IntegralCost = "on";
	const char* TerminalCost = "on";
	const char* IntegratorCost = "trapezodial";

	/* System integration */
	const char* Integrator = "heun";
	ctypeRNum IntegratorRelTol = (typeRNum)1e-6;
	ctypeRNum IntegratorAbsTol = (typeRNum)1e-8;
	ctypeRNum IntegratorMinStepSize = EPS;
	ctypeInt  IntegratorMaxSteps = (typeInt)1e8;
	ctypeInt FlagsRodas[8] = { 0, 0, 0, NX, NX, 0, NX, NX };

	/* Line search */
	const char* LineSearchType = "explicit2";
	const char* LineSearchExpAutoFallback = "on";
	ctypeRNum LineSearchMax = (typeRNum)0.75;
	ctypeRNum LineSearchMin = (typeRNum)1e-10;
	ctypeRNum LineSearchInit = (typeRNum)1e-4;
	ctypeRNum LineSearchIntervalFactor = (typeRNum)0.85;
	ctypeRNum LineSearchAdaptFactor = (typeRNum)3.0 / 2.0;
	ctypeRNum LineSearchIntervalTol = (typeRNum)1e-1;

	/* Input and or parameter optimization  */
	const char* OptimControl = "on";
	const char* OptimParam = "off";
	ctypeRNum OptimParamLineSearchFactor = (typeRNum)1.0;
	const char* OptimTime = "on";
	ctypeRNum OptimTimeLineSearchFactor = (typeRNum)1.0;

	/* Scaling Values for the states, inputs and parameters */
	const char* ScaleProblem = "off";
	ctypeRNum xScale[NX] = { 1 };
	ctypeRNum xOffset[NX] = { 0 };
	ctypeRNum uScale[NU] = { 1 };
	ctypeRNum uOffset[NU] = { 0 };
	ctypeRNum pScale[NP] = { 1 };
	ctypeRNum pOffset[NP] = { 0 };
	ctypeRNum TScale = (typeRNum)1.0;
	ctypeRNum TOffset = (typeRNum)0.0;
	ctypeRNum JScale = (typeRNum)1.0;
	ctypeRNum cScale[NC] = { 1 };

	/* Tye of constraints' consideration */
	const char* EqualityConstraints = "on";
	const char* InequalityConstraints = "on";
	const char* TerminalEqualityConstraints = "on";
	const char* TerminalInequalityConstraints = "on";
	const char* ConstraintsHandling = "auglag";
	ctypeRNum ConstraintsAbsTol[NC] = { 1e-4 };

	/* Multipliers & penalties */
	ctypeRNum MultiplierMax = (typeRNum) 1e6;
	ctypeRNum MultiplierDampingFactor = 0;
	ctypeRNum PenaltyMax = (typeRNum)1e6;
	ctypeRNum PenaltyMin = (typeRNum)1e0;
	ctypeRNum PenaltyIncreaseFactor = (typeRNum)1.05;
	ctypeRNum PenaltyDecreaseFactor = (typeRNum)0.95;
	ctypeRNum PenaltyIncreaseThreshold = (typeRNum)1.0;
	ctypeRNum AugLagUpdateGradientRelTol = (typeRNum)1e-2;

	/* Convergences tests */
	const char* ConvergenceCheck = "on";
	ctypeRNum ConvergenceGradientRelTol = (typeRNum)1e-6;

	/********* userparam *********/
	typeRNum pCost[4] = { 0, 0, 0, 0 };
	typeUSERPARAM *userparam = pCost;

	/********* grampc init *********/
	grampc_init(&grampc, userparam);


	/********* set parameters and option *********/
	grampc_setparam_real_vector(grampc, "x0", x0);
	grampc_setparam_real_vector(grampc, "xdes", xdes);

	grampc_setparam_real_vector(grampc, "u0", u0);
	grampc_setparam_real_vector(grampc, "udes", udes);
	grampc_setparam_real_vector(grampc, "umax", umax);
	grampc_setparam_real_vector(grampc, "umin", umin);

	grampc_setparam_real_vector(grampc, "p0", p0);
	grampc_setparam_real_vector(grampc, "pmax", pmax);
	grampc_setparam_real_vector(grampc, "pmin", pmin);

	grampc_setparam_real(grampc, "Thor", Thor);
	grampc_setparam_real(grampc, "Tmax", Tmax);
	grampc_setparam_real(grampc, "Tmin", Tmin);

	grampc_setparam_real(grampc, "dt", dt);
	grampc_setparam_real(grampc, "t0", t);

	/********* Option definition *********/
	grampc_setopt_int(grampc, "Nhor", Nhor);
	grampc_setopt_int(grampc, "MaxGradIter", MaxGradIter);
	grampc_setopt_int(grampc, "MaxMultIter", MaxMultIter);
	grampc_setopt_string(grampc, "ShiftControl", ShiftControl);

	grampc_setopt_string(grampc, "IntegralCost", IntegralCost);
	grampc_setopt_string(grampc, "TerminalCost", TerminalCost);
	grampc_setopt_string(grampc, "IntegratorCost", IntegratorCost);

	grampc_setopt_string(grampc, "Integrator", Integrator);
	grampc_setopt_real(grampc, "IntegratorRelTol", IntegratorRelTol);
	grampc_setopt_real(grampc, "IntegratorAbsTol", IntegratorAbsTol);
	grampc_setopt_real(grampc, "IntegratorMinStepSize", IntegratorMinStepSize);
	grampc_setopt_int(grampc, "IntegratorMaxSteps", IntegratorMaxSteps);
	grampc_setopt_int_vector(grampc, "FlagsRodas", FlagsRodas);

	grampc_setopt_string(grampc, "LineSearchType", LineSearchType);
	grampc_setopt_string(grampc, "LineSearchExpAutoFallback", LineSearchExpAutoFallback);
	grampc_setopt_real(grampc, "LineSearchMax", LineSearchMax);
	grampc_setopt_real(grampc, "LineSearchMin", LineSearchMin);
	grampc_setopt_real(grampc, "LineSearchInit", LineSearchInit);
	grampc_setopt_real(grampc, "LineSearchIntervalFactor", LineSearchIntervalFactor);
	grampc_setopt_real(grampc, "LineSearchAdaptFactor", LineSearchAdaptFactor);
	grampc_setopt_real(grampc, "LineSearchIntervalTol", LineSearchIntervalTol);

	grampc_setopt_string(grampc, "OptimControl", OptimControl);
	grampc_setopt_string(grampc, "OptimParam", OptimParam);
	grampc_setopt_real(grampc, "OptimParamLineSearchFactor", OptimParamLineSearchFactor);
	grampc_setopt_string(grampc, "OptimTime", OptimTime);
	grampc_setopt_real(grampc, "OptimTimeLineSearchFactor", OptimTimeLineSearchFactor);

	grampc_setopt_string(grampc, "ScaleProblem", ScaleProblem);
	grampc_setopt_real_vector(grampc, "xScale", xScale);
	grampc_setopt_real_vector(grampc, "xOffset", xOffset);
	grampc_setopt_real_vector(grampc, "uScale", uScale);
	grampc_setopt_real_vector(grampc, "uOffset", uOffset);
	grampc_setopt_real_vector(grampc, "pScale", pScale);
	grampc_setopt_real_vector(grampc, "pOffset", pOffset);
	grampc_setopt_real(grampc, "TScale", TScale);
	grampc_setopt_real(grampc, "TOffset", TOffset);
	grampc_setopt_real(grampc, "JScale", JScale);
	grampc_setopt_real_vector(grampc, "cScale", cScale);

	grampc_setopt_string(grampc, "EqualityConstraints", EqualityConstraints);
	grampc_setopt_string(grampc, "InequalityConstraints", InequalityConstraints);
	grampc_setopt_string(grampc, "TerminalEqualityConstraints", TerminalEqualityConstraints);
	grampc_setopt_string(grampc, "TerminalInequalityConstraints", TerminalInequalityConstraints);
	grampc_setopt_string(grampc, "ConstraintsHandling", ConstraintsHandling);
	grampc_setopt_real_vector(grampc, "ConstraintsAbsTol", ConstraintsAbsTol);

	grampc_setopt_real(grampc, "MultiplierMax", MultiplierMax);
	grampc_setopt_real(grampc, "MultiplierDampingFactor", MultiplierDampingFactor);
	grampc_setopt_real(grampc, "PenaltyMax", PenaltyMax);
	grampc_setopt_real(grampc, "PenaltyMin", PenaltyMin);
	grampc_setopt_real(grampc, "PenaltyIncreaseFactor", PenaltyIncreaseFactor);
	grampc_setopt_real(grampc, "PenaltyDecreaseFactor", PenaltyDecreaseFactor);
	grampc_setopt_real(grampc, "PenaltyIncreaseThreshold", PenaltyIncreaseThreshold);
	grampc_setopt_real(grampc, "AugLagUpdateGradientRelTol", AugLagUpdateGradientRelTol);

	grampc_setopt_string(grampc, "ConvergenceCheck", ConvergenceCheck);
	grampc_setopt_real(grampc, "ConvergenceGradientRelTol", ConvergenceGradientRelTol);

#ifdef PRINTRES
	openFile(&file_x, "resOCP/xvec.txt");
	openFile(&file_u, "resOCP/uvec.txt");
	openFile(&file_cfct, "resOCP/cfctvec.txt");
	openFile(&file_lag, "resOCP/lagvec.txt");
	openFile(&file_pen, "resOCP/penvec.txt");
	openFile(&file_p, "resOCP/piter.txt");
	openFile(&file_T, "resOCP/Thoriter.txt");
	openFile(&file_J, "resOCP/Jiter.txt");
	openFile(&file_Ncfct, "resOCP/Ncfctiter.txt");
	openFile(&file_Npen, "resOCP/Npeniter.txt");
	openFile(&file_iter, "resOCP/itervec.txt");
	openFile(&file_status, "resOCP/status.txt");
	openFile(&file_t, "resOCP/tvec.txt");
	openFile(&file_veciter, "resOCP/multitervec.txt");
#endif

	CPUtimeVec = (typeRNum*)calloc(MaxMultSimIter + 1, sizeof(*CPUtimeVec));
	if (CPUtimeVec == NULL) {
		printError("Allocating memory for computationtime measurement failed");
	}

	grampc_printopt(grampc);
	grampc_printparam(grampc);


	printf("OCP running ...\n");
	for (iOCP = 0; iOCP <= MaxMultSimIter; iOCP++) {
		/* run grampc */
		tic = clock();
		grampc_run(grampc);
		toc = clock();
		CPUtimeVec[iOCP] = (typeRNum)((toc - tic) * 1000 / CLOCKS_PER_SEC);

		/* run convergence test */
		converged_grad = convergence_test_gradient(grampc->opt->ConvergenceGradientRelTol, grampc);
		if (converged_grad) {
			converged_const = convergence_test_constraints(grampc->opt->ConstraintsAbsTol, grampc);
		}

		/* check solver status */
		if (grampc->sol->status > 0) {
			if (grampc_printstatus(grampc->sol->status, STATUS_LEVEL_ERROR)) {
				myPrint("at iteration %i:\n -----\n", iOCP);
			}
		}


		/* print information of the multiplier iterations */
#ifdef PRINTRES
		/*printNumVector2File(file_p, grampc->sol->pnext, NP);*/
		printNumVector2File(file_T, &grampc->sol->Tnext, 1);
		printNumVector2File(file_J, grampc->sol->J, 2);
		printNumVector2File(file_Ncfct, &grampc->sol->cfct, 1);
		printNumVector2File(file_Npen, &grampc->sol->pen, 1);
		printIntVector2File(file_iter, grampc->sol->iter, grampc->opt->MaxMultIter);
		printIntVector2File(file_status, &grampc->sol->status, 1);
		printIntVector2File(file_veciter, &iOCP, 1);
#endif

		if (converged_const) {
			myPrint("Converged after %i multiplier iterations\n", iOCP);
			break;
		}
	}

	/* print result trajectories */
#ifdef PRINTRES
	for (i = 0; i < Nhor; i++) {
		printNumVector2File(file_x, grampc->rws->x + i * NX, NX);
		printNumVector2File(file_u, grampc->rws->u + i * NU, NU);
		printNumVector2File(file_cfct, grampc->rws->cfct + i * NC, NC);
		printNumVector2File(file_lag, grampc->rws->mult + i * NC, NC);
		printNumVector2File(file_pen, grampc->rws->pen + i * NC, NC);
		printNumVector2File(file_t, grampc->rws->t + i, 1);
	}

	fclose(file_x);
	fclose(file_u);
	fclose(file_p);
	fclose(file_T);
	fclose(file_J);
	fclose(file_cfct);
	fclose(file_pen);
	fclose(file_iter);
	fclose(file_status);
	fclose(file_t);
#endif

	for (i = 0; i <= MaxMultSimIter; i++) {
		CPUtime = CPUtime + CPUtimeVec[i];
	}

	grampc_free(&grampc);
	free(CPUtimeVec);

	printf("OCP finished. Average computation time: %.3f ms.\n", CPUtime);

	return 0;
}
