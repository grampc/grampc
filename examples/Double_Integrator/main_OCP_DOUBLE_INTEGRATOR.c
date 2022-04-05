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

#define NX	2
#define NU	1
#define NC  3
#define NP  0

#define PI 3.1415926535

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
		fprintf(file, "%.5e ,", val[i]);
	}
	fprintf(file, "%.5e;\n", val[size - 1]); /* new line */
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
	typeInt MaxMultSimIter = 100;
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
	/* Initial values and setpoints of the states, inputs, parameters, penalties and Lagrangian multipliers, setpoints for the states and inputs */
	ctypeRNum x0[NX] = { -1, -1 };
	ctypeRNum xdes[NX] = { 0, 0 };

	/* Initial values, setpoints and limits of the inputs */
	ctypeRNum u0[NU] = { 0 };
	ctypeRNum udes[NU] = { 0 };
	ctypeRNum umax[NU] = { 1.0 };
	ctypeRNum umin[NU] = { -1.0 };

	/* Time variables */
	ctypeRNum Thor = 5.25;  /* Prediction horizon */
	ctypeRNum Tmax = (typeRNum)10;
	ctypeRNum Tmin = (typeRNum)1.0;

	ctypeRNum dt = (typeRNum)0.01;  /* Sampling time */
	typeRNum t = (typeRNum)0.0;     /* Initial time */

	/********* Option definition *********/
	/* Basic algorithmic options */
	ctypeInt MaxGradIter = 200;    /* Maximum number of gradient Iterations */
	ctypeInt MaxMultIter = 1;      /* Maximum number of augmented Lagrangian iterations */
	ctypeInt Nhor = 50;
	const char* ShiftControl = "off";

	/* System integration */
	const char* Integrator = "euler";

	/* Line search */
	ctypeRNum LineSearchMax = 1e2;

	/* Input and or parameter optimization  */
	const char* OptimTime = "off";

	/* Penalties */
	ctypeRNum PenaltyMin = 1e1;
	ctypeRNum PenaltyIncreaseFactor = 1.25;
	ctypeRNum PenaltyDecreaseFactor = 1.0;

	/* Constraints thresholds */
	ctypeRNum ConstraintsAbsTol[3] = { 1e-6, 1e-6, 1e-6 };

	/* Convergences tests */
	const char* ConvergenceCheck = "on";
	ctypeRNum ConvergenceGradientRelTol = (typeRNum)1e-9;

	/********* userparam *********/
	typeRNum pCost[2] = { 0.1, 1 };
	typeUSERPARAM *userparam = pCost;

	/********* grampc init *********/
	grampc_init(&grampc, userparam);


	/********* set parameters *********/
	grampc_setparam_real_vector(grampc, "x0", x0);
	grampc_setparam_real_vector(grampc, "xdes", xdes);

	grampc_setparam_real_vector(grampc, "u0", u0);
	grampc_setparam_real_vector(grampc, "udes", udes);
	grampc_setparam_real_vector(grampc, "umax", umax);
	grampc_setparam_real_vector(grampc, "umin", umin);

	grampc_setparam_real(grampc, "Thor", Thor);
	grampc_setparam_real(grampc, "Tmax", Tmax);
	grampc_setparam_real(grampc, "Tmin", Tmin);

	grampc_setparam_real(grampc, "dt", dt);
	grampc_setparam_real(grampc, "t0", t);

	/********* Option definition *********/
	grampc_setopt_int(grampc, "MaxGradIter", MaxGradIter);
	grampc_setopt_int(grampc, "MaxMultIter", MaxMultIter);
	grampc_setopt_int(grampc, "Nhor", Nhor);
	grampc_setopt_string(grampc, "ShiftControl", ShiftControl);


	grampc_setopt_string(grampc, "Integrator", Integrator);

	grampc_setopt_real(grampc, "LineSearchMax", LineSearchMax);

	grampc_setopt_string(grampc, "OptimTime", OptimTime);

	grampc_setopt_real(grampc, "PenaltyMin", PenaltyMin);
	grampc_setopt_real(grampc, "PenaltyIncreaseFactor", PenaltyIncreaseFactor);
	grampc_setopt_real(grampc, "PenaltyDecreaseFactor", PenaltyDecreaseFactor);

	grampc_setopt_real_vector(grampc, "ConstraintsAbsTol", ConstraintsAbsTol);

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
			myPrint("Converged after %i multiplier iterations\n", iOCP + 1);
			break;
		}
	}

	/* print result trajectories */
#ifdef PRINTRES
	for (i = 0; i < grampc->opt->Nhor; i++) {
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
