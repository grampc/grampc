/* This file is part of GRAMPC - (https://github.com/grampc/grampc)
 *
 * GRAMPC -- A software framework for embedded nonlinear model predictive
 * control using a gradient-based augmented Lagrangian approach
 *
 * Copyright 2014-2025 by Knut Graichen, Andreas Voelz, Thore Wietzke,
 * Tobias Englert (<v2.3), Felix Mesmer (<v2.3), Soenke Rhein (<v2.3),
 * Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 * All rights reserved.
 *
 * GRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
 *
 */


#include "grampc.h"
#include "timing.h"

#define N 4 /* number of chain elements */

#define NX	(3*(2*N-1))
#define NU	3
#define NC  0
#define NP  0

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
		fprintf(file, "%.5f, ", val[i]);
	}
	fprintf(file, "%.5f\n", val[size - 1]); /* new line */
}
void printIntVector2File(FILE *file, ctypeInt *const val, ctypeInt size) {
	typeInt i;
	for (i = 0; i < size - 1; i++) {
		fprintf(file, "%d, ", val[i]);
	}
	fprintf(file, "%d\n", val[size - 1]); /* new line */
}
#endif

int main()
{
	typeGRAMPC *grampc;
	typeInt iMPC, MaxSimIter, i;
	ctypeRNum Tsim = 50.0;

#ifdef PRINTRES
	FILE *file_x, *file_u, *file_p, *file_T, *file_J, *file_Ncfct, *file_Npen, *file_iter, *file_status, *file_t;
#endif

	typeTime tic, toc;
	typeRNum *CPUtimeVec;
	typeRNum CPUtime = 0;
	typeRNum rwsReferenceIntegration[2 * NX];

	/********* Parameter definition *********/
	/* Initial values and setpoints of the states, inputs, parameters, penalties and Lagrangian mmultipliers, setpoints for the states and inputs */
	ctypeRNum x0[NX] = { 1.8750,    0,   0,    3.7500,    0,    0,    5.6250,    0,    0,    7.5000,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0 };
	ctypeRNum xdes[NX] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	/* Initial values, setpoints and limits of the inputs */
	ctypeRNum u0[NU] = { 0, 0, 0 };
	ctypeRNum udes[NU] = { 0, 0, 0 };
	ctypeRNum umax[NU] = { 1.0, 1.0, 1.0 };
	ctypeRNum umin[NU] = { -1.0, -1.0, -1.0 };


	/* Time variables */
	ctypeRNum Thor = 8;  /* Prediction horizon */

	ctypeRNum dt = (typeRNum)0.1; /* Sampling time */
	typeRNum t = (typeRNum)0.0;   /* time at the current sampling step */

	/********* Option definition *********/
	/* Basic algorithmic options */
	ctypeInt MaxGradIter = (typeInt)5;        /* Maximum number of gradient Iterations */

	/* Cost integration */
	const char* TerminalCost = "off";

	/* System integration */
	const char* Integrator = "ruku45";
	ctypeRNum IntegratorRelTol = (typeRNum)1e-3;
	ctypeRNum IntegratorAbsTol = (typeRNum)1e-4;

	/********* userparam *********/
	typeRNum pCost[4] = { 25, (typeRNum)2.5, (typeRNum)0.1, 10 };
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

	grampc_setparam_real(grampc, "Thor", Thor);

	grampc_setparam_real(grampc, "dt", dt);
	grampc_setparam_real(grampc, "t0", t);

	/********* Option definition *********/
	grampc_setopt_int(grampc, "MaxGradIter", MaxGradIter);

	grampc_setopt_string(grampc, "TerminalCost", TerminalCost);

	grampc_setopt_string(grampc, "Integrator", Integrator);
	grampc_setopt_real(grampc, "IntegratorRelTol", IntegratorRelTol);
	grampc_setopt_real(grampc, "IntegratorAbsTol", IntegratorAbsTol);

#ifdef PRINTRES
	openFile(&file_x, "resNLC_4/xvec.txt");
	openFile(&file_u, "resNLC_4/uvec.txt");
	openFile(&file_p, "resNLC_4/pvec.txt");
	openFile(&file_T, "resNLC_4/Thorvec.txt");
	openFile(&file_J, "resNLC_4/Jvec.txt");
	openFile(&file_Ncfct, "resNLC_4/Ncfctvec.txt");
	openFile(&file_Npen, "resNLC_4/Npenvec.txt");
	openFile(&file_iter, "resNLC_4/itervec.txt");
	openFile(&file_status, "resNLC_4/status.txt");
	openFile(&file_t, "resNLC_4/tvec.txt");
#endif

	MaxSimIter = (int)(Tsim / grampc->param->dt);
	CPUtimeVec = (typeRNum*)calloc(MaxSimIter + 1, sizeof(*CPUtimeVec));
	if (CPUtimeVec == NULL) {
		printError("Allocating memory for computationtime measurement failed");
	}

	grampc_printopt(grampc);
	grampc_printparam(grampc);

	printf("MPC running ...\n");
	for (iMPC = 0; iMPC <= MaxSimIter; iMPC++) {
		/* run grampc */
		timer_now(&tic);
		grampc_run(grampc);
		timer_now(&toc);
		CPUtimeVec[iMPC] = timer_diff_ms(&tic, &toc);

		/* check solver status */
		if (grampc->sol->status > 0) {
			if (grampc_printstatus(grampc->sol->status, STATUS_LEVEL_ERROR)) {
				myPrint("at iteration %i:\n -----\n", iMPC);
			}
		}

		/* reference integration of the system  via heun scheme since grampc->sol->xnext is only an interpolated value */
		ffct(rwsReferenceIntegration, t, grampc->param->x0, grampc->sol->unext, grampc->sol->pnext, grampc->param, grampc->userparam);
		for (i = 0; i < NX; i++) {
			grampc->sol->xnext[i] = grampc->param->x0[i] + dt * rwsReferenceIntegration[i];
		}
		ffct(rwsReferenceIntegration + NX, t + dt, grampc->sol->xnext, grampc->sol->unext, grampc->sol->pnext, grampc->param, grampc->userparam);
		for (i = 0; i < NX; i++) {
			grampc->sol->xnext[i] = grampc->param->x0[i] + dt * (rwsReferenceIntegration[i] + rwsReferenceIntegration[i + NX]) / 2;
		}

		/* update state and time */
		t = t + dt;
		grampc_setparam_real_vector(grampc, "x0", grampc->sol->xnext);


#ifdef PRINTRES
		printNumVector2File(file_x, grampc->sol->xnext, NX);
		printNumVector2File(file_u, grampc->sol->unext, NU);
		/*printNumVector2File(file_p, grampc->sol->pnext, NP);*/
		printNumVector2File(file_T, &grampc->sol->Tnext, 1);
		printNumVector2File(file_J, grampc->sol->J, 2);
		printNumVector2File(file_Ncfct, &grampc->sol->cfct, 1);
		printNumVector2File(file_Npen, &grampc->sol->pen, 1);
		printIntVector2File(file_iter, grampc->sol->iter, grampc->opt->MaxMultIter);
		printIntVector2File(file_status, &grampc->sol->status, 1);
		printNumVector2File(file_t, &t, 1);
#endif
	}

#ifdef PRINTRES
	fclose(file_x);
	fclose(file_u);
	fclose(file_p);
	fclose(file_T);
	fclose(file_J);
	fclose(file_Ncfct);
	fclose(file_Npen);
	fclose(file_iter);
	fclose(file_status);
	fclose(file_t);
#endif

	for (i = 0; i <= MaxSimIter; i++) {
		CPUtime = CPUtime + CPUtimeVec[i];
	}
	CPUtime = CPUtime / (MaxSimIter + 1);

	grampc_free(&grampc);
	free(CPUtimeVec);

	printf("MPC finished. Average computation time: %.3f ms.\n", CPUtime);

	return 0;
}
