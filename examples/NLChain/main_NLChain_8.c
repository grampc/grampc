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

#define N 8 /* number of chain elements */

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
	typeInt iMPC, MaxSimIter, i;
	ctypeRNum Tsim = 20.0;

#ifdef PRINTRES
	FILE *file_x, *file_u, *file_p, *file_T, *file_J, *file_Ncfct, *file_Npen, *file_iter, *file_status, *file_t;
#endif

	clock_t tic, toc;
	typeRNum *CPUtimeVec;
	typeRNum CPUtime = 0;

	/********* Parameter definition *********/
	/* Initial values and setpoints of the states, inputs, parameters, penalties and Lagrangian mmultipliers, setpoints for the states and inputs */
	ctypeRNum x0[NX] = { 0.9375,    0,    0,    1.8750,    0,    0,    2.8125,    0,    0,    3.7500,    0,    0,    4.6875,    0,    0 ,
		5.6250,    0,    0,    6.5625,    0,    0,    7.5000,    0,    0,    0 ,   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0 };
	ctypeRNum xdes[NX] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 0, 0 };

	/* Initial values, setpoints and limits of the inputs */
	ctypeRNum u0[NU] = { 0, 0, 0 };
	ctypeRNum udes[NU] = { 0, 0, 0 };
	ctypeRNum umax[NU] = { 1.0, 1.0, 1.0 };
	ctypeRNum umin[NU] = { -1.0, -1.0, -1.0 };


	/* Time variables */
	ctypeRNum Thor = 8;  /* Prediction horizon */

	ctypeRNum dt = (typeRNum)0.1;  /* Sampling time */
	typeRNum t = (typeRNum)0.0;    /* time at the current sampling step */

	/********* Option definition *********/
	/* Basic algorithmic options */
	ctypeInt MaxGradIter = (typeInt)5;  /* Maximum number of gradient Iterations */

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
	openFile(&file_x, "resNLC_8/xvec.txt");
	openFile(&file_u, "resNLC_8/uvec.txt");
	openFile(&file_p, "resNLC_8/pvec.txt");
	openFile(&file_T, "resNLC_8/Thorvec.txt");
	openFile(&file_J, "resNLC_8/Jvec.txt");
	openFile(&file_Ncfct, "resNLC_8/Ncfctvec.txt");
	openFile(&file_Npen, "resNLC_8/Npenvec.txt");
	openFile(&file_iter, "resNLC_8/itervec.txt");
	openFile(&file_status, "resNLC_8/status.txt");
	openFile(&file_t, "resNLC_8/tvec.txt");
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
		tic = clock();
		grampc_run(grampc);
		toc = clock();
		CPUtimeVec[iMPC] = (typeRNum)((toc - tic) * 1000 / CLOCKS_PER_SEC);

		/* check solver status */
		if (grampc->sol->status > 0) {
			if (grampc_printstatus(grampc->sol->status, STATUS_LEVEL_ERROR)) {
				myPrint("at iteration %i:\n -----\n", iMPC);
			}
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
		CPUtime = CPUtime + CPUtimeVec[i] / (MaxSimIter + 1);
	}

	grampc_free(&grampc);
	free(CPUtimeVec);

	printf("MPC finished. Average computation time: %.3f ms.\n", CPUtime);

	return 0;
}
