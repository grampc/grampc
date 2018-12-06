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


#include "grampc.h"
#include <time.h>

#define NX	3
#define NU	2
#define NC	1

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
	ctypeRNum Tsim = 2.9;

#ifdef PRINTRES
	FILE *file_x, *file_u, *file_p, *file_T, *file_J, *file_Ncfct, *file_Npen, *file_iter, *file_status, *file_t;
#endif

	clock_t tic, toc;
	typeRNum *CPUtimeVec;
	typeRNum CPUtime = 0;

	/********* Parameter definition *********/
	/* Initial values and setpoints of the states, inputs, parameters, penalties and Lagrangian mmultipliers, setpoints for the states and inputs */
  ctypeRNum x0[NX]   = { 1.0, 0.0, 1.0 };
  ctypeRNum xdes[NX] = { 0.0, 1.0, 1.0 };
	ctypeRNum xdes2[NX] = { 0.5, 0.5, 1.0 };

	/* Initial values, setpoints and limits of the inputs */
  ctypeRNum u0[NU]   = { -2.0, 2.0 };
  ctypeRNum udes[NU] = { 0.0, 0.0 };
  ctypeRNum umax[NU] = { 2.0, 2.0 };
  ctypeRNum umin[NU] = { -2.0, -2.0 };

	/* Time variables */
	ctypeRNum Thor = 0.4;	/* Prediction horizon */

	ctypeRNum dt = (typeRNum)0.01;  /* Sampling time */
	typeRNum t = (typeRNum)0.0;   /* time at the current sampling step */

	/********* Option definition *********/
	/* Basic algorithmic options */
	ctypeInt Nhor = (typeInt)30;        /* Number of steps for the system integration */
	ctypeInt MaxGradIter = (typeInt)10;  /* Maximum number of gradient iterations */
	ctypeInt MaxMultIter = (typeInt)3;  /* Maximum number of augmented Lagrangian iterations */

	/* Cost integration */
	const char* IntegralCost = "on";
	const char* TerminalCost = "off";

	/* System integration */
	const char* Integrator = "rodas";
	ctypeRNum IntegratorRelTol = (typeRNum)1e-4;
	ctypeRNum IntegratorAbsTol = (typeRNum)1e-5;
	ctypeInt FlagsRodas[8] = { 0, 0, 0, 1, NX, NX, 0, 0 };

	/* Line search */
	ctypeRNum LineSearchMax = (typeRNum)1e-1;

	/* Tye of constraints' consideration */
	ctypeRNum ConstraintsAbsTol[NC] = { 1e-4 };

	/* Multipliers & penalties */
	ctypeRNum PenaltyIncreaseFactor = (typeRNum)1.1;

	/********* userparam *********/
	typeRNum pCost[3] = { 500, 0, 1 };
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
	grampc_setparam_real(grampc, "dt", dt);
	grampc_setparam_real(grampc, "t0", t);

	/********* Option definition *********/
	grampc_setopt_int(grampc, "Nhor", Nhor);
	grampc_setopt_int(grampc, "MaxGradIter", MaxGradIter);
	grampc_setopt_int(grampc, "MaxMultIter", MaxMultIter);

	grampc_setopt_string(grampc, "IntegralCost", IntegralCost);
	grampc_setopt_string(grampc, "TerminalCost", TerminalCost);

	grampc_setopt_string(grampc, "Integrator", Integrator);
	grampc_setopt_real(grampc, "IntegratorRelTol", IntegratorRelTol);
	grampc_setopt_real(grampc, "IntegratorAbsTol", IntegratorAbsTol);
	grampc_setopt_int_vector(grampc, "FlagsRodas", FlagsRodas);

	grampc_setopt_real(grampc, "LineSearchMax", LineSearchMax);

	grampc_setopt_real_vector(grampc, "ConstraintsAbsTol", ConstraintsAbsTol);
	grampc_setopt_real(grampc, "PenaltyIncreaseFactor", PenaltyIncreaseFactor);

	/********* estimate and set PenaltyMin *********/
	grampc_estim_penmin(grampc, 1);

#ifdef PRINTRES
	openFile(&file_x, "res/xvec.txt");
	openFile(&file_u, "res/uvec.txt");
	openFile(&file_p, "res/pvec.txt");
	openFile(&file_T, "res/Thorvec.txt");
	openFile(&file_J, "res/Jvec.txt");
	openFile(&file_Ncfct, "res/Ncfctvec.txt");
	openFile(&file_Npen, "res/Npenvec.txt");
	openFile(&file_iter, "res/itervec.txt");
	openFile(&file_status, "res/status.txt");
	openFile(&file_t, "res/tvec.txt");
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

		if (2.0 <= t) {
			grampc_setparam_real_vector(grampc, "xdes", x0);
		}
		else if (1.0 <= t) {
			grampc_setparam_real_vector(grampc, "xdes", xdes2);
		}
		
		tic = clock();
		grampc_run(grampc);
		grampc_setparam_real_vector(grampc, "x0", grampc->sol->xnext);
		toc = clock();
		CPUtimeVec[iMPC] = (typeRNum)((toc - tic) * 1000 / CLOCKS_PER_SEC);
		t = t + grampc->param->dt;

		/* check solver status */
		if (grampc->sol->status > 0) {
			if (grampc_printstatus(grampc->sol->status, STATUS_LEVEL_ERROR)) {
				myPrint("at iteration %i:\n -----\n", iMPC);
			}
		}

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
