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

#define NX	11
#define NU	1
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

int main(void)
{
    TYPE_GRAMPC_POINTER(grampc)
	typeInt iMPC, MaxSimIter, i;
	ctypeRNum Tsim = 1.0;

#ifdef PRINTRES
	FILE *file_x, *file_u, *file_p, *file_T, *file_J, *file_Ncfct, *file_Npen, *file_iter, *file_status, *file_t;
#endif

	clock_t tic, toc;
	typeRNum *CPUtimeVec;
	typeRNum CPUtime = 0;

	/********* Parameter definition *********/
	/* Initial values and setpoints of the states, inputs, parameters, penalties and Lagrangian multipliers, setpoints for the states and inputs */
	ctypeRNum x0[NX] = { 1.0000,    0.9970,    0.9876,    0.9717,    0.9491,    0.9195,    0.8829,    0.8390,    0.7878,    0.7294,    0.6637 };
	ctypeRNum xdes[NX] = { 2.0000,    1.9929,    1.9711,    1.9342,    1.8816,    1.8130,    1.7283,    1.6273,    1.5101,    1.3768,    1.2276 };


	/* Initial values, setpoints and limits of the inputs */
	ctypeRNum u0[NU] = { -0.6935 };
	ctypeRNum udes[NU] = { -1.5695 };

	ctypeRNum umax[NU] = { 2.0 };
	ctypeRNum umin[NU] = { -2.0 };

	/* Time variables */
	ctypeRNum Thor = 0.4;       /* Prediction horizon */
	ctypeRNum dt = 0.005;       /* Sampling time */
	typeRNum t = (typeRNum)0.0; /* time at the current sampling step */

	/********* Option definition *********/
	/* Basic algorithmic options */
	ctypeInt MaxGradIter = 2;   /* Maximum number of gradient iterations */
	ctypeRNum Nhor = 60;

	/* System integration */
	const char* Integrator = "rodas";
	ctypeRNum IntegratorRelTol = 1e-3;
	ctypeRNum IntegratorAbsTol = 1e-4;
	typeInt IntegratorMaxSteps = 50000;

	/*rodas specific options*/
	ctypeInt IFCN = 0;    /* 0 --> right hand side independent of time t  */
	ctypeInt IDFX = 0;    /* 0 --> DF/Dt is numerically computed */

	ctypeInt IJAC = 1;    /* 1(0) -> analytical (numerical) jacobian (partial derivatives of right hand side w.r.t. state) */
	ctypeInt MLJAC = 1;   /* no. of lower diagonals of jacobian ~	value must be between 1 and NX */
	ctypeInt MUJAC = 1;   /* no. of upper diagonals of jacobian ~	value must be between 1 and NX */

	ctypeInt IMAS = 1;    /* 1 --> mass matrix is supplied */
	ctypeInt MLMAS = 1;   /* no. of lower diagonals of mass matrix */
	ctypeInt MUMAS = 1;   /* no. of upper diagonals of mass matrix */

	ctypeInt FlagsRodas[8] = { IFCN, IDFX, IJAC, IMAS, MLJAC, MUJAC, MLMAS, MUMAS };

	/* Line search */
	ctypeRNum LineSearchMax = 2.0;
	ctypeRNum LineSearchInit = 2e-2;

	/********* userparam *********/
	typeRNum pCost[23] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
		10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
		0.01 };
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

	/********* Option definition *********/
	grampc_setopt_int(grampc, "Nhor", Nhor);
	grampc_setopt_int(grampc, "MaxGradIter", MaxGradIter);

	grampc_setopt_string(grampc, "Integrator", Integrator);
	grampc_setopt_real(grampc, "IntegratorRelTol", IntegratorRelTol);
	grampc_setopt_real(grampc, "IntegratorAbsTol", IntegratorAbsTol);
	grampc_setopt_int(grampc, "IntegratorMaxSteps", IntegratorMaxSteps);
	grampc_setopt_int_vector(grampc, "FlagsRodas", FlagsRodas);

	grampc_setopt_real(grampc, "LineSearchMax", LineSearchMax);
	grampc_setopt_real(grampc, "LineSearchInit", LineSearchInit);


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

	grampc_printopt(grampc);
	grampc_printparam(grampc);

	MaxSimIter = (int)(Tsim / grampc->param->dt);
	CPUtimeVec = (typeRNum*)calloc(MaxSimIter + 1, sizeof(*CPUtimeVec));
	if (CPUtimeVec == NULL) {
		printError("Allocating memory for computationtime measurement failed");
	}

	printf("MPC running ...\n");
	for (iMPC = 0; iMPC <= MaxSimIter; iMPC++) {
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
