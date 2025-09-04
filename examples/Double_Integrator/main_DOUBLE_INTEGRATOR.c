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

#define NX	2
#define NU	1
#define NC  2
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
		fprintf(file, "%.5e, ", val[i]);
	}
	fprintf(file, "%.5e\n", val[size - 1]); /* new line */
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
	ctypeRNum Tsim = 10.0;

#ifdef PRINTRES
	FILE *file_x, *file_u, *file_p, *file_T, *file_J, *file_Ncfct, *file_Npen, *file_iter, *file_status, *file_t;
#endif

	typeTime tic, toc;
	typeRNum *CPUtimeVec;
	typeRNum CPUtime = 0;
	typeRNum rwsReferenceIntegration[2 * NX];

	/********* userparam *********/
	typeRNum pCost[2] = { 0.1, 1 };
	typeUSERPARAM *userparam = pCost;

	/********* grampc init *********/
	grampc_init(&grampc, userparam);

	int counter = 0;

	/********* Get and set options and parameter from configuration file *********/
	const char *fileName = "config_Double_Integrator.cfg";
	grampc_get_config_from_file(grampc, fileName);

#ifdef PRINTRES
	openFile(&file_x, "resMPC/xvec.txt");
	openFile(&file_u, "resMPC/uvec.txt");
	openFile(&file_p, "resMPC/pvec.txt");
	openFile(&file_T, "resMPC/Thorvec.txt");
	openFile(&file_J, "resMPC/Jvec.txt");
	openFile(&file_Ncfct, "resMPC/Ncfctvec.txt");
	openFile(&file_Npen, "resMPC/Npenvec.txt");
	openFile(&file_iter, "resMPC/itervec.txt");
	openFile(&file_status, "resMPC/status.txt");
	openFile(&file_t, "resMPC/tvec.txt");
#endif

	grampc_printopt(grampc);
	grampc_printparam(grampc);


	/* Get time variables */
	ctypeRNum dt = grampc->param->dt;    /* Sampling time */
	typeRNum t = grampc->param->t0;	     /* time at the current sampling step */

	MaxSimIter = (int)(Tsim / grampc->param->dt);
	CPUtimeVec = (typeRNum*)calloc(MaxSimIter + 1, sizeof(*CPUtimeVec));
	if (CPUtimeVec == NULL) {
		printError("Allocating memory for computationtime measurement failed");
	}

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

		/* reference integration of the system via heun scheme since grampc->sol->xnext is only an interpolated value */
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

		/*  already minimal time -> deactivate time optimization */
		if ((grampc->sol->Tnext - grampc->param->dt) <= grampc->param->Tmin) {
			grampc_setopt_string(grampc, "OptimTime", "off");
			counter++;
		}

		if (!((grampc->sol->Tnext - grampc->param->dt) >= grampc->param->Tmin) && (counter >= 2))
		{
			break;
		}

#ifdef PRINTRES
		printNumVector2File(file_x, grampc->sol->xnext, NX);
		printNumVector2File(file_u, grampc->sol->unext, NU);
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
