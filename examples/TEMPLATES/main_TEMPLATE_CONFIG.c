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
    TYPE_GRAMPC_POINTER(grampc);
	typeInt iMPC, MaxSimIter, i;
	ctypeRNum Tsim = 20.0;

#ifdef PRINTRES
	FILE *file_x, *file_u, *file_p, *file_T, *file_J, *file_Ncfct, *file_Npen, *file_iter, *file_status, *file_t;
#endif

	typeTime tic, toc;
	typeRNum *CPUtimeVec;
	typeRNum CPUtime = 0;

	/********* userparam *********/
	typeRNum pCost[4] = { 0, 0, 0, 0 };
	typeUSERPARAM *userparam = pCost;

	/********* grampc init *********/
	grampc_init(&grampc, userparam);

	/********* Get and set options and parameter from configuration file *********/
	const char *fileName = "configFile_TEMPLATE.cfg";
	grampc_get_config_from_file(grampc, fileName);

	/********* Time Variables *********/
	typeRNum t = grampc->param->t0;	     /* time at the current sampling step */

	/********* Set additional options *********/
	ctypeRNum IntegratorMinStepSize = EPS;
	ctypeInt FlagsRodas[8] = { 0, 0, 0, grampc->param->Nx, grampc->param->Nx, 0, grampc->param->Nx, grampc->param->Nx };
	grampc_setopt_real(grampc, "IntegratorMinStepSize", IntegratorMinStepSize);
	grampc_setopt_int_vector(grampc, "FlagsRodas", FlagsRodas);


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
		timer_now(&tic);
		grampc_run(grampc);
		timer_now(&toc);
		CPUtimeVec[iMPC] = timer_diff_ms(&tic, &toc);
        
		grampc_setparam_real_vector(grampc, "x0", grampc->sol->xnext);
		t = t + grampc->param->dt;

		/* check solver status */
		if (grampc->sol->status > 0) {
			if (grampc_printstatus(grampc->sol->status, STATUS_LEVEL_ERROR)) {
				myPrint("at iteration %i:\n -----\n", iMPC);
			}
		}

#ifdef PRINTRES
		printNumVector2File(file_x, grampc->sol->xnext, grampc->param->Nx);
		printNumVector2File(file_u, grampc->sol->unext, grampc->param->Nu);
		/*printNumVector2File(file_p, grampc->sol->pnext, grampc->param->Np);*/
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
