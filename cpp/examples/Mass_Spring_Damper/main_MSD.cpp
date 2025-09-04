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

#include <chrono>
#include <vector>
#include <fstream>
#include "grampc.hpp"
#include "problem_description_MSD.hpp"

#define PRINTRES

#ifdef PRINTRES
std::ofstream openFile(const char *name) {
    std::ofstream file (name);
    file.precision(5);
    return file;
}

/* Template to use write function for floats and ints */
template <typename T>
void printVector2File(std::ofstream& file, const T* val, ctypeInt size) {
	for (typeInt i = 0; i < size - 1; i++) {
        file << val[i] << ", ";
	}
    file << val[size - 1] << std::endl;
}
#endif

int main()
{
    /* Initialize problem description and parameters */
    typeInt NN = 5;
    typeInt Nx = NN*2;
    typeInt Nu = 2;
    ctypeRNum m = 1;
    ctypeRNum c = 1;
    ctypeRNum d = 0.2;

    std::vector<typeRNum> Q(Nx, 1.0);
    std::vector<typeRNum> R(Nu, 0.01);
    std::vector<typeRNum> P(Nx, 1.0);

    MassSpringDamper problem(NN, m, c, d, Q, R, P);

#ifdef PRINTRES
	std::ofstream xOut = openFile("xvec.txt");
	std::ofstream uOut = openFile("uvec.txt");
	std::ofstream TOut = openFile("Thorvec.txt");
	std::ofstream JOut = openFile("Jvec.txt");
	std::ofstream ConstraintsOut = openFile("Ncfctvec.txt");
	std::ofstream NpenOut = openFile("Npenvec.txt");
	std::ofstream iterOut = openFile("itervec.txt");
	std::ofstream statusOut = openFile("status.txt");
	std::ofstream tOut = openFile("tvec.txt");
#endif

    /* Reference simulation definitions */
	ctypeRNum Tsim = 12.0;

	/********* Parameter definition *********/
	/* Initial values and setpoints of the states, inputs, parameters, penalties and Lagrangian mmultipliers, setpoints for the states and inputs */
	std::vector<typeRNum> x0 { 1, 0, 0, 0, 1, 0, 0, 0, 0, 0 };
	std::vector<typeRNum> xdes { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	/* Initial values, setpoints and limits of the inputs */
	std::vector<typeRNum> u0 { 0.0, 0.0 };
	std::vector<typeRNum> udes { 0.0, 0.0 };
	std::vector<typeRNum> umax { 1.0, 1.0 };
	std::vector<typeRNum> umin { -1.0, -1.0 };

	/* Time variables */
	ctypeRNum Thor = 10;  /* Prediction horizon */
	ctypeRNum dt = (typeRNum)0.005;  /* Sampling time */

	/********* Option definition *********/
	/* Basic algorithmic options */
	ctypeInt MaxGradIter = 5;  /* Maximum number of gradient iterations */

    /* Initialize Grampc with the specific problem */
	grampc::Grampc grampc(&problem);

	/********* set parameters *********/
    grampc.setparam_real_vector("x0", x0.data());
	grampc.setparam_real_vector("xdes", xdes.data());

	grampc.setparam_real_vector("u0", u0.data());
	grampc.setparam_real_vector("udes", udes.data());
	grampc.setparam_real_vector("umax", umax.data());
	grampc.setparam_real_vector("umin", umin.data());

	grampc.setparam_real("Thor", Thor);

	grampc.setparam_real("dt", dt);

	/********* Option definition *********/
	grampc.setopt_int("MaxGradIter", MaxGradIter);

	grampc.printparam();
    grampc.printopt();
    
	typeRNum t = (typeRNum)0.0;
	typeInt MaxSimIter = (int)(Tsim / dt);
	std::vector<typeRNum> CPUtimeVec(MaxSimIter + 1);
	std::vector<typeRNum> rwsReferenceIntegration(2 * Nx);
	
	printf("MPC running ...\n");
	for (typeInt iMPC = 0; iMPC <= MaxSimIter; iMPC++) {
		/* run grampc */
		auto tic = std::chrono::system_clock::now();
		grampc.run();
		auto toc = std::chrono::system_clock::now();
		CPUtimeVec[iMPC] = (typeRNum)std::chrono::duration_cast<std::chrono::nanoseconds>(toc - tic).count() * 1e-6;

        const typeGRAMPCsol *sol = grampc.getSolution();
        const typeGRAMPCparam *param = grampc.getParameters();
        const typeGRAMPCopt *opt = grampc.getOptions();

		/* check solver status */
		if (sol->status > 0) {
			if (grampc.printstatus(sol->status, STATUS_LEVEL_ERROR)) {
				myPrint("at iteration %i:\n -----\n", iMPC);
			}
		}

		/* reference integration of the system  via heun scheme since grampc->sol->xnext is only an interpolated value */
		problem.ffct(rwsReferenceIntegration.data(), t, param->x0, sol->unext, sol->pnext, param);
		for (typeInt i = 0; i < Nx; i++) {
			sol->xnext[i] = param->x0[i] + dt * rwsReferenceIntegration[i];
		}
		problem.ffct(rwsReferenceIntegration.data() + Nx, t + dt, sol->xnext, sol->unext, sol->pnext, param);
		for (typeInt i = 0; i < Nx; i++) {
			sol->xnext[i] = param->x0[i] + dt * (rwsReferenceIntegration[i] + rwsReferenceIntegration[i + Nx]) / 2;
		}

		/* update state and time */
		t = t + dt;
		grampc.setparam_real_vector("x0", sol->xnext);


#ifdef PRINTRES
		printVector2File<typeRNum>(xOut, sol->xnext, Nx);
		printVector2File<typeRNum>(uOut, sol->unext, Nu);
		printVector2File<typeRNum>(TOut, &sol->Tnext, 1);
		printVector2File<typeRNum>(JOut, sol->J, 2);
		printVector2File<typeRNum>(ConstraintsOut, &sol->cfct, 1);
		printVector2File<typeRNum>(NpenOut, &sol->pen, 1);
		printVector2File<typeInt>(iterOut, sol->iter, opt->MaxMultIter);
		printVector2File<typeInt>(statusOut, &sol->status, 1);
		printVector2File<typeRNum>(tOut, &t, 1);
#endif
	}

	typeRNum CPUtime = (typeRNum)0.0;
	for (typeInt i = 0; i <= MaxSimIter; i++) {
		CPUtime = CPUtime + CPUtimeVec[i] / (MaxSimIter + 1);
	}
	printf("MPC finished. Average computation time: %.3f ms.\n", CPUtime);

	return 0;
}
