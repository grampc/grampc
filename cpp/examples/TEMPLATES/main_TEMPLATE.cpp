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
#include "problem_description_TEMPLATE.hpp"

#define NX 1
#define NU 1
#define NC 1
#define NP 1

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
    TemplateProblemDescription problem = TemplateProblemDescription();

#ifdef PRINTRES
	std::ofstream xOut = openFile("xvec.txt");
	std::ofstream uOut = openFile("uvec.txt");
	std::ofstream pOut = openFile("vec.txt");
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
	ctypeRNum x0[NX] = {	};
	ctypeRNum xdes[NX] = { };

	/* Initial values, setpoints and limits of the inputs */
	ctypeRNum u0[NU] = { };
	ctypeRNum udes[NU] = {	};
	ctypeRNum umax[NU] = {	};
	ctypeRNum umin[NU] = {	};

	/* Initial values and limits of the parameters */
	ctypeRNum p0[NP] = {	};
	ctypeRNum pmax[NP] = {	};
	ctypeRNum pmin[NP] = {	};

	/* Time variables */
	ctypeRNum Thor = -10;	/* Prediction horizon */
	ctypeRNum Tmax = (typeRNum)1e8;
	ctypeRNum Tmin = (typeRNum)1e-8;

	ctypeRNum dt = (typeRNum)-1;  /* Sampling time */

	/********* Option definition *********/
	/* Basic algorithmic options */
	ctypeInt Nhor = (typeInt)30;        /* Number of steps for the system integration */
	ctypeInt MaxGradIter = (typeInt)2;  /* Maximum number of gradient iterations */
	ctypeInt MaxMultIter = (typeInt)1;  /* Maximum number of augmented Lagrangian iterations */
	const char* ShiftControl = "on";

	/* Cost integration */
	const char* IntegralCost = "on";
	const char* TerminalCost = "on";
	const char* IntegratorCost = "trapezoidal";

	/* System integration */
	const char* Integrator = "erk2";
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
	ctypeRNum LineSearchAdaptAbsTol = (typeRNum)1e-6;
	ctypeRNum LineSearchAdaptFactor = (typeRNum)3.0 / 2.0;
	ctypeRNum LineSearchIntervalTol = (typeRNum)1e-1;
	ctypeRNum LineSearchIntervalFactor = (typeRNum)0.85;

	/* Input and or parameter optimization	*/
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
	const char* ConvergenceCheck = "off";
	ctypeRNum ConvergenceGradientRelTol = (typeRNum)1e-6;

    /* Initialize Grampc with the specific problem */
	grampc::Grampc grampc = grampc::Grampc(&problem);

	/********* set parameters *********/
	grampc.setparam_real_vector("x0", x0);
	grampc.setparam_real_vector("xdes", xdes);

	grampc.setparam_real_vector("u0", u0);
	grampc.setparam_real_vector("udes", udes);
	grampc.setparam_real_vector("umax", umax);
	grampc.setparam_real_vector("umin", umin);

	grampc.setparam_real_vector("p0", p0);
	grampc.setparam_real_vector("pmax", pmax);
	grampc.setparam_real_vector("pmin", pmin);

	grampc.setparam_real("Thor", Thor);
	grampc.setparam_real("Tmax", Tmax);
	grampc.setparam_real("Tmin", Tmin);

	grampc.setparam_real("dt", dt);

	/********* Option definition *********/
	grampc.setopt_int("Nhor", Nhor);
	grampc.setopt_int("MaxGradIter", MaxGradIter);
	grampc.setopt_int("MaxMultIter", MaxMultIter);
	grampc.setopt_string("ShiftControl", ShiftControl);

	grampc.setopt_string("IntegralCost", IntegralCost);
	grampc.setopt_string("TerminalCost", TerminalCost);
	grampc.setopt_string("IntegratorCost", IntegratorCost);

	grampc.setopt_string("Integrator", Integrator);
	grampc.setopt_real("IntegratorRelTol", IntegratorRelTol);
	grampc.setopt_real("IntegratorAbsTol", IntegratorAbsTol);
	grampc.setopt_real("IntegratorMinStepSize", IntegratorMinStepSize);
	grampc.setopt_int("IntegratorMaxSteps", IntegratorMaxSteps);
	grampc.setopt_int_vector("FlagsRodas", FlagsRodas);

	grampc.setopt_string("LineSearchType", LineSearchType);
	grampc.setopt_string("LineSearchExpAutoFallback", LineSearchExpAutoFallback);
	grampc.setopt_real("LineSearchMax", LineSearchMax);
	grampc.setopt_real("LineSearchMin", LineSearchMin);
	grampc.setopt_real("LineSearchInit", LineSearchInit);
	grampc.setopt_real("LineSearchIntervalFactor", LineSearchIntervalFactor);
	grampc.setopt_real("LineSearchAdaptFactor", LineSearchAdaptFactor);
	grampc.setopt_real("LineSearchIntervalTol", LineSearchIntervalTol);

	grampc.setopt_string("OptimControl", OptimControl);
	grampc.setopt_string("OptimParam", OptimParam);
	grampc.setopt_real("OptimParamLineSearchFactor", OptimParamLineSearchFactor);
	grampc.setopt_string("OptimTime", OptimTime);
	grampc.setopt_real("OptimTimeLineSearchFactor", OptimTimeLineSearchFactor);

	grampc.setopt_string("ScaleProblem", ScaleProblem);
	grampc.setopt_real_vector("xScale", xScale);
	grampc.setopt_real_vector("xOffset", xOffset);
	grampc.setopt_real_vector("uScale", uScale);
	grampc.setopt_real_vector("uOffset", uOffset);
	grampc.setopt_real_vector("pScale", pScale);
	grampc.setopt_real_vector("pOffset", pOffset);
	grampc.setopt_real("TScale", TScale);
	grampc.setopt_real("TOffset", TOffset);
	grampc.setopt_real("JScale", JScale);
	grampc.setopt_real_vector("cScale", cScale);

	grampc.setopt_string("EqualityConstraints", EqualityConstraints);
	grampc.setopt_string("InequalityConstraints", InequalityConstraints);
	grampc.setopt_string("TerminalEqualityConstraints", TerminalEqualityConstraints);
	grampc.setopt_string("TerminalInequalityConstraints", TerminalInequalityConstraints);
	grampc.setopt_string("ConstraintsHandling", ConstraintsHandling);
	grampc.setopt_real_vector("ConstraintsAbsTol", ConstraintsAbsTol);

	grampc.setopt_real("MultiplierMax", MultiplierMax);
	grampc.setopt_real("MultiplierDampingFactor", MultiplierDampingFactor);
	grampc.setopt_real("PenaltyMax", PenaltyMax);
	grampc.setopt_real("PenaltyMin", PenaltyMin);
	grampc.setopt_real("PenaltyIncreaseFactor", PenaltyIncreaseFactor);
	grampc.setopt_real("PenaltyDecreaseFactor", PenaltyDecreaseFactor);
	grampc.setopt_real("PenaltyIncreaseThreshold", PenaltyIncreaseThreshold);
	grampc.setopt_real("AugLagUpdateGradientRelTol", AugLagUpdateGradientRelTol);

	grampc.setopt_string("ConvergenceCheck", ConvergenceCheck);
	grampc.setopt_real("ConvergenceGradientRelTol", ConvergenceGradientRelTol);

	grampc.printparam();
    grampc.printopt();
    
	typeRNum t = (typeRNum)0.0;
	typeInt MaxSimIter = (int)(Tsim / dt);
	std::vector<typeRNum> CPUtimeVec(MaxSimIter + 1);
	std::vector<typeRNum> rwsReferenceIntegration(2 * NX);

	printf("MPC running ...\n");
	for (typeInt iMPC = 0; iMPC <= MaxSimIter; iMPC++) {
		auto tic = std::chrono::system_clock::now();
		grampc.run();
		auto toc = std::chrono::system_clock::now();
		CPUtimeVec[iMPC] = (typeRNum)std::chrono::duration_cast<std::chrono::nanoseconds>(toc - tic).count() * 1e-6;
        
        const typeGRAMPCsol *sol = grampc.getSolution();
        const typeGRAMPCparam *param = grampc.getParameters();
        const typeGRAMPCopt *opt = grampc.getOptions();

		grampc.setparam_real_vector("x0", sol->xnext);
		t = t + param->dt;

		/* check solver status */
		if (sol->status > 0) {
			if (grampc.printstatus(sol->status, STATUS_LEVEL_ERROR)) {
				myPrint("at iteration %i:\n -----\n", iMPC);
			}
		}


#ifdef PRINTRES
		printVector2File<typeRNum>(xOut, sol->xnext, NX);
		printVector2File<typeRNum>(uOut, sol->unext, NU);
		printVector2File<typeRNum>(pOut, sol->unext, NP);
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
