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


#include "grampc_mess.h"


void grampc_error(const char *mess)
{
	printError(mess);
}

void grampc_error_addstring(const char *mess, const char *addstring)
{
	printErrorAddString(mess, addstring);
}

void grampc_warning_addstring(const char *mess, const char *addstring)
{
    printWarningAddString(mess, addstring);
}

void grampc_error_optname(const char *optname)
{
	printErrorAddString(INVALID_OPTION_NAME, optname);
}

void grampc_error_optvalue(const char *optname)
{
	printErrorAddString(INVALID_OPTION_VALUE, optname);
}

void grampc_error_paramname(const char *paramname)
{
	printErrorAddString(INVALID_PARAM_NAME, paramname);
}

void grampc_error_paramvalue(const char *paramname)
{
	printErrorAddString(INVALID_PARAM_VALUE, paramname);
}

void print_vector(const char *prefix, ctypeRNum *vector, ctypeInt size)
{
	typeInt i;
	if (vector == NULL) {
		myPrint("%s[]\n", prefix);
	}
	else if (size == 1) {
		myPrint("%s", prefix);
		myPrint("%.3f\n", vector[0]);
	}
	else {
		myPrint("%s[", prefix);
		for (i = 0; i < size - 1; i++) {
			myPrint("%.3f,", vector[i]);
		}
		myPrint("%.3f]\n", vector[size - 1]);
	}
}

typeInt grampc_printstatus(ctypeInt status, ctypeInt level)
{
	typeInt printed = 0;
	if (level & 1) {
		printed |= print_singleStatus(status, STATUS_INTEGRATOR_INPUT_NOT_CONSISTENT, STATUS_MSG_INTEGRATOR_INPUT_NOT_CONSISTENT);
		printed |= print_singleStatus(status, STATUS_INTEGRATOR_MAXSTEPS, STATUS_MSG_INTEGRATOR_MAXSTEPS);
		printed |= print_singleStatus(status, STATUS_INTEGRATOR_STEPS_TOO_SMALL, STATUS_MSG_INTEGRATOR_STEPS_TOO_SMALL);
		printed |= print_singleStatus(status, STATUS_INTEGRATOR_MATRIX_IS_SINGULAR, STATUS_MSG_INTEGRATOR_MATRIX_IS_SINGULAR);
		printed |= print_singleStatus(status, STATUS_INTEGRATOR_H_MIN, STATUS_MSG_INTEGRATOR_H_MIN);
	}
	if (level & 2) {
		printed |= print_singleStatus(status, STATUS_MULTIPLIER_MAX, STATUS_MSG_MULTIPLIER_MAX);
		printed |= print_singleStatus(status, STATUS_PENALTY_MAX, STATUS_MSG_PENALTY_MAX);
		printed |= print_singleStatus(status, STATUS_INFEASIBLE, STATUS_MSG_INFEASIBLE);
	}
	if (level & 4) {
		printed |= print_singleStatus(status, STATUS_GRADIENT_CONVERGED, STATUS_MSG_GRADIENT_CONVERGED);
		printed |= print_singleStatus(status, STATUS_CONSTRAINTS_CONVERGED, STATUS_MSG_CONSTRAINTS_CONVERGED);
		printed |= print_singleStatus(status, STATUS_LINESEARCH_INIT, STATUS_MSG_LINESEARCH_INIT);
	}
	if (level & 8) {
		printed |= print_singleStatus(status, STATUS_LINESEARCH_MIN, STATUS_MSG_LINESEARCH_MIN);
		printed |= print_singleStatus(status, STATUS_LINESEARCH_MAX, STATUS_MSG_LINESEARCH_MAX);
		printed |= print_singleStatus(status, STATUS_MULTIPLIER_UPDATE, STATUS_MSG_MULTIPLIER_UPDATE);
	}
	return printed;
}

typeInt print_singleStatus(ctypeInt status, ctypeInt statusmask, const typeChar*message) {
	if (status & statusmask) {
		myPrint("%s", message);
		return 1;
	}
	return 0;
}
