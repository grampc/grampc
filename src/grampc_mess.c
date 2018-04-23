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


#include "grampc_mess.h"


void grampc_error(const char *mess)
{
	printError(mess);
}

void grampc_error_addstring(const char *mess, const char *addstring)
{
	printErrorAddString(mess, addstring);
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
