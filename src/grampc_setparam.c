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



#include "grampc_setparam.h"


void grampc_setparam_real(const typeGRAMPC *grampc, const typeChar *paramName, ctypeRNum paramValue)
{
	/* Set prediction horizon */
	if (!strcmp(paramName, "Thor")) {
		grampc->param->Thor = paramValue;
		if (grampc->param->Thor <= 0) {
			grampc_error_addstring(INVALID_PARAM_VALUE, paramName);
		}
		init_rws_time(grampc);
	}
	/* Tmax */
	else if (!strcmp(paramName, "Tmax")) {
		grampc->param->Tmax = paramValue;
		if (grampc->param->Tmax <= 0) {
			grampc_error_addstring(INVALID_PARAM_VALUE, paramName);
		}
	}
	/* Tmin */
	else if (!strcmp(paramName, "Tmin")) {
		grampc->param->Tmin = paramValue;
		if (grampc->param->Tmin <= 0) {
			grampc_error_addstring(INVALID_PARAM_VALUE, paramName);
		}
	}
	/* Set sampling time */
	else if (!strcmp(paramName, "dt")) {
		grampc->param->dt = paramValue;
		if (grampc->param->dt <= 0) {
			grampc_error_addstring(INVALID_PARAM_VALUE, paramName);
		}
		init_rws_time(grampc);
	}
	/* Set current time */
	else if (!strcmp(paramName, "t0")) {
		grampc->param->t0 = paramValue;
	}
	/* Undefined parameter */
	else {
		grampc_error_addstring(INVALID_PARAM_NAME, paramName);
	}
}

void grampc_setparam_real_vector(const typeGRAMPC *grampc, const typeChar *paramName, ctypeRNum *paramValue)
{
	/* x0 */
	if (!strcmp(paramName, "x0")) {
		MatCopy(grampc->param->x0, paramValue, 1, grampc->param->Nx);
	}
	/* xdes */
	else if (!strcmp(paramName, "xdes")) {
		MatCopy(grampc->param->xdes, paramValue, 1, grampc->param->Nx);
	}
	/* u0 */
	else if (!strcmp(paramName, "u0")) {
		MatCopy(grampc->param->u0, paramValue, 1, grampc->param->Nu);
		init_rws_controls(grampc);
	}
	/* udes */
	else if (!strcmp(paramName, "udes")) {
		MatCopy(grampc->param->udes, paramValue, 1, grampc->param->Nu);
	}
	/* umax */
	else if (!strcmp(paramName, "umax")) {
		MatCopy(grampc->param->umax, paramValue, 1, grampc->param->Nu);
		check_ControlLimits(grampc);
	}
	/* umin */
	else if (!strcmp(paramName, "umin")) {
		MatCopy(grampc->param->umin, paramValue, 1, grampc->param->Nu);
		check_ControlLimits(grampc);
	}
	/* p0 */
	else if (!strcmp(paramName, "p0")) {
		MatCopy(grampc->param->p0, paramValue, 1, grampc->param->Np);
		init_rws_parameters(grampc);
	}
	/* pmax */
	else if (!strcmp(paramName, "pmax")) {
		MatCopy(grampc->param->pmax, paramValue, 1, grampc->param->Np);
	}
	/* pmin */
	else if (!strcmp(paramName, "pmin")) {
		MatCopy(grampc->param->pmin, paramValue, 1, grampc->param->Np);
	}
	/* Undefined parameter */
	else {
		grampc_error_addstring(INVALID_PARAM_NAME, paramName);
	}
}

void grampc_printparam(const typeGRAMPC *grampc)
{
	myPrint("%s", "-- MPC PARAMETERS --\n");
	myPrint("     Nx: %d\n", grampc->param->Nx);
	myPrint("     Nu: %d\n", grampc->param->Nu);
	myPrint("     Np: %d\n", grampc->param->Np);
	myPrint("     Nc: %d\n", grampc->param->Nc);
	myPrint("     Ng: %d\n", grampc->param->Ng);
	myPrint("     Nh: %d\n", grampc->param->Nh);
	myPrint("    NgT: %d\n", grampc->param->NgT);
	myPrint("    NhT: %d\n", grampc->param->NhT);

	print_vector("     x0: ", grampc->param->x0, grampc->param->Nx);
	print_vector("   xdes: ", grampc->param->xdes, grampc->param->Nx);

	print_vector("     u0: ", grampc->param->u0, grampc->param->Nu);
	print_vector("   udes: ", grampc->param->udes, grampc->param->Nu);
	print_vector("   umax: ", grampc->param->umax, grampc->param->Nu);
	print_vector("   umin: ", grampc->param->umin, grampc->param->Nu);

	print_vector("     p0: ", grampc->param->p0, grampc->param->Np);
	print_vector("   pmax: ", grampc->param->pmax, grampc->param->Np);
	print_vector("   pmin: ", grampc->param->pmin, grampc->param->Np);

	myPrint("   Thor: %.2f\n", grampc->param->Thor);
	myPrint("   Tmax: %.2f\n", grampc->param->Tmax);
	myPrint("   Tmin: %.2f\n", grampc->param->Tmin);

	myPrint("     dt: %.4f\n", grampc->param->dt);
	myPrint("     t0: %.4f\n", grampc->param->t0);
}
