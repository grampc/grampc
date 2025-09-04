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


#include "grampc_erk.h"


void intsysERK1(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum *t, ctypeRNum *x,
	ctypeRNum *u, ctypeRNum *p, ctypeRNum *dcdx, const typeGRAMPC *grampc, const typeSysPtr pfct)
{
	typeInt i, j;
	typeRNum h;
	typeRNum *k1 = grampc->rws->rwsGeneral + LWadjsys; /* size Nx */

	for (j = 0; j < Nint - 1; j++) {
		/* step size */
		h = t[pInt] - t[0];

		/* stage 1 */
		(*pfct)(k1, y, t, x, u, p, dcdx, grampc);

		/* update */
		for (i = 0; i < grampc->param->Nx; i++) {
			y[i + pInt * grampc->param->Nx] = y[i] + h * k1[i];
		}

		/* next pointers */
		t += pInt;
		x += pInt * grampc->param->Nx;
		u += pInt * grampc->param->Nu;
		y += pInt * grampc->param->Nx;
		dcdx += pInt * grampc->param->Nx;
	}
}


void intsysERK2(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum *t, ctypeRNum *x,
	ctypeRNum *u, ctypeRNum *p, ctypeRNum *dcdx, const typeGRAMPC *grampc, const typeSysPtr pfct)
{
	typeInt i, j;
	typeRNum h;
	typeRNum *k1 = grampc->rws->rwsGeneral + LWadjsys; /* size 3*Nx */
	typeRNum *k2 = k1 + grampc->param->Nx;
	typeRNum *ys = k2 + grampc->param->Nx;

	for (j = 0; j < Nint - 1; j++) {
		/* step size */
		h = t[pInt] - t[0];

		/* stage 1 */
		(*pfct)(k1, y, t, x, u, p, dcdx, grampc);

		/* stage 2 */
		for (i = 0; i < grampc->param->Nx; i++) {
			ys[i] = y[i] + h * k1[i];
		}
		(*pfct)(k2, ys, t + pInt, x + pInt * grampc->param->Nx, u + pInt * grampc->param->Nu, p, dcdx + pInt * grampc->param->Nx, grampc);

		/* update */
		for (i = 0; i < grampc->param->Nx; i++) {
			y[i + pInt * grampc->param->Nx] = y[i] + h * (k1[i] + k2[i]) / 2.0;
		}

		/* next pointers */
		t += pInt;
		x += pInt * grampc->param->Nx;
		u += pInt * grampc->param->Nu;
		y += pInt * grampc->param->Nx;
		dcdx += pInt * grampc->param->Nx;
	}
}


void intsysERK3(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum *t, ctypeRNum *x,
	ctypeRNum *u, ctypeRNum *p, ctypeRNum *dcdx, const typeGRAMPC *grampc, const typeSysPtr pfct)
{
	typeInt i, j;
	typeRNum h, ts;
	typeRNum *k1 = grampc->rws->rwsGeneral + LWadjsys; /* size 6*Nx+Nu */
	typeRNum *k2 = k1 + grampc->param->Nx;
	typeRNum *k3 = k2 + grampc->param->Nx;
	typeRNum *ys = k3 + grampc->param->Nx;
	typeRNum *xs = ys + grampc->param->Nx;
	typeRNum *cs = xs + grampc->param->Nx;
	typeRNum *us = cs + grampc->param->Nx;

	for (j = 0; j < Nint - 1; j++) {
		/* step size */
		h = t[pInt] - t[0];

		/* stage 1 */
		(*pfct)(k1, y, t, x, u, p, dcdx, grampc);

		/* stage 2 */
		ts = t[0] + h / 2;
		for (i = 0; i < grampc->param->Nx; i++) {
			ys[i] = y[i] + h * k1[i] / 2;
		}
		if (y != x) {
			/* interpolate x at sample point */
			for (i = 0; i < grampc->param->Nx; i++) {
				xs[i] = (x[i] + x[i + pInt * grampc->param->Nx]) / 2;
			}
			/* interpolate dcdx at sample point */
			for (i = 0; i < grampc->param->Nx; i++) {
				cs[i] = (dcdx[i] + dcdx[i + pInt * grampc->param->Nx]) / 2;
			}
		}
		/* interpolate u at sample point */
		for (i = 0; i < grampc->param->Nu; i++) {
			us[i] = (u[i] + u[i + pInt * grampc->param->Nu]) / 2;
		}
		(*pfct)(k2, ys, &ts, xs, us, p, cs, grampc);
		
		/* stage 3 */
		for (i = 0; i < grampc->param->Nx; i++) {
			ys[i] = y[i] + h * (-k1[i] + 2 * k2[i]);
		}
		(*pfct)(k3, ys, t + pInt, x + pInt * grampc->param->Nx, u + pInt * grampc->param->Nu, p, dcdx + pInt * grampc->param->Nx, grampc);

		/* update */
		for (i = 0; i < grampc->param->Nx; i++) {
			y[i + pInt * grampc->param->Nx] = y[i] + h * (k1[i] + 4 * k2[i] + k3[i]) / 6.0;
		}

		/* next pointers */
		t += pInt;
		x += pInt * grampc->param->Nx;
		u += pInt * grampc->param->Nu;
		y += pInt * grampc->param->Nx;
		dcdx += pInt * grampc->param->Nx;
	}
}


void intsysERK4(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum *t, ctypeRNum *x,
	ctypeRNum *u, ctypeRNum *p, ctypeRNum *dcdx, const typeGRAMPC *grampc, const typeSysPtr pfct)
{
	typeInt i, j;
	typeRNum h, ts;
	typeRNum *k1 = grampc->rws->rwsGeneral + LWadjsys; /* size 7*Nx+Nu */
	typeRNum *k2 = k1 + grampc->param->Nx;
	typeRNum *k3 = k2 + grampc->param->Nx;
	typeRNum *k4 = k3 + grampc->param->Nx;
	typeRNum *ys = k4 + grampc->param->Nx;
	typeRNum *xs = ys + grampc->param->Nx;
	typeRNum *cs = xs + grampc->param->Nx;
	typeRNum *us = cs + grampc->param->Nx;

	for (j = 0; j < Nint - 1; j++) {
		/* step size */
		h = t[pInt] - t[0];

		/* stage 1 */
		(*pfct)(k1, y, t, x, u, p, dcdx, grampc);

		/* stage 2 */
		ts = t[0] + h / 2;
		for (i = 0; i < grampc->param->Nx; i++) {
			ys[i] = y[i] + h * k1[i] / 2;
		}
		if (y != x) {
			/* interpolate x at sample point */
			for (i = 0; i < grampc->param->Nx; i++) {
				xs[i] = (x[i] + x[i + pInt * grampc->param->Nx]) / 2;
			}
			/* interpolate dcdx at sample point */
			for (i = 0; i < grampc->param->Nx; i++) {
				cs[i] = (dcdx[i] + dcdx[i + pInt * grampc->param->Nx]) / 2;
			}
		}
		/* interpolate u at sample point */
		for (i = 0; i < grampc->param->Nu; i++) {
			us[i] = (u[i] + u[i + pInt * grampc->param->Nu]) / 2;
		}
		(*pfct)(k2, ys, &ts, xs, us, p, cs, grampc);
		
		/* stage 3 */
		for (i = 0; i < grampc->param->Nx; i++) {
			ys[i] = y[i] + h * k2[i] / 2;
		}
		(*pfct)(k3, ys, &ts, xs, us, p, cs, grampc);

		/* stage 4 */
		for (i = 0; i < grampc->param->Nx; i++) {
			ys[i] = y[i] + h * k3[i];
		}
		(*pfct)(k4, ys, t + pInt, x + pInt * grampc->param->Nx, u + pInt * grampc->param->Nu, p, dcdx + pInt * grampc->param->Nx, grampc);

		/* update */
		for (i = 0; i < grampc->param->Nx; i++) {
			y[i + pInt * grampc->param->Nx] = y[i] + h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6.0;
		}

		/* next pointers */
		t += pInt;
		x += pInt * grampc->param->Nx;
		u += pInt * grampc->param->Nu;
		y += pInt * grampc->param->Nx;
		dcdx += pInt * grampc->param->Nx;
	}
}
