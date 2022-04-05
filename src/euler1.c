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



#include "euler1.h"

void intsysEuler(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum * t, ctypeRNum * x,
	ctypeRNum *u, ctypeRNum *p_, const typeGRAMPC *grampc, const typeffctPtr  pfct)
{
	typeInt i, j;
	typeRNum h;
	ctypeRNum *dcdx = grampc->rws->dcdx + grampc->param->Nx * (grampc->opt->Nhor - 1);
	typeRNum *s = grampc->rws->rwsGeneral + LWadjsys; /* size Nx */

	for (j = 0; j < Nint - 1; j++) {
		if (j > 0) {
			t += pInt;
			x += pInt * grampc->param->Nx;
			u += pInt * grampc->param->Nu;
			y += pInt * grampc->param->Nx;

			dcdx += (-1)*grampc->param->Nx; /* only used for integrating the adjoint system */
		}

		(*pfct)(s, y, t, x, u, p_, dcdx, grampc);

		h = t[pInt] - t[0];
		for (i = 0; i < grampc->param->Nx; i++) {
			y[i + pInt * grampc->param->Nx] = y[i] + h * s[i];
		}
	}
}
