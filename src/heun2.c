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


#include "heun2.h"

void intsysHeun(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum *t, ctypeRNum *x,
	ctypeRNum *u, ctypeRNum *p_, const typeGRAMPC *grampc, const typeffctPtr pfct)
{
	typeInt i, j;
	typeRNum h, h2;
	ctypeRNum *dcdx = grampc->rws->dcdx + grampc->param->Nx * (grampc->opt->Nhor - 1);
	typeRNum *s1 = grampc->rws->rwsGeneral + LWadjsys; /* size 3*Nx */
	typeRNum *ys1 = s1 + grampc->param->Nx;
	typeRNum *s2 = ys1 + grampc->param->Nx;

	for (j = 0; j < Nint - 1; j++) {
		h = t[pInt] - t[0];
		h2 = h / 2;

		/* s1 */
		(*pfct)(s1, y, t, x, u, p_, dcdx, grampc);
		/* s2 */
		for (i = 0; i < grampc->param->Nx; i++) {
			ys1[i] = y[i] + h * s1[i];
		}
		(*pfct)(s2, ys1, t + pInt, x + pInt * grampc->param->Nx, u + pInt * grampc->param->Nu, p_, dcdx + (-1)*grampc->param->Nx, grampc);
		/* xnew */
		for (i = 0; i < grampc->param->Nx; i++) {
			y[i + pInt * grampc->param->Nx] = y[i] + h2 * (s1[i] + s2[i]);
		}

		/* next pointers */
		t += pInt;
		x += pInt * grampc->param->Nx;
		u += pInt * grampc->param->Nu;
		y += pInt * grampc->param->Nx;

		dcdx += (-1)*grampc->param->Nx; /* only used for integrating the adjoint system */
	}
}
