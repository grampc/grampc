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


#include "trapezodial.h"


void trapezodial(typeRNum *s, ctypeRNum *t, ctypeRNum *x, ctypeRNum *u,
	ctypeRNum *p_, const typeGRAMPC *grampc)
{
	typeInt i;
	typeRNum h;

	ctypeRNum *mult = grampc->rws->mult;
	ctypeRNum *pen = grampc->rws->pen;
	ctypeRNum *cfct = grampc->rws->cfct;

	typeRNum *s1 = grampc->rws->rwsGeneral; /* size: 2*/

	s[0] = 0;
	s[1] = 0;

	/* Integration */
	for (i = 0; i < grampc->opt->Nhor; i++) {

		WintCost(s1, t[i], x + i * grampc->param->Nx, u + i * grampc->param->Nu, p_,
			mult + i * grampc->param->Nc, pen + i * grampc->param->Nc, cfct + i * grampc->param->Nc, grampc);

		if (i == 0) {
			h = (t[i + 1] - t[i]) / 2;
		}
		else if (i <= grampc->opt->Nhor - 2) {
			h = (t[i + 1] - t[i - 1]) / 2;
		}
		else {
			h = (t[i] - t[i - 1]) / 2;
		}

		s[0] = s[0] + h * s1[0];
		s[1] = s[1] + h * s1[1];
	}
}
