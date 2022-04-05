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
