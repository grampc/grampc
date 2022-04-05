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


#include "simpson.h"


void simpson(typeRNum *s, ctypeRNum *t, ctypeRNum *x, ctypeRNum *u,
	ctypeRNum *p_, const typeGRAMPC *grampc)
{
	typeInt i, j;
	typeRNum h;

	ctypeRNum *mult = grampc->rws->mult;
	ctypeRNum *pen = grampc->rws->pen;
	ctypeRNum *cfct = grampc->rws->cfct;

	typeRNum *s1 = grampc->rws->rwsGeneral; /* size: 5+Nx+Nu+3*Nc */
	typeRNum *s2 = s1 + 2;
	typeRNum *ts = s2 + 2;
	typeRNum *xs = ts + 1;
	typeRNum *us = xs + grampc->param->Nx;
	typeRNum *lags = us + grampc->param->Nu;
	typeRNum *pens = lags + grampc->param->Nc;
	typeRNum *cfcts = pens + grampc->param->Nc;

	s[0] = 0;
	s[1] = 0;

	for (i = 0; i < grampc->opt->Nhor; i++) {

		/* Integration */
		s1[0] = 0;
		s1[1] = 0;
		WintCost(s1, t[i], x + i * grampc->param->Nx, u + i * grampc->param->Nu, p_,
			mult + i * grampc->param->Nc, pen + i * grampc->param->Nc, cfct + i * grampc->param->Nc, grampc);

		if (i == 0) {
			h = (t[i + 1] - t[i]) / 6;
		}
		else if (i <= grampc->opt->Nhor - 2) {
			h = (t[i + 1] - t[i - 1]) / 6;
		}
		else {
			h = (t[i] - t[i - 1]) / 6;
		}
		s[0] = s[0] + h * s1[0];
		s[1] = s[1] + h * s1[1];

		/* Intermediate points */
		if (i < grampc->opt->Nhor - 1) {

			/* Interpolation */
			ts[0] = (t[i] + t[i + 1]) / 2;
			for (j = 0; j < grampc->param->Nx; j++) {
				xs[j] = (x[i * grampc->param->Nx + j] + x[(i + 1) * grampc->param->Nx + j]) / 2;
			}
			for (j = 0; j < grampc->param->Nu; j++) {
				us[j] = (u[i * grampc->param->Nu + j] + u[(i + 1) * grampc->param->Nu + j]) / 2;
			}
			for (j = 0; j < grampc->param->Nc; j++) {
				lags[j] = (mult[i * grampc->param->Nc + j] + mult[(i + 1) * grampc->param->Nc + j]) / 2;
				pens[j] = (pen[i * grampc->param->Nc + j] + pen[(i + 1) * grampc->param->Nc + j]) / 2;
				cfcts[j] = (cfct[i * grampc->param->Nc + j] + cfct[(i + 1) * grampc->param->Nc + j]) / 2;
			}

			/* Integration */
			s2[0] = 0;
			s2[1] = 0;
			WintCost(s2, ts[0], xs, us, p_, lags, pens, cfcts, grampc);

			h = 4 * (t[i + 1] - t[i]) / 6;
			s[0] = s[0] + h * s2[0];
			s[1] = s[1] + h * s2[1];
		}
	}
}
