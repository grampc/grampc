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


#ifndef TRAPEZODIAL_H
#define TRAPEZODIAL_H

#include "grampc_init.h"
#include "grampc_run.h"


void trapezodial(typeRNum *s, ctypeRNum *t, ctypeRNum *x, ctypeRNum *u,
	ctypeRNum *p_, const typeGRAMPC *grampc);


#endif
