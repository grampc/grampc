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


#ifndef RUKU45_H_
#define RUKU45_H_

#include "grampc_init.h"
#include "grampc_util.h"
#include "grampc_mess.h"

void intsysRuKu45(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum *t, ctypeRNum *x,
	ctypeRNum *u, ctypeRNum *p_, const typeGRAMPC *grampc, const typeffctPtr pfct);

#endif /* RUKU45_H_ */
