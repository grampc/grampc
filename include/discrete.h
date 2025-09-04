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
#ifndef DISCRETE_H_
#define DISCRETE_H_

#include "grampc_init.h"
#include "grampc_run.h"

void intsysDiscrete(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum *t, ctypeRNum *x,
    ctypeRNum *u, ctypeRNum *p, ctypeRNum *dcdx, const typeGRAMPC *grampc, const typeSysPtr pfct);

void intcostDiscrete(typeRNum *s, ctypeRNum *t, ctypeRNum *x, ctypeRNum *u,
                     ctypeRNum *p, const typeGRAMPC *grampc);

#endif /* DISCRETE_H_ */
