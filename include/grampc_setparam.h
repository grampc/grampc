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



#ifndef GRAMPC_SETPARAM_H_
#define GRAMPC_SETPARAM_H_

#include "grampc_init.h"
#include "grampc_mess.h"
#include "grampc_util.h"

void grampc_setparam_real(const typeGRAMPC *grampc, const typeChar *paramName, ctypeRNum paramValue);
void grampc_setparam_real_vector(const typeGRAMPC *grampc, const typeChar *paramName, ctypeRNum *paramValue);
void grampc_printparam(const typeGRAMPC *grampc);

#endif /* GRAMPC_SETPARAM_H_ */
