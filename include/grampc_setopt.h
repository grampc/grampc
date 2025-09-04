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



#ifndef GRAMPC_SETOPT_H_
#define GRAMPC_SETOPT_H_

#include "grampc_init.h"
#include "grampc_mess.h"
#include "grampc_util.h"

void grampc_setopt_real(const typeGRAMPC *grampc, const typeChar *optName, ctypeRNum optValue);
void grampc_setopt_int(const typeGRAMPC *grampc, const typeChar *optName, ctypeInt optValue);
void grampc_setopt_string(const typeGRAMPC *grampc, const typeChar *optName, const typeChar *optValue);
void grampc_setopt_real_vector(const typeGRAMPC *grampc, const typeChar *optName, ctypeRNum *optValue);
void grampc_setopt_int_vector(const typeGRAMPC *grampc, const typeChar *optName, ctypeInt *optValue);
void grampc_printopt(const typeGRAMPC *grampc);

#endif /* GRAMPC_SETOPT_H_ */
