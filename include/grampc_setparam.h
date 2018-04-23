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



#ifndef GRAMPC_SETPARAM_H_
#define GRAMPC_SETPARAM_H_

#include "grampc_init.h"
#include "grampc_mess.h"
#include "grampc_util.h"

void grampc_setparam_real(const typeGRAMPC *grampc, const typeChar *paramName, ctypeRNum paramValue);
void grampc_setparam_real_vector(const typeGRAMPC *grampc, const typeChar *paramName, ctypeRNum *paramValue);
void grampc_printparam(const typeGRAMPC *grampc);

#endif /* GRAMPC_SETPARAM_H_ */
