/*
 *
 * This file is part of GRAMPC.
 *
 * GRAMPC - a gradient-based MPC software for real-time applications
 *
 * Copyright (C) 2014 by Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Developed at the Institute of Measurement, Control, and
 * Microtechnology, University of Ulm. All rights reserved.
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
 * License along with GRAMPC. If not, see <http://www.gnu.org/licenses/>.
 *
 */


/*
 *
 * File: grampc_setparam.h
 * Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Date: February 2014
 * Version: v1.0
 *
 * HEADER FILE
 * Parameter setting file for GRAMPC.
 *
 */


#ifndef GRAMPC_SETPARAM_H_
#define GRAMPC_SETPARAM_H_

#include "grampc_init.h"

void grampc_setparam_real(typeGRAMPC *grampc, typeChar paramName[], typeRNum paramValue);
void grampc_setparam_int(typeGRAMPC *grampc, typeChar paramName[], typeInt paramValue);
void grampc_setparam_vector(typeGRAMPC *grampc, typeChar paramName[], typeRNum *paramValue);
void grampc_printparam(typeGRAMPC *grampc);

#endif /* GRAMPC_SETPARAM_H_ */
