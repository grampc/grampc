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


#ifndef GRAMPC_ERK_H_
#define GRAMPC_ERK_H_

#include "grampc_init.h"

/* First-order explicit Runge-Kutta integration (Euler's method) */
void intsysERK1(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum * t, ctypeRNum *x,
	ctypeRNum *u, ctypeRNum *p, ctypeRNum *dcdx, const typeGRAMPC *grampc, const typeSysPtr pfct);

/* Second-order explicit Runge-Kutta integration (Heun's method) */
void intsysERK2(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum * t, ctypeRNum *x,
	ctypeRNum *u, ctypeRNum *p, ctypeRNum *dcdx, const typeGRAMPC *grampc, const typeSysPtr pfct);

/* Third-order explicit Runge-Kutta integration (Kutta's third-order method) */
void intsysERK3(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum * t, ctypeRNum *x,
	ctypeRNum *u, ctypeRNum *p, ctypeRNum *dcdx, const typeGRAMPC *grampc, const typeSysPtr pfct);

/* Fourth-order explicit Runge-Kutta integration (Kutta's classic fourth-order method )*/
void intsysERK4(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum * t, ctypeRNum *x,
	ctypeRNum *u, ctypeRNum *p, ctypeRNum *dcdx, const typeGRAMPC *grampc, const typeSysPtr pfct);
	
#endif /* GRAMPC_ERK_H_ */
