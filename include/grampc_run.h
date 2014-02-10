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
 * File: grampc_run.h
 * Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Date: February 2014
 * Version: v1.0
 *
 * HEADER FILE
 * MPC running file for GRAMPC.
 *
 */


#ifndef GRAMPC_RUN_H_
#define GRAMPC_RUN_H_

#include "grampc_init.h"
#include "euler1.h"
#include "eulermod2.h"
#include "heun2.h"
#include "ruku45.h"

/* Definitions of functions */
void grampc_run(typeGRAMPC *grampc);
void dHdufct(typeGRAMPC *grampc);
void inputproj(typeRNum *u, typeGRAMPC *grampc);
void lsearch_fit2(typeRNum *kfit, typeRNum *Jfit, typeRNum *k, typeRNum *J);
void interplin(typeRNum *varint, typeRNum *tvec, typeRNum *varvec, typeRNum tint,
    typeInt Nvar, typeInt Nvec, typeInt searchdir);
void MatAdd(typeRNum *C, typeRNum *A, typeRNum *B, typeInt n1, typeInt n2);
void MatMult(typeRNum *C, typeRNum *A, typeRNum *B, typeInt n1, typeInt n2, typeInt n3);
void minfct(typeRNum *amin, typeInt *amini, typeRNum *a, typeInt Na);
void maxfct(typeRNum *amax, typeInt *amaxi, typeRNum *a, typeInt Na);
void Wadjsys(typeRNum *s, typeRNum *adj, typeRNum *t, typeRNum *x, typeRNum *u, typeGRAMPC *grampc);
void Wsys(typeRNum *s, typeRNum *x, typeRNum *t, typeRNum *dummy, typeRNum *u, typeGRAMPC *grampcs);
void intCostTrapezodial(typeRNum *s, typeRNum *t, typeRNum *x, typeRNum *u, typeGRAMPC *grampc);
void intCostSimpson(typeRNum *s, typeRNum *t, typeRNum *x, typeRNum *u, typeGRAMPC *grampc);

#endif /* GRAMPC_RUN_H_ */
