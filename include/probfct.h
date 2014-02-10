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
 * File: probfct.h
 * Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Date: February 2014
 * Version: v1.0
 *
 * HEADER FILE
 * Header file of problem function for GRAMPC.
 *
 */


#ifndef PROBFCT_C_
#define PROBFCT_C_

#include "grampc_init.h"

/* Function definitions */
void sysdim(typeInt *Nx, typeInt *Nu);
void sysfct(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *pSys);
void sysjacxadj(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *adj,
    typeRNum *u, typeRNum *pSys);
void sysjacuadj(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *adj,
    typeRNum *u, typeRNum *pSys);
void sysjacx(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *pSys);
void sysjacu(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *pSys);
void icostfct(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *xdes,
    typeRNum *udes, typeRNum *pCost);
void icostjacx(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *xdes,
    typeRNum *udes, typeRNum *pCost);
void icostjacu(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *u, typeRNum *xdes,
    typeRNum *udes, typeRNum *pCost);
void fcostfct(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *xdes, typeRNum *pCost);
void fcostjacx(typeRNum *out, typeRNum t, typeRNum *x, typeRNum *xdes, typeRNum *pCost);

#endif /* PROBFCT_C_ */
