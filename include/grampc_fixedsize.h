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
#ifndef GRAMPC_FIXEDSIZE_H
#define GRAMPC_FIXEDSIZE_H

#ifdef FIXEDSIZE

/* The user is responsible that this header file is found by the compiler */
#include "fixedsize_settings.h"

/* Validate parameters and options */
#ifndef NX
#error "The number of states NX is not defined."
#endif

#ifndef NU
#error "The number of controls NU is not defined."
#endif

#ifndef NP
#error "The number of parameters NP is not defined."
#endif

#ifndef NG
#error "The number of equality constraints NG is not defined."
#endif

#ifndef NH
#error "The number of inequality constraints NH is not defined."
#endif

#ifndef NGT
#error "The number of terminal equality constraints NGT is not defined."
#endif

#ifndef NHT
#error "The number of terminal inequality constraints NHT is not defined."
#endif

#ifndef MAXGRADITER
#error "The number of gradient iterations MAXGRADITER is not defined."
#endif

#ifndef MAXMULTITER
#error "The number of multiplier iterations MAXMULTITER is not defined."
#endif

#ifndef NHOR
#error "The number of discretization points NHOR is not defined."
#endif

#ifndef INTEGRATOR
#error "The integrator INTEGRATOR is not defined."
#endif

#ifndef INTEGRATORCOST
#error "The cost integrator INTEGRATORCOST is not defined."
#endif

#ifndef LINESEARCHTYPE
#error "The line search LINESEARCHTYPE is not defined."
#endif

/* Total number of constraints */
#define NC (NG + NH + NGT + NHT)

/* Compute required size of field lsAdapt */
#if (LINESEARCHTYPE == INT_ADAPTIVELS)
#define SIZE_LSADAPTIVE (2 * (NALS + 1) * (1 + MAXGRADITER))
#else
#define SIZE_LSADAPTIVE 0
#endif

/* Compute required size of field lsExplicit */
#if (LINESEARCHTYPE == INT_EXPLS1 || LINESEARCHTYPE == INT_EXPLS2)
#define SIZE_LSEXPLICIT (NELS)
#else
#define SIZE_LSEXPLICIT 0
#endif

/* Compute required size of field rwsScale */
#define SIZE_RWSSCALE (2 * (NX + NU + NP))

/* Compute required size of field rwsGeneral */
#define SIZE_WADJSYS (NX)
#define SIZE_EULER (NX)
#define SIZE_MODEULER (5 * NX + NU + NC)
#define SIZE_HEUN (3 * NX)
#define SIZE_RUKU45 (18 * NX + NU)
#define SIZE_RODAS (2 * NX + NU)
#define SIZE_SIMPSON (NX + NU + 3 * NC + 5)
#define SIZE_TRAPEZOIDAL (2)
#define SIZE_GRADU (2 * NU)
#define SIZE_GRADP (3 * NP)
#define SIZE_GRADT (NX)

/* Size for integration depends on integrator type */
#if (INTEGRATOR == INT_EULER)
#define SIZE_INTEGRATOR (SIZE_WADJSYS + SIZE_EULER)
#elif (INTEGRATOR == INT_MODEULER)
#define SIZE_INTEGRATOR (SIZE_WADJSYS + SIZE_MODEULER)
#elif (INTEGRATOR == INT_HEUN)
#define SIZE_INTEGRATOR (SIZE_WADJSYS + SIZE_HEUN)
#elif (INTEGRATOR == INT_RUKU45)
#define SIZE_INTEGRATOR (SIZE_WADJSYS + SIZE_RUKU45)
#elif (INTEGRATOR == INT_RODAS)
#define SIZE_INTEGRATOR (SIZE_WADJSYS + SIZE_RODAS)
#else
#error "Invalid value for INTEGRATOR, define as one of INT_EULER, INT_MODEULER, INT_HEUN, INT_RUKU45, INT_RODAS."
#endif

/* Size for cost integration depends on integrator type */
#if (INTEGRATORCOST == INT_TRAPZ)
#define SIZE_INTEGRATORCOST SIZE_TRAPEZOIDAL
#elif (INTEGRATORCOST == INT_SIMPSON)
#define SIZE_INTEGRATORCOST SIZE_SIMPSON
#else
#error "Invalid value for INTEGRATORCOST, define as one of INT_TRAPZ, INT_SIMPSON."
#endif

/* Size for constraints depends on number of constraints */
#if (NC > 0)
#define SIZE_CONSTRAINTS (NC + 2 * (NX + NU + NP))
#else
#define SIZE_CONSTRAINTS 0
#endif

#define SIZE_RWSGENERAL (MAX(MAX(MAX(MAX(MAX(SIZE_GRADP, SIZE_GRADU), SIZE_GRADT), SIZE_INTEGRATORCOST), SIZE_INTEGRATOR), SIZE_CONSTRAINTS))

/* Compute required size of rodas-specific fields */
#if (INTEGRATOR == INT_RODAS)

#if (RODAS_MLJAC < NX)
#define RODAS_LJAC (RODAS_MLJAC + RODAS_MUJAC + 1)
#define RODAS_LE1 (2 * RODAS_MLJAC + RODAS_MUJAC + 1)
#else
#define RODAS_LJAC (NX)
#define RODAS_LE1 (NX)
#endif

#if (RODAS_IMAS == 0)
#define RODAS_LMAS 0
#elif (RODAS_MLMAS == NX)
#define RODAS_LMAS (NX)
#else
#define RODAS_LMAS (RODAS_MLMAS + RODAS_MUMAS + 1)
#endif

#define SIZE_RPARRODAS (NX * NHOR)
#define SIZE_IPARRODAS 20
#define SIZE_RWORKRODAS (NX * (RODAS_LJAC + RODAS_LMAS + RODAS_LE1 + 14) + 20)
#define SIZE_IWORKRODAS (20 + NX)

#else

#define SIZE_RPARRODAS 0
#define SIZE_IPARRODAS 0
#define SIZE_RWORKRODAS 0
#define SIZE_IWORKRODAS 0

#endif /* INTEGRATOR == INT_RODAS */

/* Macro for creating grampc structure and pointer 'name' without dynamic memory allocation */
#define TYPE_GRAMPC_POINTER(name) \
    typeGRAMPCparam name##_param = {}; \
    typeGRAMPCopt name##_opt = {}; \
    typeGRAMPCsol name##_sol = {}; \
    typeGRAMPCrws name##_rws = {}; \
    typeGRAMPC name##_struct = {}; \
    name##_struct.param = &name##_param; \
    name##_struct.opt = &name##_opt; \
    name##_struct.sol = &name##_sol; \
    name##_struct.rws = &name##_rws; \
    typeGRAMPC* name = &name##_struct;

#endif /* FIXEDSIZE */

#endif /* GRAMPC_FIXEDSIZE_H */
