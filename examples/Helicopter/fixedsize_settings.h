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
#ifndef FIXEDSIZE_SETTINGS_H
#define FIXEDSIZE_SETTINGS_H

/* Define this macro if GRAMPC should be compiled without dynamic memory allocation */
#ifdef FIXEDSIZE

/* Required header for definitions of integrators and line search types */
#include "grampc_macro.h"

/* Define this macro as the number of states Nx */
#define NX 6
/* Define this macro as the number of controls Nu */
#define NU 2
/* Define this macro as the number of parameters Np */
#define NP 0
/* Define this macro as the number of equality constraints Ng */
#define NG 0
/* Define this macro as the number of inequality constraints Nh */
#define NH 2
/* Define this macro as the number of terminal equality constraints NgT */
#define NGT 0
/* Define this macro as the number of terminal inequality constraints NhT */
#define NHT 0

/* Define this macro as the number of gradient iterations MaxGradIter */
#define MAXGRADITER 3
/* Define this macro as the number of multiplier iterations MaxMultIter */
#define MAXMULTITER 1
/* Define this macro as the number of discretization points Nhor */
#define NHOR 30
/* Define this macro as one of the integrators INT_EULER, INT_HEUN, INT_MODEULER, INT_RUKU45, INT_RODAS */
#define INTEGRATOR INT_HEUN
/* Define this macro as one of the cost integrators INT_TRAPZ, INT_SIMPSON */
#define INTEGRATORCOST INT_TRAPZ
/* Define this macro as one the line search strategies INT_ADAPTIVELS, INT_EXPLS1, INT_EXPLS2 */
#define LINESEARCHTYPE INT_EXPLS2

/* The following macros are only relevant if the integrator RODAS is used */
/* Define this macro as 1 if right hand side depends on time t */
#define RODAS_IFCN 0
/* Define this macro as 1 if partial derivatives dfdt and dHdxdt are supplied */
#define RODAS_IDFX 0
/* Define this macro as 1 if Jacobians dfdx and dfdxtrans are supplied */
#define RODAS_IJAC 0
/* Define this macro as 1 if mass matrix Mfct is supplied */
#define RODAS_IMAS 0
/* Define this macro as the number of lower diagonals of Jacobian */
#define RODAS_MLJAC NX
/* Define this macro as the number of upper diagonals of Jacobian */
#define RODAS_MUJAC NX
/* Define this macro as the number of lower diagonals of mass matrix */
#define RODAS_MLMAS NX
/* Define this macro as the number of upper diagonals of mass matrix */
#define RODAS_MUMAS NX

#endif /* FIXEDSIZE */

#endif /* FIXEDSIZE_SETTINGS_H */
