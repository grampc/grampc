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
#ifndef GRAMPC_ALLOC_H
#define GRAMPC_ALLOC_H

#include "grampc_init.h"


/* macro for creating typeGRAMPC* pointer that is replaced in fixed-size mode */
#ifndef FIXEDSIZE
#define TYPE_GRAMPC_POINTER(name) typeGRAMPC* name;
#endif

void createNumMatrix(typeRNum **cs, const size_t size);
void createIntMatrix(typeInt **cs, const size_t size);
void resizeNumMatrix(typeRNum **cs, const size_t size);
void resizeIntMatrix(typeInt **cs, const size_t size);

void resize_rwsLinesearch(const typeGRAMPC *grampc);
void resize_rwsGeneral(const typeGRAMPC *grampc);
void resize_rwsRodas(const typeGRAMPC *grampc);
void setLWorkRodas(const typeGRAMPC *grampc);

typeInt CastDvec2Intvec(typeInt** Intvec, const double* Numvec, const size_t size);
typeInt CastDvec2Numvec(typeRNum** Realvec, const double* Numvec, const size_t size);

void grampc_alloc_structs(typeGRAMPC **grampc, typeUSERPARAM *userparam);
void grampc_alloc_fields(typeGRAMPC **grampc, typeUSERPARAM *userparam);
void grampc_free_param(typeGRAMPCparam **param);
void grampc_free_opt(typeGRAMPCopt **opt);
void grampc_free_sol(typeGRAMPCsol **sol);
void grampc_free_rws(typeGRAMPCrws **rws);
void grampc_free(typeGRAMPC **grampc);

#endif /* GRAMPC_ALLOC_H */
