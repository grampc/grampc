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


#ifndef RODAS_H_
#define RODAS_H_

#include "grampc_init.h"
#include "grampc_util.h"
#include "grampc_mess.h"
#include "probfct.h"



void intsysRodas(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum *tvec, ctypeRNum *xvec, ctypeRNum *uvec, ctypeRNum *pvec, const typeGRAMPC *grampc, const typeffctPtr pfct);
void ffctRodas(int *N, typeRNum *t, typeRNum *x, typeRNum *rhs, ctypeRNum *tvec, ctypeRNum *xvec, ctypeRNum *uvec, ctypeRNum *pvec, ctypeRNum *dcdxvec, const typeGRAMPC *grampc, const typeffctPtr pfct);
void dfdxRodas(int *N, typeRNum *t, typeRNum *x, typeRNum *dfdx_val, int *ldfdy, ctypeRNum *tvec, ctypeRNum *xvec, ctypeRNum *uvec, ctypeRNum *pvec, ctypeRNum *dcdxvec, const typeGRAMPC *grampc, const typeffctPtr pfct);
void dfdtRodas(int *N, typeRNum *t, typeRNum *x, typeRNum *dfdt_val, int *ldfdy, ctypeRNum *tvec, ctypeRNum *xvec, ctypeRNum *uvec, ctypeRNum *pvec, ctypeRNum *dcdxvec, const typeGRAMPC *grampc, const typeffctPtr pfct);
void MfctRodas(int *N, typeRNum *out, int *LMAS, const typeGRAMPC *grampc, const typeffctPtr pfct);
void solout(int *nr, typeRNum *xold, typeRNum *x, typeRNum *h, typeRNum *y, typeRNum *cont, int *lrc, int *n, const typeGRAMPC *grampc, const typeffctPtr pfct, int *irtrn);

#endif /* RODAS_H_ */
