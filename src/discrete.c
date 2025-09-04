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
#include "discrete.h"

void intsysDiscrete(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum *t, ctypeRNum *x,
    ctypeRNum *u, ctypeRNum *p, ctypeRNum *dcdx, const typeGRAMPC *grampc, const typeSysPtr  pfct)
{
    typeInt i, j;
    typeRNum *tmp = grampc->rws->rwsGeneral + LWadjsys; /* size Nx */

    if (pInt == FWINT) {
        /* forward propagation of system dynamics */
        /* x_{k+1} = f(t_k, x_k, u_k), k = 0, ..., N-2 */
        for (i = 0; i < Nint - 1; i++) {
            (*pfct)(y + grampc->param->Nx, y, t, NULL, u, p, NULL, grampc);

            t += 1;
            x += grampc->param->Nx;
            u += grampc->param->Nu;
            y += grampc->param->Nx;
            /* note that the last control u_{N-1} is unused in the state propagation */
        }
    }
    else {
        /* backward propagation of adjoint dynamics */
        /* \lambda_k = H_x(t_k, x_k, u_k, \lambda_{k+1}), k = N-2, ..., 0 */
        for (i = 0; i < Nint - 1; i++) {
            t -= 1;
            x -= grampc->param->Nx;
            u -= grampc->param->Nu;
            y -= grampc->param->Nx;
            dcdx -= grampc->param->Nx;
            /* note that the last control u_{N-1} is unused in the adjoint state propagation */
            /* note that the last state x_{N-1} affects \lambda_{N-1} only through the terminal condition */

            (*pfct)(tmp, y + grampc->param->Nx, t, x, u, p, dcdx, grampc);

            /* change sign with respect to Wadjsys */
            for (j = 0; j < grampc->param->Nx; ++j) {
                y[j] = -tmp[j];
            }
        }
    }
}

void intcostDiscrete(typeRNum *s, ctypeRNum *t, ctypeRNum *x, ctypeRNum *u,
                     ctypeRNum *p, const typeGRAMPC *grampc)
{
    typeInt i;
    ctypeRNum *mult = grampc->rws->mult;
    ctypeRNum *pen = grampc->rws->pen;
    ctypeRNum *cfct = grampc->rws->cfct;
    typeRNum *tmp = grampc->rws->rwsGeneral; /* size: 2*/

    /* Integration */
    s[0] = 0;
    s[1] = 0;
    for (i = 0; i < grampc->opt->Nhor; i++) {

        WintCost(tmp, t[0], x, u, p, mult, pen, cfct, grampc);

        s[0] = s[0] + tmp[0];
        s[1] = s[1] + tmp[1];

        t += 1;
        x += grampc->param->Nx;
        u += grampc->param->Nu;
        mult += grampc->param->Nc;
        pen += grampc->param->Nc;
        cfct += grampc->param->Nc;
    }
}
