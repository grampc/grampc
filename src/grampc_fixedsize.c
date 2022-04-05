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

/* suppress empty translation unit warning */
typedef int make_iso_compilers_happy;

#ifdef FIXEDSIZE

#include "grampc_alloc.h"
#include "grampc_mess.h"
#include "probfct.h"


void grampc_alloc_structs(typeGRAMPC **grampc, typeUSERPARAM *userparam)
{
    /* this replaces grampc_alloc_structs from grampc_alloc.c in fixed-size mode */

    /* set userparam for compatibility with grampc_alloc.c */
    (*grampc)->userparam = userparam;

    /* call ocp_dim once during initialization for compatibility with grampc_alloc.c */
    ocp_dim(&(*grampc)->param->Nx, &(*grampc)->param->Nu, &(*grampc)->param->Np, &(*grampc)->param->Ng, &(*grampc)->param->Nh, &(*grampc)->param->NgT, &(*grampc)->param->NhT, (*grampc)->userparam);

    /* check problem dimensions */
    if((*grampc)->param->Nx != NX)
    {
        grampc_error(INVALID_NX);
    }
    if((*grampc)->param->Nu != NU)
    {
        grampc_error(INVALID_NU);
    }
    if((*grampc)->param->Np != NP)
    {
        grampc_error(INVALID_NP);
    }
    if((*grampc)->param->Ng != NG)
    {
        grampc_error(INVALID_Ng);
    }
    if((*grampc)->param->Nh != NH)
    {
        grampc_error(INVALID_Nh);
    }
    if((*grampc)->param->NgT != NGT)
    {
        grampc_error(INVALID_NgT);
    }
    if((*grampc)->param->NhT != NHT)
    {
        grampc_error(INVALID_NhT);
    }
    (*grampc)->param->Nc = (*grampc)->param->Ng + (*grampc)->param->Nh + (*grampc)->param->NgT + (*grampc)->param->NhT;
}


void grampc_alloc_fields(typeGRAMPC **grampc, typeUSERPARAM *userparam)
{
    /* this replaces grampc_alloc_fields from grampc_alloc.c in fixed-size mode */
}


void grampc_free(typeGRAMPC **grampc)
{
    /* this replaces grampc_free from grampc_alloc.c in fixed-size mode */
}

#endif /* FIXEDSIZE */
