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
#ifndef FIXEDSIZE

#include "grampc_alloc.h"
#include "grampc_mess.h"
#include "probfct.h"


void createNumMatrix(typeRNum **cs, const size_t size) {
    if (size == 0) {
        *cs = NULL;
    }
    else {
        *cs = (typeRNum *)calloc(size, sizeof(typeRNum));
        if (*cs == NULL) {
            grampc_error(PARAM_ALLOC_FAILED);
        }
    }
}


void createIntMatrix(typeInt **cs, const size_t size) {
    if (size == 0) {
        *cs = NULL;
    }
    else {
        *cs = (typeInt *)calloc(size, sizeof(typeInt));
        if (*cs == NULL) {
            grampc_error(PARAM_ALLOC_FAILED);
        }
    }
}


void resizeNumMatrix(typeRNum **cs, const size_t size) {
    free(*cs);
    if (size == 0) {
        *cs = NULL;
    }
    else {
        *cs = (typeRNum *)calloc(size, sizeof(typeRNum));
        if (*cs == NULL) {
            grampc_error(RWS_ELEMENT_ALLOC_FAILED);
        }
    }
}


void resizeIntMatrix(typeInt **cs, const size_t size) {
    free(*cs);
    if (size == 0) {
        *cs = NULL;
    }
    else {
        *cs = (typeInt *)calloc(size, sizeof(typeInt));
        if (*cs == NULL) {
            grampc_error(SOL_ALLOC_FAILED);
        }
    }
}


void resize_rwsLinesearch(const typeGRAMPC *grampc)
{
    if (grampc->opt->LineSearchType == INT_ADAPTIVELS) {
        resizeNumMatrix(&grampc->rws->lsAdapt, 2 * (NALS + 1)*(1 + grampc->opt->MaxGradIter));
        resizeNumMatrix(&grampc->rws->lsExplicit, 0);
    }
    else {
        resizeNumMatrix(&grampc->rws->lsAdapt, 0);
        resizeNumMatrix(&grampc->rws->lsExplicit, NELS);
    }
}


void resize_rwsGeneral(const typeGRAMPC *grampc) {
    typeInt LInt = LWadjsys;
    typeInt LIntCost = 0;
    typeInt LConst = 0;

    switch (grampc->opt->Integrator) {
    case INT_SYSDISC:  LInt += Lerk1; break;
    case INT_ERK1:     LInt += Lerk1; break;
    case INT_ERK2:     LInt += Lerk2; break;
    case INT_ERK3:     LInt += Lerk3; break;
    case INT_ERK4:     LInt += Lerk4; break;
    case INT_RUKU45:   LInt += Lruku45; break;
    case INT_RODAS:	   LInt += Lrodas; break;
    }
    switch (grampc->opt->IntegratorCost) {
    case INT_COSTDISC: LIntCost = LIntCostTrapezoidal; break;
    case INT_TRAPZ:    LIntCost = LIntCostTrapezoidal; break;
    case INT_SIMPSON:  LIntCost = LIntCostSimpson; break;
    }
    if (grampc->param->Nc > 0) {
        LConst = LevaluateConstraints;
    }

    grampc->rws->lrwsGeneral = MAX(MAX(MAX(MAX(MAX(Lgradp, Lgradu), LgradT), LIntCost), LInt), LConst);
    resizeNumMatrix(&grampc->rws->rwsGeneral, grampc->rws->lrwsGeneral);
}


void resize_rwsRodas(const typeGRAMPC *grampc) {
    if (grampc->opt->Integrator == INT_RODAS) {
        resizeNumMatrix(&grampc->rws->rparRodas, grampc->param->Nx*grampc->opt->Nhor);
        resizeIntMatrix(&grampc->rws->iparRodas, 20);
        resizeNumMatrix(&grampc->rws->workRodas, grampc->rws->lworkRodas);
        resizeIntMatrix(&grampc->rws->iworkRodas, grampc->rws->liworkRodas);
    }
    else {
        resizeNumMatrix(&grampc->rws->rparRodas, 0);
        resizeIntMatrix(&grampc->rws->iparRodas, 0);
        resizeNumMatrix(&grampc->rws->workRodas, 0);
        resizeIntMatrix(&grampc->rws->iworkRodas, 0);
    }
}


void setLWorkRodas(const typeGRAMPC *grampc)
{
    typeInt LJAC, LMAS, LE1;							/* length of rodas work space*/
/* ctypeInt IFCN = grampc->opt->FlagsRodas[0];*/		/* 0 --> right hand side independent of time t */
/* ctypeInt IDFX = grampc->opt->FlagsRodas[1];*/		/* 0 --> DF/DX is numerically computed */
/* typeInt IJAC = grampc->opt->FlagsRodas[2];*/			/* 1(0) -> analytical (numerical) jacobian (partial derivatives of right hand side w.r.t. state) */

    ctypeInt MLJAC = grampc->opt->FlagsRodas[4];		/* no. of lower diagonals of jacobian */
    ctypeInt MUJAC = grampc->opt->FlagsRodas[5];		/* no. of upper diagonals of jacobian */

    typeInt IMAS = grampc->opt->FlagsRodas[3];			/* 1 --> mass matrix is supplied */
    ctypeInt MLMAS = grampc->opt->FlagsRodas[6];		/* no. of lower diagonals of mass matrix */
    ctypeInt MUMAS = grampc->opt->FlagsRodas[7];		/* no. of upper diagonals of mass matrix */

    if (MLJAC < grampc->param->Nx)
    {
        LJAC = MLJAC + MUJAC + 1;
        LE1 = 2 * MLJAC + MUJAC + 1;
    }
    else
    {
        LJAC = grampc->param->Nx;
        LE1 = grampc->param->Nx;
    }

    if (IMAS == 0)
    {
        LMAS = 0;
    }
    else
    {
        if (MLMAS == grampc->param->Nx)
        {
            LMAS = grampc->param->Nx;
        }
        else
        {
            LMAS = MLMAS + MUMAS + 1;
        }
    }
    grampc->rws->lworkRodas = grampc->param->Nx*(LJAC + LMAS + LE1 + 14) + 20;
    resizeNumMatrix(&grampc->rws->workRodas, grampc->rws->lworkRodas);
}


typeInt CastDvec2Intvec(typeInt** Intvec, const double* Numvec, const size_t size) {
    unsigned typeInt i;
    *Intvec = (typeInt*)malloc(size * sizeof(typeInt));
    if (*Intvec != NULL) {
        for (i = 0; i < size; i++) {
            (*Intvec)[i] = (typeInt)Numvec[i];
        }
        return 1;
    }
    return -1;
}


typeInt CastDvec2Numvec(typeRNum** Realvec, const double* Numvec, const size_t size) {
    unsigned typeInt i;
    *Realvec = (typeRNum*)malloc(size * sizeof(typeRNum));
    if (*Realvec != NULL) {
        for (i = 0; i < size; i++) {
            (*Realvec)[i] = (typeRNum)Numvec[i];
        }
        return 1;
    }
    return -1;
}


void grampc_alloc_structs(typeGRAMPC **grampc, typeUSERPARAM *userparam)
{
    /* STRUCTURE MEMORY ALLOCATION **********************************************/
    *grampc = (typeGRAMPC *)calloc(1, sizeof(**grampc));
    if (*grampc == NULL) {
        grampc_error(GRAMPC_ALLOC_FAILED);
    }
    (*grampc)->param = (typeGRAMPCparam *)calloc(1, sizeof(*(*grampc)->param));
    if ((*grampc)->param == NULL) {
        grampc_error(PARAM_ALLOC_FAILED);
    }
    (*grampc)->sol = (typeGRAMPCsol *)calloc(1, sizeof(*(*grampc)->sol));
    if ((*grampc)->sol == NULL) {
        grampc_error(SOL_ALLOC_FAILED);
    }
    (*grampc)->rws = (typeGRAMPCrws *)calloc(1, sizeof(*(*grampc)->rws));
    if ((*grampc)->rws == NULL) {
        grampc_error(RWS_ALLOC_FAILED);
    }
    (*grampc)->opt = (typeGRAMPCopt *)calloc(1, sizeof(*(*grampc)->opt));
    if ((*grampc)->opt == NULL) {
        grampc_error(OPT_ALLOC_FAILED);
    }

    /* USERPARAM ****************************************************************/
    (*grampc)->userparam = userparam;

    /* PARAM STRUCTURE **********************************************************/
    ocp_dim(&(*grampc)->param->Nx, &(*grampc)->param->Nu, &(*grampc)->param->Np, &(*grampc)->param->Ng, &(*grampc)->param->Nh, &(*grampc)->param->NgT, &(*grampc)->param->NhT, (*grampc)->userparam);
    /* check number of states and control */
    if ((*grampc)->param->Nx <= 0) {
        grampc_error(INVALID_NX);
    }
    if ((*grampc)->param->Nu < 0) {
        grampc_error(INVALID_NU);
    }
    if ((*grampc)->param->Np < 0) {
        grampc_error(INVALID_NP);
    }
    if ((*grampc)->param->Ng < 0) {
        grampc_error(INVALID_Ng);
    }
    if ((*grampc)->param->Nh < 0) {
        grampc_error(INVALID_Nh);
    }
    if ((*grampc)->param->NgT < 0) {
        grampc_error(INVALID_NgT);
    }
    if ((*grampc)->param->NhT < 0) {
        grampc_error(INVALID_NhT);
    }
    (*grampc)->param->Nc = (*grampc)->param->Ng + (*grampc)->param->Nh + (*grampc)->param->NgT + (*grampc)->param->NhT;
}


void grampc_alloc_fields(typeGRAMPC **grampc, typeUSERPARAM *userparam)
{
    createNumMatrix(&(*grampc)->param->x0, (*grampc)->param->Nx);
    createNumMatrix(&(*grampc)->param->xdes, (*grampc)->param->Nx);
    createNumMatrix(&(*grampc)->param->u0, (*grampc)->param->Nu);
    createNumMatrix(&(*grampc)->param->udes, (*grampc)->param->Nu);
    createNumMatrix(&(*grampc)->param->umax, (*grampc)->param->Nu);
    createNumMatrix(&(*grampc)->param->umin, (*grampc)->param->Nu);
    createNumMatrix(&(*grampc)->param->p0, (*grampc)->param->Np);
    createNumMatrix(&(*grampc)->param->pmax, (*grampc)->param->Np);
    createNumMatrix(&(*grampc)->param->pmin, (*grampc)->param->Np);

    /* OPTIONS STRUCTURE *******************************************************/
    createNumMatrix(&(*grampc)->opt->xScale, (*grampc)->param->Nx);
    createNumMatrix(&(*grampc)->opt->xOffset, (*grampc)->param->Nx);
    createNumMatrix(&(*grampc)->opt->uScale, (*grampc)->param->Nu);
    createNumMatrix(&(*grampc)->opt->uOffset, (*grampc)->param->Nu);
    createNumMatrix(&(*grampc)->opt->pScale, (*grampc)->param->Np);
    createNumMatrix(&(*grampc)->opt->pOffset, (*grampc)->param->Np);
    createNumMatrix(&(*grampc)->opt->cScale, (*grampc)->param->Nc);
    createNumMatrix(&(*grampc)->opt->ConstraintsAbsTol, (*grampc)->param->Nc);

    /* SOLUTION STRUCTURE ******************************************************/
    createNumMatrix(&(*grampc)->sol->xnext, (*grampc)->param->Nx);
    createNumMatrix(&(*grampc)->sol->unext, (*grampc)->param->Nu);
    createNumMatrix(&(*grampc)->sol->pnext, (*grampc)->param->Np);
    createIntMatrix(&(*grampc)->sol->iter, (*grampc)->opt->MaxMultIter);

    /* RWS STRUCTURE **********************************************************/
    createNumMatrix(&(*grampc)->rws->t, (*grampc)->opt->Nhor);
    createNumMatrix(&(*grampc)->rws->tls, (*grampc)->opt->Nhor);
    createNumMatrix(&(*grampc)->rws->x, (*grampc)->param->Nx*(*grampc)->opt->Nhor);
    createNumMatrix(&(*grampc)->rws->adj, (*grampc)->param->Nx*(*grampc)->opt->Nhor);
    createNumMatrix(&(*grampc)->rws->dcdx, (*grampc)->param->Nx*((*grampc)->opt->Nhor + 1));
    createNumMatrix(&(*grampc)->rws->u, (*grampc)->param->Nu*(*grampc)->opt->Nhor);
    createNumMatrix(&(*grampc)->rws->uls, (*grampc)->param->Nu*(*grampc)->opt->Nhor);
    createNumMatrix(&(*grampc)->rws->uprev, (*grampc)->param->Nu*(*grampc)->opt->Nhor);
    createNumMatrix(&(*grampc)->rws->gradu, (*grampc)->param->Nu*(*grampc)->opt->Nhor);
    createNumMatrix(&(*grampc)->rws->graduprev, (*grampc)->param->Nu*(*grampc)->opt->Nhor);
    createNumMatrix(&(*grampc)->rws->dcdu, (*grampc)->param->Nu*(*grampc)->opt->Nhor);
    createNumMatrix(&(*grampc)->rws->p, (*grampc)->param->Np);
    createNumMatrix(&(*grampc)->rws->pls, (*grampc)->param->Np);
    createNumMatrix(&(*grampc)->rws->pprev, (*grampc)->param->Np);
    createNumMatrix(&(*grampc)->rws->gradp, (*grampc)->param->Np);
    createNumMatrix(&(*grampc)->rws->gradpprev, (*grampc)->param->Np);
    createNumMatrix(&(*grampc)->rws->dcdp, (*grampc)->param->Np * ((*grampc)->opt->Nhor + 1));
    createNumMatrix(&(*grampc)->rws->mult, (*grampc)->param->Nc*(*grampc)->opt->Nhor);
    createNumMatrix(&(*grampc)->rws->pen, (*grampc)->param->Nc*(*grampc)->opt->Nhor);
    createNumMatrix(&(*grampc)->rws->cfct, (*grampc)->param->Nc*(*grampc)->opt->Nhor);
    createNumMatrix(&(*grampc)->rws->cfctprev, (*grampc)->param->Nc*(*grampc)->opt->Nhor);
    createNumMatrix(&(*grampc)->rws->cfctAbsTol, (*grampc)->param->Nc);
    createNumMatrix(&(*grampc)->rws->rwsScale, 2 * ((*grampc)->param->Nx + (*grampc)->param->Nu + (*grampc)->param->Np));
}

void grampc_free_param(typeGRAMPCparam **param)
{
    free((*param)->x0);
    free((*param)->xdes);
    free((*param)->u0);
    free((*param)->udes);
    free((*param)->umax);
    free((*param)->umin);
    free((*param)->p0);
    free((*param)->pmax);
    free((*param)->pmin);

    free(*param);
}

void grampc_free_opt(typeGRAMPCopt **opt)
{
    free((*opt)->xScale);
    free((*opt)->xOffset);
    free((*opt)->uScale);
    free((*opt)->uOffset);
    free((*opt)->pScale);
    free((*opt)->pOffset);
    free((*opt)->cScale);
    free((*opt)->ConstraintsAbsTol);

    free(*opt);
}

void grampc_free_sol(typeGRAMPCsol **sol)
{
    free((*sol)->xnext);
    free((*sol)->unext);
    free((*sol)->pnext);
    free((*sol)->iter);

    free(*sol);
}

void grampc_free_rws(typeGRAMPCrws **rws)
{
    free((*rws)->t);
    free((*rws)->tls);
    free((*rws)->x);
    free((*rws)->adj);
    free((*rws)->dcdx);
    free((*rws)->u);
    free((*rws)->uls);
    free((*rws)->uprev);
    free((*rws)->gradu);
    free((*rws)->graduprev);
    free((*rws)->dcdu);
    free((*rws)->p);
    free((*rws)->pls);
    free((*rws)->pprev);
    free((*rws)->gradp);
    free((*rws)->gradpprev);
    free((*rws)->dcdp);
    free((*rws)->mult);
    free((*rws)->pen);
    free((*rws)->cfct);
    free((*rws)->cfctprev);
    free((*rws)->cfctAbsTol);
    free((*rws)->lsAdapt);
    free((*rws)->lsExplicit);
    free((*rws)->rwsScale);
    free((*rws)->rwsGeneral);
    free((*rws)->rparRodas);
    free((*rws)->iparRodas);
    free((*rws)->workRodas);
    free((*rws)->iworkRodas);

    free(*rws);
}

void grampc_free(typeGRAMPC **grampc)
{
    grampc_free_param(&(*grampc)->param);
    grampc_free_opt(&(*grampc)->opt);
    grampc_free_sol(&(*grampc)->sol);
    grampc_free_rws(&(*grampc)->rws);

    free(*grampc);
}

#endif /* FIXEDSIZE */
