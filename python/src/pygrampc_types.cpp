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
#include "pygrampc_types.hpp"

namespace grampc
{
    GrampcParam::GrampcParam()
        : x0(NULL, 0), xdes(NULL, 0), u0(NULL, 0), udes(NULL, 0),
        umax(NULL, 0), umin(NULL, 0), p0(NULL, 0), pmax(NULL, 0), pmin(NULL, 0)
    {}

    GrampcParam::GrampcParam(const typeGRAMPCparam* param)
        : Nx(&param->Nx),
          Nu(&param->Nu),
          Np(&param->Np),
          Ng(&param->Ng),
          Nh(&param->Nh),
          NgT(&param->NgT),
          NhT(&param->NhT),
          Nc(&param->Nc),
          x0(param->x0, *Nx), 
          xdes(param->xdes, *Nx), 
          u0(param->u0, *Nu), 
          udes(param->udes, *Nu),
          umax(param->umax, *Nu), 
          umin(param->umin, *Nu), 
          p0(param->p0, *Np), 
          pmax(param->pmax, *Np), 
          pmin(param->pmin, *Np),
          Thor(&param->Thor),
          Tmax(&param->Tmax),
          Tmin(&param->Tmin),
          dt(&param->dt),
          t0(&param->t0)
    {}

    void GrampcParam::remap_memory(const typeGRAMPC *grampc)
    {
        Nx = &grampc->param->Nx;
        Nu = &grampc->param->Nu;
        Np = &grampc->param->Np;
        Ng = &grampc->param->Ng;
        Nh = &grampc->param->Nh;
        NgT = &grampc->param->NgT;
        NhT = &grampc->param->NhT;
        Nc = &grampc->param->Nc;

        // placement new uses the preallocated memory of the Eigen::Map types on the stack, so no delete has to be called
        // this is the way to go according to the Eigen documentation: https://eigen.tuxfamily.org/dox/classEigen_1_1Map.html
        new (&x0) Eigen::Map<Vector> (grampc->param->x0, grampc->param->Nx);
        new (&xdes) Eigen::Map<Vector>(grampc->param->xdes, grampc->param->Nx);

        new (&u0) Eigen::Map<Vector>(grampc->param->u0, grampc->param->Nu);
        new (&udes) Eigen::Map<Vector>(grampc->param->udes, grampc->param->Nu);
        new (&umax) Eigen::Map<Vector>(grampc->param->umax, grampc->param->Nu);
        new (&umin) Eigen::Map<Vector>(grampc->param->umin, grampc->param->Nu);

        new (&p0) Eigen::Map<Vector>(grampc->param->p0, grampc->param->Np);
        new (&pmax) Eigen::Map<Vector>(grampc->param->pmax, grampc->param->Np);
        new (&pmin) Eigen::Map<Vector>(grampc->param->pmin, grampc->param->Np);

        Thor = &grampc->param->Thor;
        Tmax = &grampc->param->Tmax;
        Tmin = &grampc->param->Tmin;

        dt = &grampc->param->dt;
        t0 = &grampc->param->t0;
    }
    
    // Eigen::Map needs to be explicitly constructed
    GrampcOpt::GrampcOpt()
        : xScale(NULL, 0), xOffset(NULL, 0), uScale(NULL, 0), uOffset(NULL, 0),
        pScale(NULL, 0), pOffset(NULL, 0), cScale(NULL, 0), ConstraintsAbsTol(NULL, 0)
    {
    }

    void GrampcOpt::remap_memory(const typeGRAMPC *grampc)
    {
        Nhor = &grampc->opt->Nhor;
        MaxGradIter = &grampc->opt->MaxGradIter;
        MaxMultIter = &grampc->opt->MaxMultIter;
        ShiftControl = &grampc->opt->ShiftControl;

        TimeDiscretization = &grampc->opt->TimeDiscretization;

        IntegralCost = &grampc->opt->IntegralCost;
        TerminalCost = &grampc->opt->TerminalCost;
        IntegratorCost = &grampc->opt->IntegratorCost;
        Integrator = &grampc->opt->Integrator;
        IntegratorRelTol = &grampc->opt->IntegratorRelTol;
        IntegratorAbsTol = &grampc->opt->IntegratorAbsTol;
        IntegratorMinStepSize = &grampc->opt->IntegratorMinStepSize;
        IntegratorMaxSteps = &grampc->opt->IntegratorMaxSteps;
        FlagsRodas = std::vector<int>(grampc->opt->FlagsRodas, grampc->opt->FlagsRodas + 8);

        LineSearchType = &grampc->opt->LineSearchType;
        LineSearchExpAutoFallback = &grampc->opt->LineSearchExpAutoFallback;
        LineSearchMax = &grampc->opt->LineSearchMax;
        LineSearchMin = &grampc->opt->LineSearchMin;
        LineSearchInit = &grampc->opt->LineSearchInit;
        LineSearchAdaptAbsTol = &grampc->opt->LineSearchAdaptAbsTol;
        LineSearchAdaptFactor = &grampc->opt->LineSearchAdaptFactor;
        LineSearchIntervalTol = &grampc->opt->LineSearchIntervalTol;
        LineSearchIntervalFactor = &grampc->opt->LineSearchIntervalFactor;

        OptimControl = &grampc->opt->OptimControl;
        OptimParam = &grampc->opt->OptimParam;
        OptimParamLineSearchFactor = &grampc->opt->OptimParamLineSearchFactor;
        OptimTime = &grampc->opt->OptimTime;
        OptimTimeLineSearchFactor = &grampc->opt->OptimTimeLineSearchFactor;

        ScaleProblem = &grampc->opt->ScaleProblem;
        
        // placement new uses the preallocated memory of the Eigen::Map types on the stack, so no delete has to be called
        // this is the way to go according to the Eigen documentation: https://eigen.tuxfamily.org/dox/classEigen_1_1Map.html
        new (&xScale) Eigen::Map<Vector>(grampc->opt->xScale, grampc->param->Nx);
        new (&xOffset) Eigen::Map<Vector>(grampc->opt->xOffset, grampc->param->Nx);
        new (&uScale) Eigen::Map<Vector>(grampc->opt->uScale, grampc->param->Nu);
        new (&uOffset) Eigen::Map<Vector>(grampc->opt->uOffset, grampc->param->Nu);
        new (&pScale) Eigen::Map<Vector>(grampc->opt->pScale, grampc->param->Np);
        new (&pOffset) Eigen::Map<Vector>(grampc->opt->pOffset, grampc->param->Np);
        TScale = &grampc->opt->TScale;
        TOffset = &grampc->opt->TOffset;
        JScale = &grampc->opt->JScale;
        new (&cScale) Eigen::Map<Vector>(grampc->opt->cScale, grampc->param->Nc);

        EqualityConstraints = &grampc->opt->EqualityConstraints;
        InequalityConstraints = &grampc->opt->InequalityConstraints;
        TerminalEqualityConstraints = &grampc->opt->TerminalEqualityConstraints;
        TerminalInequalityConstraints = &grampc->opt->TerminalInequalityConstraints;
        ConstraintsHandling = &grampc->opt->ConstraintsHandling;
        new (&ConstraintsAbsTol) Eigen::Map<Vector>(grampc->opt->ConstraintsAbsTol, grampc->param->Nc);

        MultiplierMax = &grampc->opt->MultiplierMax;
        MultiplierDampingFactor = &grampc->opt->MultiplierDampingFactor;
        PenaltyMax = &grampc->opt->PenaltyMax;
        PenaltyMin = &grampc->opt->PenaltyMin;
        PenaltyIncreaseFactor = &grampc->opt->PenaltyIncreaseFactor;
        PenaltyDecreaseFactor = &grampc->opt->PenaltyDecreaseFactor;
        PenaltyIncreaseThreshold = &grampc->opt->PenaltyIncreaseThreshold;
        AugLagUpdateGradientRelTol = &grampc->opt->AugLagUpdateGradientRelTol;

        ConvergenceCheck = &grampc->opt->ConvergenceCheck;
        ConvergenceGradientRelTol = &grampc->opt->ConvergenceGradientRelTol;
    }

    // Eigen::Map needs to be explicitly constructed
    GrampcSol::GrampcSol()
        : xnext(NULL, 0), unext(NULL, 0), pnext(NULL, 0), J(NULL, 0)
    {
    }

    void GrampcSol::remap_memory(const typeGRAMPC *grampc)
    {
        // placement new uses the preallocated memory of the Eigen::Map types on the stack, so no delete has to be called
        // this is the way to go according to the Eigen documentation: https://eigen.tuxfamily.org/dox/classEigen_1_1Map.html
        new (&xnext) Eigen::Map<Vector>(grampc->sol->xnext, grampc->param->Nx);
        new (&unext) Eigen::Map<Vector>(grampc->sol->unext, grampc->param->Nu);
        new (&pnext) Eigen::Map<Vector>(grampc->sol->pnext, grampc->param->Np);
        Tnext = &grampc->sol->Tnext;
        new (&J) Eigen::Map<Vector>((typeRNum*)&grampc->sol->J[0], 2);
        cfct = &grampc->sol->cfct;
        pen = &grampc->sol->pen;
        iter = std::vector<int>(grampc->sol->iter, (grampc->sol->iter) + (grampc->opt->MaxMultIter));
        status = grampc->sol->status;
    }

    // Eigen::Map needs to be explicitly constructed
    GrampcRWS::GrampcRWS()
        : t(NULL, 0), tls(NULL, 0), x(NULL, 0, 0), adj(NULL, 0, 0), dcdx(NULL, 0, 0), u(NULL, 0, 0), uls(NULL, 0, 0), 
        uprev(NULL, 0, 0), gradu(NULL, 0, 0), graduprev(NULL, 0, 0), dcdu(NULL, 0, 0), p(NULL, 0), pls(NULL, 0), 
        pprev(NULL, 0), gradp(NULL, 0), gradpprev(NULL, 0), dcdp(NULL, 0, 0), mult(NULL, 0, 0), pen(NULL, 0, 0), 
        cfct(NULL, 0, 0), cfctprev(NULL, 0, 0), cfctAbsTol(NULL, 0), lsAdapt(NULL, 0), lsExplicit(NULL, 0), 
        rwsScale(NULL, 0), rwsGeneral(NULL, 0), rparRodas(NULL, 0), workRodas(NULL, 0)
    {
    }

    void GrampcRWS::remap_memory(const typeGRAMPC *grampc)
    {
        // placement new uses the preallocated memory of the Eigen::Map types on the stack, so no delete has to be called
        // this is the way to go according to the Eigen documentation: https://eigen.tuxfamily.org/dox/classEigen_1_1Map.html

        // dim = Nhor;
        new (&t) Eigen::Map<Vector>(grampc->rws->t, grampc->opt->Nhor);
        new (&tls) Eigen::Map<Vector>(grampc->rws->tls, grampc->opt->Nhor);

        // dim = Nx x Nhor;
        new (&x) Eigen::Map<Matrix>(grampc->rws->x, grampc->param->Nx, grampc->opt->Nhor);
        new (&adj) Eigen::Map<Matrix>(grampc->rws->adj, grampc->param->Nx, grampc->opt->Nhor);
        new (&dcdx) Eigen::Map<Matrix>(grampc->rws->dcdx, grampc->param->Nx, grampc->opt->Nhor + 1);

        // dim = Nu x Nhor;
        new (&u) Eigen::Map<Matrix>(grampc->rws->u, grampc->param->Nu, grampc->opt->Nhor);
        new (&uls) Eigen::Map<Matrix>(grampc->rws->uls, grampc->param->Nu, grampc->opt->Nhor);
        new (&uprev) Eigen::Map<Matrix>(grampc->rws->uprev, grampc->param->Nu, grampc->opt->Nhor);
        new (&gradu) Eigen::Map<Matrix>(grampc->rws->gradu, grampc->param->Nu, grampc->opt->Nhor);
        new (&graduprev) Eigen::Map<Matrix>(grampc->rws->graduprev, grampc->param->Nu, grampc->opt->Nhor);
        new (&dcdu) Eigen::Map<Matrix>(grampc->rws->dcdu, grampc->param->Nu, grampc->opt->Nhor);

        // dim = Np
        new (&p) Eigen::Map<Vector>(grampc->rws->p, grampc->param->Np);
        new (&pls) Eigen::Map<Vector>(grampc->rws->pls, grampc->param->Np);
        new (&pprev) Eigen::Map<Vector>(grampc->rws->pprev, grampc->param->Np);
        new (&gradp) Eigen::Map<Vector>(grampc->rws->gradp, grampc->param->Np);
        new (&gradpprev) Eigen::Map<Vector>(grampc->rws->gradpprev, grampc->param->Np);
        // dim = Np x (Nhor + 1)
        new (&dcdp) Eigen::Map<Matrix>(grampc->rws->dcdp, grampc->param->Np, grampc->opt->Nhor + 1);

        // dim = 1
        T = &grampc->rws->T;
        Tprev = &grampc->rws->Tprev;
        gradT = &grampc->rws->gradT;
        gradTprev = &grampc->rws->gradTprev;
        dcdt = &grampc->rws->dcdt;

        // dim = Nc x Nhor
        new (&mult) Eigen::Map<Matrix>(grampc->rws->mult, grampc->param->Nc, grampc->opt->Nhor);
        new (&pen) Eigen::Map<Matrix>(grampc->rws->pen, grampc->param->Nc, grampc->opt->Nhor);
        new (&cfct) Eigen::Map<Matrix>(grampc->rws->cfct, grampc->param->Nc, grampc->opt->Nhor);
        new (&cfctprev) Eigen::Map<Matrix>(grampc->rws->cfctprev, grampc->param->Nc, grampc->opt->Nhor);

        // dim = 1 x Nc
        new (&cfctAbsTol) Eigen::Map<Vector>(grampc->rws->cfctAbsTol, grampc->param->Nc);

        // dim = DYN!
        // Only memory for the selected linesearch is allocated
        if (grampc->opt->LineSearchType == INT_ADAPTIVELS) 
        {
            const int dim = 2 * (NALS + 1) * (1 + grampc->opt->MaxGradIter);
            new (&lsAdapt) Eigen::Map<Vector>(grampc->rws->lsAdapt, dim);
        }
        else // if(grampc->opt->LineSearchType == INT_EXPLICIT)
        {
            new (&lsExplicit) Eigen::Map<Vector>(grampc->rws->lsExplicit, NELS);
        }

        // dim = 2*(Nx+Nu+Np)
        const int dim = 2 * (grampc->param->Nx + grampc->param->Nu + grampc->param->Np);
        new (&rwsScale) Eigen::Map<Vector>(grampc->rws->rwsScale, dim);
        lrwsGeneral = &grampc->rws->lrwsGeneral;
        // dim = lrwsGeneral
        new (&rwsGeneral) Eigen::Map<Vector>(grampc->rws->rwsGeneral, grampc->rws->lrwsGeneral);

        lworkRodas = &grampc->rws->lworkRodas;
        liworkRodas = &grampc->rws->liworkRodas;

        if (grampc->opt->Integrator == INT_RODAS)
        {
            // dim = Nhor;
            new (&rparRodas) Eigen::Map<Vector>(grampc->rws->rparRodas, grampc->opt->Nhor);
            // dim = 20;
            iparRodas = std::vector<int>(grampc->rws->iparRodas, (grampc->rws->iparRodas)+20);
            // dim = lworkRodas
            new (&workRodas) Eigen::Map<Vector>(grampc->rws->workRodas, grampc->rws->lworkRodas);
            // dim = liworkRodas
            iworkRodas = std::vector<int>(grampc->rws->iworkRodas, (grampc->rws->iworkRodas) + (grampc->rws->liworkRodas));
        }
    }
}