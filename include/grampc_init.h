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



#ifndef GRAMPC_INIT_H_
#define GRAMPC_INIT_H_


/* Required Headers */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "grampc_macro.h"
#include "grampc_fixedsize.h"


/* Definition of typeGRAMPCparam containing all problem-specific parameters */
typedef struct
{
    typeInt Nx;
	typeInt Nu;
	typeInt Np;
	typeInt Ng;
	typeInt Nh;
	typeInt NgT;
	typeInt NhT;
	typeInt Nc;

#if defined(FIXEDSIZE) && (NX > 0)
    typeRNum x0[NX];
    typeRNum xdes[NX];
#else
	typeRNum *x0;
	typeRNum *xdes;
#endif

#if defined(FIXEDSIZE) && (NU > 0)
    typeRNum u0[NU];
    typeRNum udes[NU];
    typeRNum umax[NU];
    typeRNum umin[NU];
#else
	typeRNum *u0;
	typeRNum *udes;
	typeRNum *umax;
	typeRNum *umin;
#endif

#if defined(FIXEDSIZE) && (NP > 0)
    typeRNum p0[NP];
    typeRNum pmax[NP];
    typeRNum pmin[NP];
#else
	typeRNum *p0;
	typeRNum *pmax;
	typeRNum *pmin;
#endif

	typeRNum Thor;
	typeRNum Tmax;
	typeRNum Tmin;

	typeRNum dt;
	typeRNum t0;
} typeGRAMPCparam;


/* Definition of typeGRAMPCopt containing all algorithm-related options */
typedef struct
{
	typeInt Nhor;
	typeInt MaxGradIter;
	typeInt MaxMultIter;
	typeInt ShiftControl;

	typeInt TimeDiscretization;

	typeInt IntegralCost;
	typeInt TerminalCost;
	typeInt IntegratorCost;

	typeInt  Integrator;
	typeRNum IntegratorRelTol;
	typeRNum IntegratorAbsTol;
	typeRNum IntegratorMinStepSize;
	typeInt  IntegratorMaxSteps;
    typeInt  FlagsRodas[8];

	typeInt  LineSearchType;
	typeInt  LineSearchExpAutoFallback;
	typeRNum LineSearchMax;
	typeRNum LineSearchMin;
	typeRNum LineSearchInit;
	typeRNum LineSearchAdaptAbsTol;
	typeRNum LineSearchAdaptFactor;
	typeRNum LineSearchIntervalTol;
	typeRNum LineSearchIntervalFactor;
	
	typeInt  OptimControl;
	typeInt  OptimParam;
	typeRNum OptimParamLineSearchFactor;
	typeInt  OptimTime;
	typeRNum OptimTimeLineSearchFactor;

	typeInt  ScaleProblem;
#if defined(FIXEDSIZE) && (NX > 0)
    typeRNum xScale[NX];
    typeRNum xOffset[NX];
#else
	typeRNum *xScale;
	typeRNum *xOffset;
#endif
#if defined(FIXEDSIZE) && (NU > 0)
    typeRNum uScale[NU];
    typeRNum uOffset[NU];
#else
	typeRNum *uScale;
	typeRNum *uOffset;
#endif
#if defined(FIXEDSIZE) && (NP > 0)
    typeRNum pScale[NP];
    typeRNum pOffset[NP];
#else
	typeRNum *pScale;
	typeRNum *pOffset;
#endif
	typeRNum TScale;
	typeRNum TOffset;
	typeRNum JScale;
#if defined(FIXEDSIZE) && (NC > 0)
    typeRNum cScale[NC];
#else
	typeRNum *cScale;
#endif

	typeInt  EqualityConstraints;
	typeInt  InequalityConstraints;
	typeInt  TerminalEqualityConstraints;
	typeInt  TerminalInequalityConstraints;
	typeInt  ConstraintsHandling;
#if defined(FIXEDSIZE) && (NC > 0)
    typeRNum ConstraintsAbsTol[NC];
#else
	typeRNum *ConstraintsAbsTol;
#endif

	typeRNum MultiplierMax;
	typeRNum MultiplierDampingFactor;
	typeRNum PenaltyMax;
	typeRNum PenaltyMin;
	typeRNum PenaltyIncreaseFactor;
	typeRNum PenaltyDecreaseFactor;
	typeRNum PenaltyIncreaseThreshold;
	typeRNum AugLagUpdateGradientRelTol;

	typeInt  ConvergenceCheck;
	typeRNum ConvergenceGradientRelTol;
} typeGRAMPCopt;


/* Definition of typeGRAMPCsol containing the (public) solution data */
typedef struct
{
#if defined(FIXEDSIZE) && (NX > 0)
    typeRNum xnext[NX];
#else
	typeRNum *xnext;
#endif
#if defined(FIXEDSIZE) && (NU > 0)
    typeRNum unext[NU];
#else
    typeRNum *unext;
#endif
#if defined(FIXEDSIZE) && (NP > 0)
    typeRNum pnext[NP];
#else
	typeRNum *pnext;
#endif
	typeRNum Tnext;
    typeRNum J[2];
	typeRNum cfct;
	typeRNum pen;
#if defined(FIXEDSIZE) && (MAXMULTITER > 0)
    typeInt iter[MAXMULTITER];
#else
	typeInt  *iter;
#endif
	typeInt  status;
} typeGRAMPCsol;


/* Definition of typeGRAMPCrws containing the (private) workspace data */
typedef struct
{
#if defined(FIXEDSIZE) && (NHOR > 0)
    typeRNum t[NHOR];
    typeRNum tls[NHOR];
#else
	typeRNum *t;
	typeRNum *tls;
#endif

#if defined(FIXEDSIZE) && (NX > 0)
    typeRNum x[NX*NHOR];
    typeRNum adj[NX*NHOR];
    typeRNum dcdx[NX*(NHOR+1)];
#else
	typeRNum *x;
	typeRNum *adj;
	typeRNum *dcdx;
#endif

#if defined(FIXEDSIZE) && (NU > 0)
    typeRNum u[NU*NHOR];
    typeRNum uls[NU*NHOR];
    typeRNum uprev[NU*NHOR];
    typeRNum gradu[NU*NHOR];
    typeRNum graduprev[NU*(NHOR)];
    typeRNum dcdu[NU*(NHOR+1)];
#else
    typeRNum *u;
    typeRNum *uls;
    typeRNum *uprev;
    typeRNum *gradu;
    typeRNum *graduprev;
    typeRNum *dcdu;
#endif

#if defined(FIXEDSIZE) && (NP > 0)
    typeRNum p[NP];
    typeRNum pls[NP];
    typeRNum pprev[NP];
    typeRNum gradp[NP];
    typeRNum gradpprev[NP];
    typeRNum dcdp[NP*(NHOR+1)];
#else
    typeRNum *p;
    typeRNum *pls;
    typeRNum *pprev;
    typeRNum *gradp;
    typeRNum *gradpprev;
    typeRNum *dcdp;
#endif

	typeRNum T;
	typeRNum Tprev;
	typeRNum gradT;
	typeRNum gradTprev;
	typeRNum dcdt;

#if defined(FIXEDSIZE) && (NC > 0)
    typeRNum mult[NC*NHOR];
    typeRNum pen[NC*NHOR];
    typeRNum cfct[NC*NHOR];
    typeRNum cfctprev[NC*NHOR];
    typeRNum cfctAbsTol[NC];
#else
    typeRNum *mult;
    typeRNum *pen;
    typeRNum *cfct;
    typeRNum *cfctprev;
    typeRNum *cfctAbsTol;
#endif

#if defined(FIXEDSIZE) && (SIZE_LSADAPTIVE > 0)
    typeRNum lsAdapt[SIZE_LSADAPTIVE];
#else
	typeRNum *lsAdapt;
#endif
#if defined(FIXEDSIZE) && (SIZE_LSEXPLICIT > 0)
    typeRNum lsExplicit[SIZE_LSEXPLICIT];
#else
	typeRNum *lsExplicit;
#endif
#if defined(FIXEDSIZE) && (SIZE_RWSSCALE > 0)
    typeRNum rwsScale[SIZE_RWSSCALE];
#else
	typeRNum *rwsScale;
#endif
	typeInt  lrwsGeneral;
#if defined(FIXEDSIZE) && (SIZE_RWSGENERAL > 0)
    typeRNum rwsGeneral[SIZE_RWSGENERAL];
#else
	typeRNum *rwsGeneral;
#endif

	typeInt  lworkRodas;
    typeInt  liworkRodas;
#if defined(FIXEDSIZE) && (INTEGRATOR == INT_RODAS)
    typeRNum rparRodas[SIZE_RPARRODAS];
    typeInt iparRodas[SIZE_IPARRODAS];
    typeRNum workRodas[SIZE_RWORKRODAS];
    typeInt iworkRodas[SIZE_IWORKRODAS];
#else
	typeRNum *rparRodas;
	typeInt  *iparRodas;
	typeRNum *workRodas;
    typeInt  *iworkRodas;
#endif
} typeGRAMPCrws;


/* Defintition of typeGRAMPC containing the structs above and the userparameters */
typedef struct
{
	typeGRAMPCparam *param;
	typeGRAMPCopt *opt;
	typeGRAMPCsol *sol;
    typeGRAMPCrws *rws;
	typeUSERPARAM *userparam;
} typeGRAMPC;


/* Function pointer defintions */
typedef void(*typeffctPtr)(typeRNum *s, ctypeRNum *y, ctypeRNum *t, ctypeRNum *x,
	ctypeRNum *u, ctypeRNum *p_, ctypeRNum *dcdx, const typeGRAMPC *grampc);
typedef void(*typeIntffctPtr)(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum *t,
	ctypeRNum *x, ctypeRNum *u, ctypeRNum *p_, const typeGRAMPC *grampc, const typeffctPtr pfct);

typedef void(*typeInVfctPtr)(typeRNum *s, ctypeRNum *t, ctypeRNum *x, ctypeRNum *u,
	ctypeRNum *p, const typeGRAMPC *grampc);

/* Utility macros for computing size of lrwsGeneral */
#define LWadjsys (grampc->param->Nx)
#define Leuler (grampc->param->Nx)
#define Lmodeuler (5*grampc->param->Nx+grampc->param->Nu+grampc->param->Nc)
#define Lheun (3*grampc->param->Nx)
#define Lruku45 (18*grampc->param->Nx+grampc->param->Nu)
#define Lrodas (2*grampc->param->Nx+grampc->param->Nu)

#define LIntCostSimpson (grampc->param->Nx+grampc->param->Nu+3*grampc->param->Nc+5)
#define LIntCostTrapezoidal (2)
#define Lgradu (2*grampc->param->Nu)
#define Lgradp (3*grampc->param->Np)
#define LgradT (grampc->param->Nx)
#define LevaluateConstraints (grampc->param->Nc+2*(grampc->param->Nx+grampc->param->Np+grampc->param->Nu))

/* Definition of functions */
void init_rws_time(const typeGRAMPC *grampc);
void init_rws_controls(const typeGRAMPC *grampc);
void init_rws_parameters(const typeGRAMPC *grampc);
void init_rws_linesearch(const typeGRAMPC *grampc);
void init_rws_multipliers(const typeGRAMPC *grampc);
void init_rws_constraints(const typeGRAMPC *grampc);

void grampc_init(typeGRAMPC **grampc, typeUSERPARAM *userparam);

#endif /* GRAMPC_INIT_H_ */
