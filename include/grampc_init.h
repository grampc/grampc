/* This file is part of GRAMPC - (https://sourceforge.net/projects/grampc/)
 *
 * GRAMPC -- A software framework for embedded nonlinear model predictive
 * control using a gradient-based augmented Lagrangian approach
 *
 * Copyright (C) 2014-2018 by Tobias Englert, Knut Graichen, Felix Mesmer,
 * Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 * Developed at the Institute of Measurement, Control, and Microtechnology,
 * Ulm University. All rights reserved.
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
 * License along with GRAMPC. If not, see <http://www.gnu.org/licenses/>
 *
 */



#ifndef GRAMPC_INIT_H_
#define GRAMPC_INIT_H_


/* Required Headers */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/* Constants for definitions below */
#define USE_DOUBLE 0
#define USE_FLOAT  1

/* type independent min, max, abs operations */
#ifndef MAX
#define	MAX(a,b)  ((a) > (b) ? (a) : (b))
#endif /* max */
#ifndef MIN
#define	MIN(a,b)  ((a) > (b) ? (b) : (a))
#endif /* min */
#ifndef ABS
#define	ABS(a)    ((a) >= 0 ? (a) : -(a))
#endif /* min */

/* Some definitions */
#define USE_typeRNum   USE_DOUBLE /* USE_FLOAT */
#define typeInt        int
#define typeLInt       int      /* for rodas */
#define typeLogical    long int /* for rodas */
#define typeBoolean    int
#define typeChar       char
#define typeUSERPARAM  void
#define INF    (typeRNum)1e20
#define FWINT  1
#define BWINT -1
#define NALS   3
#define NELS   4

#if USE_typeRNum == USE_FLOAT
#define mxtypeRNum_CLASS  mxSINGLE_CLASS
#define typeRNum          float
#define SS_TYPERNUM       SS_SINGLE
#define EPS               1.1e-8f
#define SQRT(x)           sqrtf(x)
#define POW(x,y)          powf(x,y)
#else
#define mxtypeRNum_CLASS  mxDOUBLE_CLASS
#define typeRNum          double
#define SS_TYPERNUM       SS_DOUBLE
#define EPS               2.2e-16
#define SQRT(x)           sqrt(x)
#define POW(x,y)          pow(x,y)
#endif

#define mxtypeInt_CLASS  mxINT32_CLASS

#define ctypeInt   const int
#define ctypeRNum  const typeRNum

/* encoding of the options */
#define INT_OFF         0
#define INT_ON          1

#define INT_UNIFORM     0
#define INT_NONUNIFORM  1

#define INT_EULER       0
#define INT_MODEULER    1
#define INT_HEUN        2
#define INT_RODAS       3
#define INT_RUKU45      4

#define INT_TRAPZ       0
#define INT_SIMPSON     1

#define INT_ADAPTIVELS  0
#define INT_EXPLS1      1
#define INT_EXPLS2      2

#define INT_EXTPEN      0
#define INT_AUGLAG      1

/* Parameter for adaptive line search fitting */
#define aEPS  1e-5
#define JEPS  1e-6

/* Definition of new datatypes */
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

	typeRNum *x0;
	typeRNum *xdes;

	typeRNum *u0;
	typeRNum *udes;
	typeRNum *umax;
	typeRNum *umin;

	typeRNum *p0;
	typeRNum *pmax;
	typeRNum *pmin;

	typeRNum Thor;
	typeRNum Tmax;
	typeRNum Tmin;

	typeRNum dt;
	typeRNum t0;
} typeGRAMPCparam;

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
	typeInt  *FlagsRodas;

	typeInt  LineSearchType;
	typeInt  LineSearchExpAutoFallback;
	typeRNum LineSearchMax;
	typeRNum LineSearchMin;
	typeRNum LineSearchInit;
	typeRNum LineSearchIntervalFactor;
	typeRNum LineSearchAdaptFactor;
	typeRNum LineSearchIntervalTol;

	typeInt  OptimControl;
	typeInt  OptimParam;
	typeRNum OptimParamLineSearchFactor;
	typeInt  OptimTime;
	typeRNum OptimTimeLineSearchFactor;

	typeInt  ScaleProblem;
	typeRNum *xScale;
	typeRNum *xOffset;
	typeRNum *uScale;
	typeRNum *uOffset;
	typeRNum *pScale;
	typeRNum *pOffset;
	typeRNum TScale;
	typeRNum TOffset;
	typeRNum JScale;
	typeRNum *cScale;

	typeInt  EqualityConstraints;
	typeInt  InequalityConstraints;
	typeInt  TerminalEqualityConstraints;
	typeInt  TerminalInequalityConstraints;
	typeInt  ConstraintsHandling;
	typeRNum *ConstraintsAbsTol;

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

typedef struct
{
	typeRNum *xnext;
	typeRNum *unext;
	typeRNum *pnext;
	typeRNum Tnext;
	typeRNum *J;
	typeRNum cfct;
	typeRNum pen;
	typeInt  *iter;
	typeInt  status;
} typeGRAMPCsol;

typedef struct
{
	typeRNum *t;
	typeRNum *tls;

	/* Nx  */
	typeRNum *x;
	typeRNum *adj;
	typeRNum *dcdx;

	/* Nu  */
	typeRNum *u;
	typeRNum *uls;
	typeRNum *uprev;
	typeRNum *gradu;
	typeRNum *graduprev;
	typeRNum *dcdu;

	/* Np */
	typeRNum *p;
	typeRNum *pls;
	typeRNum *pprev;
	typeRNum *gradp;
	typeRNum *gradpprev;
	typeRNum *dcdp;

	/* 1 */
	typeRNum T;
	typeRNum Tprev;
	typeRNum gradT;
	typeRNum gradTprev;
	typeRNum dcdt;

	/* Nc  */
	typeRNum *mult;
	typeRNum *pen;
	typeRNum *cfct;
	typeRNum *cfctprev;
	typeRNum *cfctAbsTol;

	typeRNum *lsAdapt;
	typeRNum *lsExplicit;
	typeRNum *rwsScale;
	typeInt  lrwsGeneral;
	typeRNum *rwsGeneral;

	/* rodas*/
	typeInt  lworkRodas;
	typeInt  liworkRodas;
	typeRNum *rparRodas;
	typeInt  *iparRodas;
	typeRNum *workRodas;
	typeInt  *iworkRodas;
} typeGRAMPCrws;

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


/* lrwsGeneral */
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
#define LevaluateConstraints (grampc->param->Nc+MAX(grampc->param->Nu,grampc->param->Nx))

/* Definition of functions */
void init_rws_time(const typeGRAMPC *grampc);
void init_rws_controls(const typeGRAMPC *grampc);
void init_rws_parameters(const typeGRAMPC *grampc);
void init_rws_linesearch(const typeGRAMPC *grampc);
void init_rws_multipliers(const typeGRAMPC *grampc);
void init_rws_constraints(const typeGRAMPC *grampc);

void grampc_init(typeGRAMPC **grampc, typeUSERPARAM *userparam);
void grampc_free(typeGRAMPC **grampc);
#endif /* GRAMPC_INIT_H_ */
