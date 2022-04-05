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


#ifndef GRAMPC_RUN_H_
#define GRAMPC_RUN_H_

#include "grampc_init.h"
#include "grampc_mess.h"
#include "grampc_util.h"
#include "probfct.h"
#include "euler1.h"
#include "eulermod2.h"
#include "heun2.h"
#include "rodas.h"
#include "ruku45.h"
#include "trapezodial.h"
#include "simpson.h"


 /* main function, run MaxGradIter gradient and MaxMultIter multiplier iterations */
void grampc_run(const typeGRAMPC *grampc);

/* evaluate constraints, evaluate jacobians w.r.t. x, u and p and updates multipliers */
void evaluate_constraints(ctypeRNum *t, ctypeRNum *u, ctypeRNum *p, const typeBoolean evaljac, const typeBoolean updatemultiplier, const typeGRAMPC *grampc);
/* update lagrange multiplier and penalty parameters for (terminal) equality constraints */
void update_multiplier_eqc(typeRNum *mult, typeRNum *pen, ctypeRNum *cfct, typeRNum *cfctprev, ctypeRNum *thresholds, ctypeInt Ncon, typeBoolean converged_grad, const typeGRAMPC *grampc);
/* update lagrange multiplier and penalty parameters for (terminal) inequality constraints */
void update_multiplier_ieqc(typeRNum *mult, typeRNum *pen, ctypeRNum *cfct, typeRNum *cfctprev, ctypeRNum *thresholds, ctypeInt Ncon, typeBoolean converged_grad, const typeGRAMPC *grampc);
/* compute cfct^+ = max(cfct, -mult/pen) necessary for (terminal) inequality constraints */
void update_cfct_for_ieqc(ctypeRNum *mult, ctypeRNum *pen, typeRNum *cfct, ctypeInt Ncon);
/* compute multiplier for constraint jacobians */
void compute_jacobian_multiplier(typeRNum *c, ctypeRNum *mult, ctypeRNum *pen, ctypeRNum *cfct, ctypeInt Ncon);

/* forward integration of system */
void evaluate_sys(ctypeRNum *t, ctypeRNum *u, ctypeRNum *p, const typeGRAMPC *grampc);
/* backward integration of adjoint system */
void evaluate_adjsys(ctypeRNum *t, ctypeRNum *u, ctypeRNum *p, const typeGRAMPC *grampc);
/* single step of system integration */
void Wsys(typeRNum *s, ctypeRNum *x, ctypeRNum *t, ctypeRNum *dummy, ctypeRNum *u,
	ctypeRNum *p_, ctypeRNum *dcdx, const typeGRAMPC *grampc);
/* single step of adjoint system integration */
void Wadjsys(typeRNum *s, ctypeRNum *adj, ctypeRNum *t, ctypeRNum *x,
	ctypeRNum *u, ctypeRNum *p_, ctypeRNum *dcdx, const typeGRAMPC *grampc);

/* evaluate gradient w.r.t. u */
void evaluate_gradu(const typeGRAMPC *grampc);
/* projection of u on admissible set [umin, umax] */
void inputproj(typeRNum *u, const typeGRAMPC *grampc);

/* evaluate gradient w.r.t. p */
void evaluate_gradp(const typeGRAMPC *grampc);
/* single step of param gradient integration */
void WintParam(typeRNum *s, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u,
	ctypeRNum *p_, ctypeRNum *dcdp, const typeGRAMPC *grampc);
/* terminal gradient */
void WtermParam(typeRNum *s, ctypeRNum t, ctypeRNum *x, ctypeRNum *p_, ctypeRNum *dcdp, const typeGRAMPC *grampc);
/* projection of p on admissible set [pmin, pmax] */
void paramproj(typeRNum *p, const typeGRAMPC *grampc);

/* evaluate gradient w.r.t. T */
void evaluate_gradT(const typeGRAMPC *grampc);
/* projection of T on admissible set [Tmin, Tmax] */
void timeproj(typeRNum *T, const typeGRAMPC *grampc);

/* adaptive line search strategy */
void linesearch_adaptive(typeRNum *alpha, ctypeInt igrad, const typeGRAMPC *grampc);
/* explicit line search strategy */
void linesearch_explicit(typeRNum *alpha, const typeGRAMPC *grampc);
/* explicit line search formula */
void update_lsExplicit(typeRNum* lsExplicit, ctypeRNum* a, ctypeRNum* aprev, ctypeRNum* dHda, ctypeRNum* dHdaprev, ctypeInt length, const typeGRAMPC *grampc);

/* evaluate cost function */
void evaluate_cost(typeRNum *s, ctypeRNum *t, ctypeRNum *u, ctypeRNum *p, const typeGRAMPC *grampc);
/* single step of integral cost integration */
void WintCost(typeRNum *s, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p,
	ctypeRNum *mult, ctypeRNum *pen, ctypeRNum *cfct, const typeGRAMPC *grampc);
/* terminal cost */
void WtermCost(typeRNum *s, ctypeRNum t, ctypeRNum *x, ctypeRNum *p,
	ctypeRNum *mult, ctypeRNum *pen, ctypeRNum *cfct, const typeGRAMPC *grampc);

/* test if gradient is converged ||u-uprev|| / ||u|| < rel_tol */
typeBoolean convergence_test_gradient(typeRNum rel_tol, const typeGRAMPC *grampc);
/* test if constraints are converged cfct[i] < abs_tol[i] */
typeBoolean convergence_test_constraints(ctypeRNum *abs_tol, const typeGRAMPC *grampc);

/* shifts trajectory by the sampling time */
void shiftTrajectory(typeRNum* trajectory, ctypeInt Nhor, ctypeInt Nvar, ctypeInt Nshiftrows, ctypeRNum dt, ctypeRNum*t);
/* shortens trajectories from the vetor t to the vector linspace(0,t[Nhor-1]-dt,Nhor) */
void shortenTrajectory(typeRNum *trajectory, ctypeInt Nhor, ctypeInt Nrows, ctypeInt Nshortenrows, ctypeRNum dt, ctypeRNum *t);

#endif /* GRAMPC_RUN_H_ */
