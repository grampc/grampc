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


#ifndef GRAMPC_MESS_H_
#define GRAMPC_MESS_H_


 /* Required Headers */
#include <stdlib.h>
#include <stdio.h>
#include "grampc_init.h"


/* Macro definitions */
#define GRAMPC_ALLOC_FAILED       "Memory allocation for grampc structure failed.\n"
#define PARAM_ALLOC_FAILED        "Memory allocation for parameters structure failed.\n"
#define OPT_ALLOC_FAILED          "Memory allocation for MPC options failed.\n"
#define SOL_ALLOC_FAILED          "Memory allocation for solution structure failed.\n"
#define RWS_ALLOC_FAILED          "Memory allocation for rws structure failed.\n"
#define RWS_ELEMENT_ALLOC_FAILED  "Memory allocation rws elements failed.\n"
#define INVALID_NO_ELEMENTS       "Input vector has an invalid number of elements.\n"
#define DT_NOT_VALID              "Sampling time dt is not valid. dt must be greater than zero.\n"  
#define THOR_NOT_VALID            "Prediction horizon Thor is not valid. Thor must be greater than the sampling time dt.\n"
#define INVALID_OPTION_NAME       "Invalid option name.\n"
#define INVALID_OPTION_VALUE      "Invalid value for option.\n"
#define INVALID_OPTION_DATATYP    "Invalid datatyp of parameter. \n"
#define INVALID_OPTION_FIXEDSIZE  "This option cannot be changed in the fixed-size mode of GRAMPC.\n"
#define INVALID_PARAM_NAME        "Invalid parameter.\n"
#define INVALID_PARAM_VALUE       "Invalid value for parameter.\n"
#define INVALID_PARAM_DATATYP     "Invalid datatyp of parameter. \n"
#define INVALID_SET_OPERATION     "Invalid setting of option or parameter.\n"
#define INVALID_NX                "Invalid number of states Nx.\n"
#define INVALID_NU                "Invalid number of inputs Nu.\n"
#define INVALID_NP                "Invalid number of parameters Np.\n"
#define INVALID_Nh                "Invalid number of inequality constraints Nh.\n"
#define INVALID_Ng                "Invalid number of equality constraints Ng.\n"
#define INVALID_NgT               "Invalid number of terminal equality constraints NgT.\n"
#define INVALID_NhT               "Invalid number of terminal inequality constraints NhT.\n"


/* status codes */
#define STATUS_NONE                               0
#define STATUS_GRADIENT_CONVERGED                 1
#define STATUS_CONSTRAINTS_CONVERGED              2
#define STATUS_LINESEARCH_MIN                     4
#define STATUS_LINESEARCH_MAX                     8
#define STATUS_LINESEARCH_INIT                   16
#define STATUS_MULTIPLIER_UPDATE                 32
#define STATUS_MULTIPLIER_MAX                    64
#define STATUS_PENALTY_MAX                      128
#define STATUS_INFEASIBLE                       256
#define STATUS_INTEGRATOR_INPUT_NOT_CONSISTENT  512
#define STATUS_INTEGRATOR_MAXSTEPS             1024
#define STATUS_INTEGRATOR_STEPS_TOO_SMALL      2048
#define STATUS_INTEGRATOR_MATRIX_IS_SINGULAR   4096
#define STATUS_INTEGRATOR_H_MIN                8192

/* status levels */
#define STATUS_LEVEL_ERROR  1
#define STATUS_LEVEL_WARN   3
#define STATUS_LEVEL_INFO   7
#define STATUS_LEVEL_DEBUG 15

/* status messages */
#define STATUS_MSG_GRADIENT_CONVERGED               "ConvergenceGradientRelTol satisfied for u, p and T.\n" 
#define STATUS_MSG_CONSTRAINTS_CONVERGED            "ConstraintsAbsTol satisfied for all constraints.\n" 
#define STATUS_MSG_LINESEARCH_MIN                   "Line search used LineSearchMin.\n" 
#define STATUS_MSG_LINESEARCH_MAX                   "Line search used LineSearchMax.\n" 
#define STATUS_MSG_LINESEARCH_INIT                  "Line search used LineSearchInit.\n" 
#define STATUS_MSG_MULTIPLIER_UPDATE                "AugLagUpdateGradientRelTol satisfied for u, p and T, updating multipliers.\n" 
#define STATUS_MSG_MULTIPLIER_MAX                   "Lagrange multiplier reached MultiplierMax.\n" 
#define STATUS_MSG_PENALTY_MAX                      "Penalty parameter reached PenaltyMax.\n" 
#define STATUS_MSG_INFEASIBLE                       "Constraints not improved, problem may be infeasible.\n"
#define STATUS_MSG_INTEGRATOR_INPUT_NOT_CONSISTENT  "Input is not consistent for integrator.\n" 
#define STATUS_MSG_INTEGRATOR_MAXSTEPS              "Integrator needs larger Nmax.\n" 
#define STATUS_MSG_INTEGRATOR_STEPS_TOO_SMALL       "Step size becomes too small for integrator.\n" 
#define STATUS_MSG_INTEGRATOR_MATRIX_IS_SINGULAR    "Integration not successful - Matrix is repeatedly singular.\n" 
#define STATUS_MSG_INTEGRATOR_H_MIN                 "Integrator ruku45 used hmin.\n" 


#ifdef MEXCOMPILE
#include "mex.h"
#define myPrint(x,y)   mexPrintf((x),(y))
#define printError(x)  mexErrMsgTxt((x))
#define printErrorAddString(mess,addstring)  mexPrintf("%s: %s",(addstring),(mess)); mexErrMsgTxt(INVALID_SET_OPERATION)
#define printWarningAddString(mess,addstring)  mexPrintf("%s: %s",(addstring),(mess));
#else 
#define myPrint(x,y)   printf((x),(y))
#define printError(x)  printf("%s",(x)); exit(EXIT_FAILURE)
#define printErrorAddString(mess,addstring)  printf("%s: %s",(addstring),(mess)); exit(EXIT_FAILURE)
#define printWarningAddString(mess,addstring)  printf("%s: %s",(addstring),(mess));
#endif



/* Function definitions */
void grampc_error(const char *mess);
void grampc_error_addstring(const char *mess, const char *addstring);
void grampc_warning_addstring(const char *mess, const char *addstring);
void grampc_error_optname(const char *optname);
void grampc_error_optvalue(const char *optname);
void grampc_error_paramname(const char *paramname);
void grampc_error_paramvalue(const char *paramname);

void print_vector(const char *prefix, ctypeRNum *vector, ctypeInt size);

typeInt grampc_printstatus(ctypeInt status, ctypeInt level);
typeInt print_singleStatus(ctypeInt status, ctypeInt statusmask, const typeChar*message);

#endif /* GRAMPC_MESS_H_ */
