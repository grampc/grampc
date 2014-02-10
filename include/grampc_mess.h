/*
 *
 * This file is part of GRAMPC.
 *
 * GRAMPC - a gradient-based MPC software for real-time applications
 *
 * Copyright (C) 2014 by Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Developed at the Institute of Measurement, Control, and
 * Microtechnology, University of Ulm. All rights reserved.
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
 * License along with GRAMPC. If not, see <http://www.gnu.org/licenses/>.
 *
 */


/*
 *
 * File: grampc_mess.h
 * Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Date: February 2014
 * Version: v1.0
 *
 * HEADER FILE
 * Error printing file for GRAMPC.
 *
 */

#ifndef GRAMPC_MESS_H_
#define GRAMPC_MESS_H_


/* Required Headers */
#include <stdlib.h>
#include <stdio.h>


/* Macro definitions */
#define GRAMPC_ALLOC_FAILED       "Memory allocation for grampc structure failed.\n"
#define PARAM_ALLOC_FAILED			  "Memory allocation for parameters structure failed.\n"
#define SOL_ALLOC_FAILED			    "Memory allocation for solution structure failed.\n"
#define RWS_ALLOC_FAILED			    "Memory allocation for rws structure failed.\n"
#define RWS_ELEMENT_ALLOC_FAILED	"Memory allocation rws elements failed.\n"
#define STATESCALE_VALUE_ZERO		  "At least one scaling value for states is zero.\n"
#define CONTROLSCALE_VALUE_ZERO		"At least one scaling value for controls is zero.\n"
#define OPT_ALLOC_FAILED			    "Memory allocation for MPC options failed.\n"
#define INVALID_NO_ELEMENTS			"Input vector has an invalid number of elements.\n"
#define U0_NOT_DEFINED				"Initpoint u0 not defined.\n"
#define XK_NOT_DEFINED				"Initpoint xk not defined.\n"
#define UDES_NOT_DEFINED			"Setpoint udes not defined.\n"
#define XDES_NOT_DEFINED			"Setpoint xdes not defined.\n"
#define DT_NOT_DEFINED                          "Sampling time dt is not defined.\n"  
#define THOR_NOT_DEFINED                        "Prediction horizon Thor is not defined.\n"
#define NPCOST_ZERO					"Invalid operation. Number of elements NpCost is zero.\n"
#define NPSYS_ZERO					"Invalid operation. Number of elements NpSys is zero.\n"
#define PCOST_NULL					"NpCost unequal zero but pCost is empty.\n"
#define PSYS_NULL					"NpSys unequal zero but pSys is empty.\n"
#define INVALID_OPTION_NAME			"Invalid option name.\n"
#define INVALID_OPTION_VALUE		"Invalid value for option.\n"
#define INVALID_PARAM_NAME			"Invalid parameter.\n"
#define INVALID_PARAM_VALUE			"Invalid value for parameter.\n"
#define INVALID_SET_OPERATION		"Invalid setting of option or parameter.\n"
#define INVALID_NX					"Invalid number of states Nx.\n"
#define INVALID_NU					"Invalid number of states Nu.\n"

#ifdef MEXCOMPILE
#define myPrint(x,y)								mexPrintf((x),(y))
#define printError(x)								mexErrMsgTxt((x))
#define printErrorAddString(mess,addstring)			mexPrintf("%s: %s",(addstring),(mess)); mexErrMsgTxt(INVALID_SET_OPERATION)
#else
#define myPrint(x,y)								printf((x),(y))
#define printError(x)								printf("%s",(x)); exit(EXIT_FAILURE)
#define printErrorAddString(mess,addstring)			printf("%s: %s",(addstring),(mess)); exit(EXIT_FAILURE)
#endif


/* Function definitions */
void grampc_error(char *mess);
void grampc_error_addstring(char *mess, char *addstring);
void grampc_error_optname(char *optname);
void grampc_error_optvalue(char *optname);
void grampc_error_paramname(char *paramname);
void grampc_error_paramvalue(char *paramname);


#endif /* GRAMPC_MESS_H_ */
