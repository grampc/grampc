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
 * File: grampc_init.h
 * Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Date: February 2014
 * Version: v1.0
 *
 * HEADER FILE
 * Initialization file for GRAMPC.
 *
 */


#ifndef GRAMPC_INIT_H_
#define GRAMPC_INIT_H_


/* Required Headers */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>


/* type independent min, max, abs operations */
#ifndef MAX
#define	MAX(a,b)	((a) > (b) ? (a) : (b))
#endif /* max */
#ifndef MIN
#define	MIN(a,b)	((a) > (b) ? (b) : (a))
#endif /* min */
#ifndef ABS
#define ABS(a)  fabs(a)
#endif /* abs */
#ifndef ROUND
#define	ROUND(a)	((a) > (floor(a)+0.5) ? (ceil(a)) : (floor(a)))
#endif /* round */


/* Definition of calloc and free function for general use in c
 * and for use in Matlab mex
 */

/*
#ifdef MEXCOMPILE
#define myCalloc(x,y)								mxCalloc((x),(y))
#define myFree(x)									mxFree((x))
#define myRealloc(x,y)								mxRealloc((x),(y))
#else
#define myCalloc(x,y)								calloc((x),(y))
#define myFree(x)									free((x))
#define myRealloc(x,y)								realloc((x),(y))
#endif
 */
#define myCalloc(x,y)	  calloc((x),(y))
#define myRealloc(x,y)	realloc((x),(y))
#define myFree(x)		    free((x))

/* Some definitions */
#define typeRNum		double
#define typeInt			int
#define typeBoolean		int
#define typeChar		char
#define TRUE			1
#define FALSE			0
#define INF				1e20
#define FWINT			1
#define BWINT			-1
#define NLS 			3


/* Maximal number of characters for option name */
#define NAME_MAXITER			(7  + 1)
#define NAME_SHIFTCONTROL		(12 + 1)
#define NAME_SCALEPROBLEM		(12 + 1)
#define NAME_COSTINTMETHOD		(21 + 1)
#define NAME_INTEGRATOR			(10 + 1)
#define NAME_INTRELTOL			(16 + 1)
#define NAME_INTABSTOL			(16 + 1)
#define NAME_LSTYPE				(14 + 1)
#define NAME_LSMAX				(13 + 1)
#define NAME_LSMIN				(13 + 1)
#define NAME_LSINIT				(14 + 1)
#define NAME_LSFORMULA			(17 + 1)
#define NAME_LSINTFACTOR		(24 + 1)
#define NAME_LSADAPTFACTOR		(21 + 1)
#define NAME_LSINTTOL			(21 + 1)
#define NAME_JACOBIANX		    (9  + 1)
#define NAME_JACOBIANU		    (9  + 1)
#define NAME_INTEGRALCOST		(12 + 1)
#define NAME_FINALCOST		    (9 + 1)
#define NAME_ERROR				(23 + 1)
#define NAME_THOR				(4  + 1)
#define NAME_NHOR				(4  + 1)
#define NAME_DT					(2  + 1)
#define NAME_TK					(2  + 1)
#define NAME_NPCOST				(6  + 1)
#define NAME_NPSYS				(5  + 1)
#define NAME_XK					(2  + 1)
#define NAME_U0					(2  + 1)
#define NAME_XDES				(4  + 1)
#define NAME_UDES				(4  + 1)
#define NAME_PCOST				(5  + 1)
#define NAME_PSYS				(4  + 1)
#define NAME_UMAX				(4  + 1)
#define NAME_UMIN				(4  + 1)
#define NAME_XSCALE				(6  + 1)
#define NAME_XOFFSET			(7  + 1)
#define NAME_USCALE				(6  + 1)
#define NAME_UOFFSET			(7  + 1)
/* Maximum number of characters for option value */
#define VALUE_ONOFF				(3  + 1)
#define VALUE_COSTINTMETHOD		(11 + 1)
#define VALUE_INTEGRATOR		(8  + 1)
#define VALUE_LSTYPE			(9  + 1)
#define VALUE_JACOBIANX		    (10 + 1)
#define VALUE_JACOBIANU		    (10 + 1)


/* Parameter for adaptive line search fitting */
#define aEPS		1e-5
#define JEPS		1e-6


/* Definition of new datatypes */
typedef struct
{
  typeInt MaxIter;
  typeChar ShiftControl[VALUE_ONOFF];
  typeChar ScaleProblem[VALUE_ONOFF];
  typeChar CostIntegrator[VALUE_COSTINTMETHOD];
  typeChar Integrator[VALUE_INTEGRATOR];
  typeRNum IntegratorRelTol;
  typeRNum IntegratorAbsTol;
  typeChar LineSearchType[VALUE_LSTYPE];
  typeRNum LineSearchMax;
  typeRNum LineSearchMin;
  typeRNum LineSearchInit;
  typeRNum LineSearchIntervalFactor;
  typeRNum LineSearchAdaptFactor;
  typeRNum LineSearchIntervalTol;
  typeChar JacobianX[VALUE_JACOBIANX];
  typeChar JacobianU[VALUE_JACOBIANU];
  typeChar IntegralCost[VALUE_ONOFF];
  typeChar FinalCost[VALUE_ONOFF];
} typeGRAMPCopt;

typedef struct
{
  typeInt Nx;
  typeInt Nu;
  typeRNum *xk;
  typeRNum *u0;
  typeRNum *xdes;
  typeRNum *udes;
  typeRNum Thor;
  typeRNum dt;
  typeRNum tk;
  typeInt Nhor;
  typeRNum *pCost;
  typeRNum *pSys;
  typeRNum *umax;
  typeRNum *umin;
  typeRNum *xScale;
  typeRNum *xOffset;
  typeRNum *uScale;
  typeRNum *uOffset;
  typeInt NpSys;
  typeInt NpCost;
} typeGRAMPCparam;

typedef struct
{
  typeRNum *t;
  typeRNum *x;
  typeRNum *adj;
  typeRNum *u;
  typeRNum *dHdu;
  typeRNum *lsAdapt;
  typeRNum *uls;
  typeRNum *lsExplicit;
  typeRNum *uprev;
  typeRNum *dHduprev;
  typeRNum *J;
  typeRNum *rwsScale;
  typeRNum *rwsGradient;
  typeRNum *rwsCostIntegration;
  typeRNum *rwsAdjIntegration;
  typeRNum *rwsIntegration;
} typeGRAMPCrws;

typedef struct
{
  typeRNum *xnext;
  typeRNum *unext;
  typeRNum *J;
} typeGRAMPCsol;

typedef struct
{
  typeGRAMPCopt *opt;
  typeGRAMPCparam *param;
  typeGRAMPCrws *rws;
  typeGRAMPCsol *sol;
} typeGRAMPC;


/* Definition of functions */
void grampc_init(typeGRAMPC **grampc);
void grampc_free(typeGRAMPC **grampc);

#endif /* GRAMPC_INIT_H_ */
