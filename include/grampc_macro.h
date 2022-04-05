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
#ifndef GRAMPC_MACRO_H
#define GRAMPC_MACRO_H

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
#endif /* abs */

/* Some definitions */
#define USE_typeRNum   USE_DOUBLE /* USE_FLOAT */
#define typeInt        int
#define typeLInt       int      /* for rodas */
#define typeLogical    long int /* for rodas */
#define typeBoolean    int
#define typeChar       char
#define typeUSERPARAM  void
#define INF    (typeRNum)1e20

/* Integer values corresponding to forward / backward integration */
#define FWINT  1
#define BWINT -1

/* Number of fields for adaptive line search */
#define NALS   3
/* Number of fields for explicit line search */
#define NELS   4

/* Definitions depending on floating-point type */
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
#define ctypeInt         const int
#define ctypeRNum        const typeRNum

/* Integer values corresponding to binary options (off/on) */
#define INT_OFF         0
#define INT_ON          1

/* Integer values corresponding to option TimeDiscretization */
#define INT_UNIFORM     0
#define INT_NONUNIFORM  1

/* Integer values corresponding to option Integrator */
#define INT_EULER       0
#define INT_MODEULER    1
#define INT_HEUN        2
#define INT_RODAS       3
#define INT_RUKU45      4

/* Integer values corresponding to option IntegratorCost */
#define INT_TRAPZ       0
#define INT_SIMPSON     1

/* Integer values corresponding to option LineSearchType */
#define INT_ADAPTIVELS  0
#define INT_EXPLS1      1
#define INT_EXPLS2      2

/* Integer values corresponding to option ConstraintsHandling */
#define INT_EXTPEN      0
#define INT_AUGLAG      1

/* Parameter for adaptive line search fitting */
#define aEPS  1e-5

#endif /* GRAMPC_MACRO_H */
