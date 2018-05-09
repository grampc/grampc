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


 /*

 Copyright (c) 2004, Ernst Hairer

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 - Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.

 - Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS
 IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

 /* f2cmod.h  --  Standard Fortran to C header file */



#ifndef F2C_INCLUDE
#define F2C_INCLUDE

#ifdef INTEGER_STAR_8	/* Adjust for typeLInt*8. */
typedef long long longint;		/* system-dependent */
typedef unsigned long long ulongint;	/* system-dependent */
#define qbit_clear(a,b)	((a) & ~((ulongint)1 << (b)))
#define qbit_set(a,b)	((a) |  ((ulongint)1 << (b)))
#endif

#define TRUE_ (1)
#define FALSE_ (0)

/* I/O stuff */
#ifdef f2c_i2
/* for -i2 */
typedef short flag;
typedef short ftnlen;
typedef short ftnint;
#else
typedef long int flag;
typedef long int ftnlen;
typedef long int ftnint;
#endif

/*external read, write*/
typedef struct
{
	flag cierr;
	ftnint ciunit;
	flag ciend;
	char *cifmt;
	ftnint cirec;
} cilist;

#ifndef DABS
#define DABS(x) ((x) >= 0. ? (x) : -(x))
#endif

#ifdef __cplusplus
typedef int /* Unknown procedure type */(*U_fp)(.);
typedef /* Subroutine */ int(*S_fp)(.);
#else
typedef int /* Unknown procedure type */(*U_fp)();
typedef /* Subroutine */ int(*S_fp)();
#endif



/* this code was originally in f2c_extractedfcts.h */

/* ----------------------
pow_dd
---------------------- */
#ifdef KR_headers
typeRNum pow_dd(ap, bp) typeRNum *ap, *bp;
#else
#undef abs
#include "math.h"
#ifdef __cplusplus
extern "C" {
#endif
	typeRNum pow_dd(typeRNum *ap, typeRNum *bp)
#endif
	{
		return(POW(*ap, *bp));
	}
#ifdef __cplusplus
}
#endif

/* ----------------------
d_sign
---------------------- */
#ifdef __cplusplus
extern "C" {
#endif

#ifdef KR_headers
	typeRNum d_sign(a, b) typeRNum *a, *b;
#else
	typeRNum d_sign(typeRNum *a, typeRNum *b)
#endif
	{
		typeRNum x;
		x = (*a >= 0 ? *a : -*a);
		return(*b >= 0 ? x : -x);
	}
#ifdef __cplusplus
}
#endif

/* ----------------------
dummy output functions
---------------------- */
typeLInt e_wsle(void)
{
	return(-1);
}

typeLInt s_wsle(cilist *dummy)
{
	return(-1);
}

typeLInt do_lio(typeLInt *dummy1, typeLInt *dummy2, char *dummy3, ftnlen dummy4)
{
	return(-1);
}

typeLInt e_wsfe(void)
{
	return(-1);
}

typeLInt do_fio(typeLInt *dummy1, char *dummy2, ftnlen dummy3)
{
	return(-1);
}

typeLInt s_wsfe(cilist *dummy)
{
	return(-1);
};



#endif


