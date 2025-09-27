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

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS �AS
 IS� AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
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

#define TRUE_ (1)
#define FALSE_ (0)

#ifndef DABS
#define DABS(x) ((x) >= 0. ? (x) : -(x))
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


#endif


