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



/* ==================================================================
	 rodas_decsol_f2c.h

	 -  merged version of RODAS & DECSOL (algebra fcts. required for RODAS)
	 -> rodas_decsol.f
	 -  C code translation with f2c
	 -  additional modification: abs replaced by DABS (see f2cmod.h)
	 -  all necessary f2c-functions are gathered in f2c_mod.h
	 -> linkage with f2c library not required
	 -  IO routines of f2c are replaced by dummy routines
	 -> no screen output during call of RODAS!

	 Knut Graichen, 2011/10/26
	 ================================================================== */

	 /* rodas_decsol.f -- translated by f2c (version 20060506).
			You must link the resulting object file with libf2c:
		 on Microsoft Windows system, link with libf2c.lib;
		 on Linux or Unix systems, link with ../path/to/libf2c.a -lm
		 or, if you install libf2c.a in a standard place, with -lf2c -lm
		 -- in that order, at the end of the command line, as in
			 cc *.o -lf2c -lm
		 Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

			 http://www.netlib.org/f2c/libf2c.zip
	 */


#ifndef RODAS_DECSOL_F2C_H_
#define RODAS_DECSOL_F2C_H_

	 /* Common Block Declarations */

struct {
	typeLInt mle, mue, mbjac, mbb, mdiag, mdiff, mbdiag;
} linal_;

#define linal_1 linal_

union {
	struct {
		typeRNum xold, hout;
		typeLInt nn;
	} _1;
	struct {
		typeRNum xold, h__;
		typeLInt n;
	} _2;
} conros_;

#define conros_1 (conros_._1)
#define conros_2 (conros_._2)

/* Table of constant values */

static typeLInt c__1 = 1;
static typeLInt c__9 = 9;
static typeLInt c__3 = 3;
static typeLInt c__5 = 5;
static typeRNum c_b361 = 1;
static typeLogical c_false = FALSE_;
static typeLogical c_true = TRUE_;
static typeRNum c_b374 = 0;
static typeRNum c_b378 = .25;

/* ---------------------------------------------------------- */
/* Merged source file of RODAS.F and DECSOL.F and DC_DECSOL.F */

/*     Knut Graichen, 2011/10/25 */
/* ---------------------------------------------------------- */
/* Subroutine */ int dec_(typeLInt *n, typeLInt *ndim, typeRNum *a, typeLInt *
	ip, typeLInt *ier)
{
	/* System generated locals */
	typeLInt a_dim1, a_offset, i__1, i__2, i__3;
	typeRNum d__1, d__2;

	/* Local variables */
	static typeLInt i__, j, k, m;
	static typeRNum t;
	static typeLInt nm1, kp1;

	/* VERSION REAL DOUBLE PRECISION */
	/* ----------------------------------------------------------------------- */
	/*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION. */
	/*  INPUT.. */
	/*     N = ORDER OF MATRIX. */
	/*     NDIM = DECLARED DIMENSION OF ARRAY  A . */
	/*     A = MATRIX TO BE TRIANGULARIZED. */
	/*  OUTPUT.. */
	/*     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U . */
	/*     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
	/*     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW. */
	/*     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O . */
	/*     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE */
	/*           SINGULAR AT STAGE K. */
	/*  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
	/*  DETERM(A) = IP(N)*A(1,1)*A(2,2)*.*A(N,N). */
	/*  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO. */

	/*  REFERENCE.. */
	/*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
	/*     C.A.C.M. 15 (1972), P. 274. */
	/* ----------------------------------------------------------------------- */
			/* Parameter adjustments */
	--ip;
	a_dim1 = *ndim;
	a_offset = 1 + a_dim1;
	a -= a_offset;

	/* Function Body */
	*ier = 0;
	ip[*n] = 1;
	if (*n == 1) {
		goto L70;
	}
	nm1 = *n - 1;
	i__1 = nm1;
	for (k = 1; k <= i__1; ++k) {
		kp1 = k + 1;
		m = k;
		i__2 = *n;
		for (i__ = kp1; i__ <= i__2; ++i__) {
			if ((d__1 = a[i__ + k * a_dim1], DABS(d__1)) > (d__2 = a[m + k *
				a_dim1], DABS(d__2))) {
				m = i__;
			}
			/* L10: */
		}
		ip[k] = m;
		t = a[m + k * a_dim1];
		if (m == k) {
			goto L20;
		}
		ip[*n] = -ip[*n];
		a[m + k * a_dim1] = a[k + k * a_dim1];
		a[k + k * a_dim1] = t;
	L20:
		if (t == 0) {
			goto L80;
		}
		t = (typeRNum) 1. / t;
		i__2 = *n;
		for (i__ = kp1; i__ <= i__2; ++i__) {
			/* L30: */
			a[i__ + k * a_dim1] = -a[i__ + k * a_dim1] * t;
		}
		i__2 = *n;
		for (j = kp1; j <= i__2; ++j) {
			t = a[m + j * a_dim1];
			a[m + j * a_dim1] = a[k + j * a_dim1];
			a[k + j * a_dim1] = t;
			if (t == 0) {
				goto L45;
			}
			i__3 = *n;
			for (i__ = kp1; i__ <= i__3; ++i__) {
				/* L40: */
				a[i__ + j * a_dim1] += a[i__ + k * a_dim1] * t;
			}
		L45:
			/* L50: */
			;
		}
		/* L60: */
	}
L70:
	k = *n;
	if (a[*n + *n * a_dim1] == 0) {
		goto L80;
	}
	return 0;
L80:
	*ier = k;
	ip[*n] = 0;
	return 0;
	/* ----------------------- END OF SUBROUTINE DEC ------------------------- */
} /* dec_ */



/* Subroutine */ int sol_(typeLInt *n, typeLInt *ndim, typeRNum *a,
	typeRNum *b, typeLInt *ip)
{
	/* System generated locals */
	typeLInt a_dim1, a_offset, i__1, i__2;

	/* Local variables */
	static typeLInt i__, k, m;
	static typeRNum t;
	static typeLInt kb, km1, nm1, kp1;

	/* VERSION REAL DOUBLE PRECISION */
	/* ----------------------------------------------------------------------- */
	/*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
	/*  INPUT.. */
	/*    N = ORDER OF MATRIX. */
	/*    NDIM = DECLARED DIMENSION OF ARRAY  A . */
	/*    A = TRIANGULARIZED MATRIX OBTAINED FROM DEC. */
	/*    B = RIGHT HAND SIDE VECTOR. */
	/*    IP = PIVOT VECTOR OBTAINED FROM DEC. */
	/*  DO NOT USE IF DEC HAS SET IER .NE. 0. */
	/*  OUTPUT.. */
	/*    B = SOLUTION VECTOR, X . */
	/* ----------------------------------------------------------------------- */
			/* Parameter adjustments */
	--ip;
	--b;
	a_dim1 = *ndim;
	a_offset = 1 + a_dim1;
	a -= a_offset;

	/* Function Body */
	if (*n == 1) {
		goto L50;
	}
	nm1 = *n - 1;
	i__1 = nm1;
	for (k = 1; k <= i__1; ++k) {
		kp1 = k + 1;
		m = ip[k];
		t = b[m];
		b[m] = b[k];
		b[k] = t;
		i__2 = *n;
		for (i__ = kp1; i__ <= i__2; ++i__) {
			/* L10: */
			b[i__] += a[i__ + k * a_dim1] * t;
		}
		/* L20: */
	}
	i__1 = nm1;
	for (kb = 1; kb <= i__1; ++kb) {
		km1 = *n - kb;
		k = km1 + 1;
		b[k] /= a[k + k * a_dim1];
		t = -b[k];
		i__2 = km1;
		for (i__ = 1; i__ <= i__2; ++i__) {
			/* L30: */
			b[i__] += a[i__ + k * a_dim1] * t;
		}
		/* L40: */
	}
L50:
	b[1] /= a[a_dim1 + 1];
	return 0;
	/* ----------------------- END OF SUBROUTINE SOL ------------------------- */
} /* sol_ */



/* Subroutine */ int dech_(typeLInt *n, typeLInt *ndim, typeRNum *a, typeLInt *
	lb, typeLInt *ip, typeLInt *ier)
{
	/* System generated locals */
	typeLInt a_dim1, a_offset, i__1, i__2, i__3;
	typeRNum d__1, d__2;

	/* Local variables */
	static typeLInt i__, j, k, m;
	static typeRNum t;
	static typeLInt na, nm1, kp1;

	/* VERSION REAL DOUBLE PRECISION */
	/* ----------------------------------------------------------------------- */
	/*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A HESSENBERG */
	/*  MATRIX WITH LOWER BANDWIDTH LB */
	/*  INPUT.. */
	/*     N = ORDER OF MATRIX A. */
	/*     NDIM = DECLARED DIMENSION OF ARRAY  A . */
	/*     A = MATRIX TO BE TRIANGULARIZED. */
	/*     LB = LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED, LB.GE.1). */
	/*  OUTPUT.. */
	/*     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U . */
	/*     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
	/*     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW. */
	/*     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O . */
	/*     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE */
	/*           SINGULAR AT STAGE K. */
	/*  USE  SOLH  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
	/*  DETERM(A) = IP(N)*A(1,1)*A(2,2)*.*A(N,N). */
	/*  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO. */

	/*  REFERENCE.. */
	/*     THIS IS A SLIGHT MODIFICATION OF */
	/*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
	/*     C.A.C.M. 15 (1972), P. 274. */
	/* ----------------------------------------------------------------------- */
			/* Parameter adjustments */
	--ip;
	a_dim1 = *ndim;
	a_offset = 1 + a_dim1;
	a -= a_offset;

	/* Function Body */
	*ier = 0;
	ip[*n] = 1;
	if (*n == 1) {
		goto L70;
	}
	nm1 = *n - 1;
	i__1 = nm1;
	for (k = 1; k <= i__1; ++k) {
		kp1 = k + 1;
		m = k;
		/* Computing MIN */
		i__2 = *n, i__3 = *lb + k;
		na = MIN(i__2, i__3);
		i__2 = na;
		for (i__ = kp1; i__ <= i__2; ++i__) {
			if ((d__1 = a[i__ + k * a_dim1], DABS(d__1)) > (d__2 = a[m + k *
				a_dim1], DABS(d__2))) {
				m = i__;
			}
			/* L10: */
		}
		ip[k] = m;
		t = a[m + k * a_dim1];
		if (m == k) {
			goto L20;
		}
		ip[*n] = -ip[*n];
		a[m + k * a_dim1] = a[k + k * a_dim1];
		a[k + k * a_dim1] = t;
	L20:
		if (t == 0) {
			goto L80;
		}
		t = (typeRNum) 1. / t;
		i__2 = na;
		for (i__ = kp1; i__ <= i__2; ++i__) {
			/* L30: */
			a[i__ + k * a_dim1] = -a[i__ + k * a_dim1] * t;
		}
		i__2 = *n;
		for (j = kp1; j <= i__2; ++j) {
			t = a[m + j * a_dim1];
			a[m + j * a_dim1] = a[k + j * a_dim1];
			a[k + j * a_dim1] = t;
			if (t == 0) {
				goto L45;
			}
			i__3 = na;
			for (i__ = kp1; i__ <= i__3; ++i__) {
				/* L40: */
				a[i__ + j * a_dim1] += a[i__ + k * a_dim1] * t;
			}
		L45:
			/* L50: */
			;
		}
		/* L60: */
	}
L70:
	k = *n;
	if (a[*n + *n * a_dim1] == 0) {
		goto L80;
	}
	return 0;
L80:
	*ier = k;
	ip[*n] = 0;
	return 0;
	/* ----------------------- END OF SUBROUTINE DECH ------------------------ */
} /* dech_ */



/* Subroutine */ int solh_(typeLInt *n, typeLInt *ndim, typeRNum *a, typeLInt *
	lb, typeRNum *b, typeLInt *ip)
{
	/* System generated locals */
	typeLInt a_dim1, a_offset, i__1, i__2, i__3;

	/* Local variables */
	static typeLInt i__, k, m;
	static typeRNum t;
	static typeLInt kb, na, km1, nm1, kp1;

	/* VERSION REAL DOUBLE PRECISION */
	/* ----------------------------------------------------------------------- */
	/*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
	/*  INPUT.. */
	/*    N = ORDER OF MATRIX A. */
	/*    NDIM = DECLARED DIMENSION OF ARRAY  A . */
	/*    A = TRIANGULARIZED MATRIX OBTAINED FROM DECH. */
	/*    LB = LOWER BANDWIDTH OF A. */
	/*    B = RIGHT HAND SIDE VECTOR. */
	/*    IP = PIVOT VECTOR OBTAINED FROM DEC. */
	/*  DO NOT USE IF DECH HAS SET IER .NE. 0. */
	/*  OUTPUT.. */
	/*    B = SOLUTION VECTOR, X . */
	/* ----------------------------------------------------------------------- */
			/* Parameter adjustments */
	--ip;
	--b;
	a_dim1 = *ndim;
	a_offset = 1 + a_dim1;
	a -= a_offset;

	/* Function Body */
	if (*n == 1) {
		goto L50;
	}
	nm1 = *n - 1;
	i__1 = nm1;
	for (k = 1; k <= i__1; ++k) {
		kp1 = k + 1;
		m = ip[k];
		t = b[m];
		b[m] = b[k];
		b[k] = t;
		/* Computing MIN */
		i__2 = *n, i__3 = *lb + k;
		na = MIN(i__2, i__3);
		i__2 = na;
		for (i__ = kp1; i__ <= i__2; ++i__) {
			/* L10: */
			b[i__] += a[i__ + k * a_dim1] * t;
		}
		/* L20: */
	}
	i__1 = nm1;
	for (kb = 1; kb <= i__1; ++kb) {
		km1 = *n - kb;
		k = km1 + 1;
		b[k] /= a[k + k * a_dim1];
		t = -b[k];
		i__2 = km1;
		for (i__ = 1; i__ <= i__2; ++i__) {
			/* L30: */
			b[i__] += a[i__ + k * a_dim1] * t;
		}
		/* L40: */
	}
L50:
	b[1] /= a[a_dim1 + 1];
	return 0;
	/* ----------------------- END OF SUBROUTINE SOLH ------------------------ */
} /* solh_ */


/* Subroutine */ int decc_(typeLInt *n, typeLInt *ndim, typeRNum *ar,
	typeRNum *ai, typeLInt *ip, typeLInt *ier)
{
	/* System generated locals */
	typeLInt ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2, i__3;
	typeRNum d__1, d__2, d__3, d__4;

	/* Local variables */
	static typeLInt i__, j, k, m;
	static typeRNum ti, tr;
	static typeLInt nm1, kp1;
	static typeRNum den, prodi, prodr;

	/* VERSION COMPLEX DOUBLE PRECISION */
	/* ----------------------------------------------------------------------- */
	/*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION */
	/*  ------ MODIFICATION FOR COMPLEX MATRICES -------- */
	/*  INPUT.. */
	/*     N = ORDER OF MATRIX. */
	/*     NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI . */
	/*     (AR, AI) = MATRIX TO BE TRIANGULARIZED. */
	/*  OUTPUT.. */
	/*     AR(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; REAL PART. */
	/*     AI(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; IMAGINARY PART. */
	/*     AR(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
	/*                                                    REAL PART. */
	/*     AI(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
	/*                                                    IMAGINARY PART. */
	/*     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW. */
	/*     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O . */
	/*     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE */
	/*           SINGULAR AT STAGE K. */
	/*  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
	/*  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO. */

	/*  REFERENCE.. */
	/*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
	/*     C.A.C.M. 15 (1972), P. 274. */
	/* ----------------------------------------------------------------------- */
			/* Parameter adjustments */
	--ip;
	ai_dim1 = *ndim;
	ai_offset = 1 + ai_dim1;
	ai -= ai_offset;
	ar_dim1 = *ndim;
	ar_offset = 1 + ar_dim1;
	ar -= ar_offset;

	/* Function Body */
	*ier = 0;
	ip[*n] = 1;
	if (*n == 1) {
		goto L70;
	}
	nm1 = *n - 1;
	i__1 = nm1;
	for (k = 1; k <= i__1; ++k) {
		kp1 = k + 1;
		m = k;
		i__2 = *n;
		for (i__ = kp1; i__ <= i__2; ++i__) {
			if ((d__1 = ar[i__ + k * ar_dim1], DABS(d__1)) + (d__2 = ai[i__ +
				k * ai_dim1], DABS(d__2)) > (d__3 = ar[m + k * ar_dim1],
					DABS(d__3)) + (d__4 = ai[m + k * ai_dim1], DABS(d__4))) {
				m = i__;
			}
			/* L10: */
		}
		ip[k] = m;
		tr = ar[m + k * ar_dim1];
		ti = ai[m + k * ai_dim1];
		if (m == k) {
			goto L20;
		}
		ip[*n] = -ip[*n];
		ar[m + k * ar_dim1] = ar[k + k * ar_dim1];
		ai[m + k * ai_dim1] = ai[k + k * ai_dim1];
		ar[k + k * ar_dim1] = tr;
		ai[k + k * ai_dim1] = ti;
	L20:
		if (DABS(tr) + DABS(ti) == 0) {
			goto L80;
		}
		den = tr * tr + ti * ti;
		tr /= den;
		ti = -ti / den;
		i__2 = *n;
		for (i__ = kp1; i__ <= i__2; ++i__) {
			prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
			prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
			ar[i__ + k * ar_dim1] = -prodr;
			ai[i__ + k * ai_dim1] = -prodi;
			/* L30: */
		}
		i__2 = *n;
		for (j = kp1; j <= i__2; ++j) {
			tr = ar[m + j * ar_dim1];
			ti = ai[m + j * ai_dim1];
			ar[m + j * ar_dim1] = ar[k + j * ar_dim1];
			ai[m + j * ai_dim1] = ai[k + j * ai_dim1];
			ar[k + j * ar_dim1] = tr;
			ai[k + j * ai_dim1] = ti;
			if (DABS(tr) + DABS(ti) == 0) {
				goto L48;
			}
			if (ti == 0) {
				i__3 = *n;
				for (i__ = kp1; i__ <= i__3; ++i__) {
					prodr = ar[i__ + k * ar_dim1] * tr;
					prodi = ai[i__ + k * ai_dim1] * tr;
					ar[i__ + j * ar_dim1] += prodr;
					ai[i__ + j * ai_dim1] += prodi;
					/* L40: */
				}
				goto L48;
			}
			if (tr == 0) {
				i__3 = *n;
				for (i__ = kp1; i__ <= i__3; ++i__) {
					prodr = -ai[i__ + k * ai_dim1] * ti;
					prodi = ar[i__ + k * ar_dim1] * ti;
					ar[i__ + j * ar_dim1] += prodr;
					ai[i__ + j * ai_dim1] += prodi;
					/* L45: */
				}
				goto L48;
			}
			i__3 = *n;
			for (i__ = kp1; i__ <= i__3; ++i__) {
				prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] *
					ti;
				prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] *
					ti;
				ar[i__ + j * ar_dim1] += prodr;
				ai[i__ + j * ai_dim1] += prodi;
				/* L47: */
			}
		L48:
			/* L50: */
			;
		}
		/* L60: */
	}
L70:
	k = *n;
	if ((d__1 = ar[*n + *n * ar_dim1], DABS(d__1)) + (d__2 = ai[*n + *n *
		ai_dim1], DABS(d__2)) == 0) {
		goto L80;
	}
	return 0;
L80:
	*ier = k;
	ip[*n] = 0;
	return 0;
	/* ----------------------- END OF SUBROUTINE DECC ------------------------ */
} /* decc_ */



/* Subroutine */ int solc_(typeLInt *n, typeLInt *ndim, typeRNum *ar,
	typeRNum *ai, typeRNum *br, typeRNum *bi, typeLInt *ip)
{
	/* System generated locals */
	typeLInt ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2;

	/* Local variables */
	static typeLInt i__, k, m, kb;
	static typeRNum ti, tr;
	static typeLInt km1, nm1, kp1;
	static typeRNum den, prodi, prodr;

	/* VERSION COMPLEX DOUBLE PRECISION */
	/* ----------------------------------------------------------------------- */
	/*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
	/*  INPUT.. */
	/*    N = ORDER OF MATRIX. */
	/*    NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI. */
	/*    (AR,AI) = TRIANGULARIZED MATRIX OBTAINED FROM DEC. */
	/*    (BR,BI) = RIGHT HAND SIDE VECTOR. */
	/*    IP = PIVOT VECTOR OBTAINED FROM DEC. */
	/*  DO NOT USE IF DEC HAS SET IER .NE. 0. */
	/*  OUTPUT.. */
	/*    (BR,BI) = SOLUTION VECTOR, X . */
	/* ----------------------------------------------------------------------- */
			/* Parameter adjustments */
	--ip;
	--bi;
	--br;
	ai_dim1 = *ndim;
	ai_offset = 1 + ai_dim1;
	ai -= ai_offset;
	ar_dim1 = *ndim;
	ar_offset = 1 + ar_dim1;
	ar -= ar_offset;

	/* Function Body */
	if (*n == 1) {
		goto L50;
	}
	nm1 = *n - 1;
	i__1 = nm1;
	for (k = 1; k <= i__1; ++k) {
		kp1 = k + 1;
		m = ip[k];
		tr = br[m];
		ti = bi[m];
		br[m] = br[k];
		bi[m] = bi[k];
		br[k] = tr;
		bi[k] = ti;
		i__2 = *n;
		for (i__ = kp1; i__ <= i__2; ++i__) {
			prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
			prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
			br[i__] += prodr;
			bi[i__] += prodi;
			/* L10: */
		}
		/* L20: */
	}
	i__1 = nm1;
	for (kb = 1; kb <= i__1; ++kb) {
		km1 = *n - kb;
		k = km1 + 1;
		den = ar[k + k * ar_dim1] * ar[k + k * ar_dim1] + ai[k + k * ai_dim1]
			* ai[k + k * ai_dim1];
		prodr = br[k] * ar[k + k * ar_dim1] + bi[k] * ai[k + k * ai_dim1];
		prodi = bi[k] * ar[k + k * ar_dim1] - br[k] * ai[k + k * ai_dim1];
		br[k] = prodr / den;
		bi[k] = prodi / den;
		tr = -br[k];
		ti = -bi[k];
		i__2 = km1;
		for (i__ = 1; i__ <= i__2; ++i__) {
			prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
			prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
			br[i__] += prodr;
			bi[i__] += prodi;
			/* L30: */
		}
		/* L40: */
	}
L50:
	den = ar[ar_dim1 + 1] * ar[ar_dim1 + 1] + ai[ai_dim1 + 1] * ai[ai_dim1 +
		1];
	prodr = br[1] * ar[ar_dim1 + 1] + bi[1] * ai[ai_dim1 + 1];
	prodi = bi[1] * ar[ar_dim1 + 1] - br[1] * ai[ai_dim1 + 1];
	br[1] = prodr / den;
	bi[1] = prodi / den;
	return 0;
	/* ----------------------- END OF SUBROUTINE SOLC ------------------------ */
} /* solc_ */



/* Subroutine */ int dechc_(typeLInt *n, typeLInt *ndim, typeRNum *ar,
	typeRNum *ai, typeLInt *lb, typeLInt *ip, typeLInt *ier)
{
	/* System generated locals */
	typeLInt ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2, i__3;
	typeRNum d__1, d__2, d__3, d__4;

	/* Local variables */
	static typeLInt i__, j, k, m, na;
	static typeRNum ti, tr;
	static typeLInt nm1, kp1;
	static typeRNum den, prodi, prodr;

	/* VERSION COMPLEX DOUBLE PRECISION */
	/* ----------------------------------------------------------------------- */
	/*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION */
	/*  ------ MODIFICATION FOR COMPLEX MATRICES -------- */
	/*  INPUT.. */
	/*     N = ORDER OF MATRIX. */
	/*     NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI . */
	/*     (AR, AI) = MATRIX TO BE TRIANGULARIZED. */
	/*  OUTPUT.. */
	/*     AR(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; REAL PART. */
	/*     AI(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; IMAGINARY PART. */
	/*     AR(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
	/*                                                    REAL PART. */
	/*     AI(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
	/*                                                    IMAGINARY PART. */
	/*     LB = LOWER BANDWIDTH OF A (DIAGONAL NOT COUNTED), LB.GE.1. */
	/*     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW. */
	/*     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O . */
	/*     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE */
	/*           SINGULAR AT STAGE K. */
	/*  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
	/*  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO. */

	/*  REFERENCE.. */
	/*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
	/*     C.A.C.M. 15 (1972), P. 274. */
	/* ----------------------------------------------------------------------- */
			/* Parameter adjustments */
	--ip;
	ai_dim1 = *ndim;
	ai_offset = 1 + ai_dim1;
	ai -= ai_offset;
	ar_dim1 = *ndim;
	ar_offset = 1 + ar_dim1;
	ar -= ar_offset;

	/* Function Body */
	*ier = 0;
	ip[*n] = 1;
	if (*lb == 0) {
		goto L70;
	}
	if (*n == 1) {
		goto L70;
	}
	nm1 = *n - 1;
	i__1 = nm1;
	for (k = 1; k <= i__1; ++k) {
		kp1 = k + 1;
		m = k;
		/* Computing MIN */
		i__2 = *n, i__3 = *lb + k;
		na = MIN(i__2, i__3);
		i__2 = na;
		for (i__ = kp1; i__ <= i__2; ++i__) {
			if ((d__1 = ar[i__ + k * ar_dim1], DABS(d__1)) + (d__2 = ai[i__ +
				k * ai_dim1], DABS(d__2)) > (d__3 = ar[m + k * ar_dim1],
					DABS(d__3)) + (d__4 = ai[m + k * ai_dim1], DABS(d__4))) {
				m = i__;
			}
			/* L10: */
		}
		ip[k] = m;
		tr = ar[m + k * ar_dim1];
		ti = ai[m + k * ai_dim1];
		if (m == k) {
			goto L20;
		}
		ip[*n] = -ip[*n];
		ar[m + k * ar_dim1] = ar[k + k * ar_dim1];
		ai[m + k * ai_dim1] = ai[k + k * ai_dim1];
		ar[k + k * ar_dim1] = tr;
		ai[k + k * ai_dim1] = ti;
	L20:
		if (DABS(tr) + DABS(ti) == 0) {
			goto L80;
		}
		den = tr * tr + ti * ti;
		tr /= den;
		ti = -ti / den;
		i__2 = na;
		for (i__ = kp1; i__ <= i__2; ++i__) {
			prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
			prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
			ar[i__ + k * ar_dim1] = -prodr;
			ai[i__ + k * ai_dim1] = -prodi;
			/* L30: */
		}
		i__2 = *n;
		for (j = kp1; j <= i__2; ++j) {
			tr = ar[m + j * ar_dim1];
			ti = ai[m + j * ai_dim1];
			ar[m + j * ar_dim1] = ar[k + j * ar_dim1];
			ai[m + j * ai_dim1] = ai[k + j * ai_dim1];
			ar[k + j * ar_dim1] = tr;
			ai[k + j * ai_dim1] = ti;
			if (DABS(tr) + DABS(ti) == 0) {
				goto L48;
			}
			if (ti == 0) {
				i__3 = na;
				for (i__ = kp1; i__ <= i__3; ++i__) {
					prodr = ar[i__ + k * ar_dim1] * tr;
					prodi = ai[i__ + k * ai_dim1] * tr;
					ar[i__ + j * ar_dim1] += prodr;
					ai[i__ + j * ai_dim1] += prodi;
					/* L40: */
				}
				goto L48;
			}
			if (tr == 0) {
				i__3 = na;
				for (i__ = kp1; i__ <= i__3; ++i__) {
					prodr = -ai[i__ + k * ai_dim1] * ti;
					prodi = ar[i__ + k * ar_dim1] * ti;
					ar[i__ + j * ar_dim1] += prodr;
					ai[i__ + j * ai_dim1] += prodi;
					/* L45: */
				}
				goto L48;
			}
			i__3 = na;
			for (i__ = kp1; i__ <= i__3; ++i__) {
				prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] *
					ti;
				prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] *
					ti;
				ar[i__ + j * ar_dim1] += prodr;
				ai[i__ + j * ai_dim1] += prodi;
				/* L47: */
			}
		L48:
			/* L50: */
			;
		}
		/* L60: */
	}
L70:
	k = *n;
	if ((d__1 = ar[*n + *n * ar_dim1], DABS(d__1)) + (d__2 = ai[*n + *n *
		ai_dim1], DABS(d__2)) == 0) {
		goto L80;
	}
	return 0;
L80:
	*ier = k;
	ip[*n] = 0;
	return 0;
	/* ----------------------- END OF SUBROUTINE DECHC ----------------------- */
} /* dechc_ */



/* Subroutine */ int solhc_(typeLInt *n, typeLInt *ndim, typeRNum *ar,
	typeRNum *ai, typeLInt *lb, typeRNum *br, typeRNum *bi, typeLInt *
	ip)
{
	/* System generated locals */
	typeLInt ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2, i__3, i__4;

	/* Local variables */
	static typeLInt i__, k, m, kb;
	static typeRNum ti, tr;
	static typeLInt km1, nm1, kp1;
	static typeRNum den, prodi, prodr;

	/* VERSION COMPLEX DOUBLE PRECISION */
	/* ----------------------------------------------------------------------- */
	/*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
	/*  INPUT.. */
	/*    N = ORDER OF MATRIX. */
	/*    NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI. */
	/*    (AR,AI) = TRIANGULARIZED MATRIX OBTAINED FROM DEC. */
	/*    (BR,BI) = RIGHT HAND SIDE VECTOR. */
	/*    LB = LOWER BANDWIDTH OF A. */
	/*    IP = PIVOT VECTOR OBTAINED FROM DEC. */
	/*  DO NOT USE IF DEC HAS SET IER .NE. 0. */
	/*  OUTPUT.. */
	/*    (BR,BI) = SOLUTION VECTOR, X . */
	/* ----------------------------------------------------------------------- */
			/* Parameter adjustments */
	--ip;
	--bi;
	--br;
	ai_dim1 = *ndim;
	ai_offset = 1 + ai_dim1;
	ai -= ai_offset;
	ar_dim1 = *ndim;
	ar_offset = 1 + ar_dim1;
	ar -= ar_offset;

	/* Function Body */
	if (*n == 1) {
		goto L50;
	}
	nm1 = *n - 1;
	if (*lb == 0) {
		goto L25;
	}
	i__1 = nm1;
	for (k = 1; k <= i__1; ++k) {
		kp1 = k + 1;
		m = ip[k];
		tr = br[m];
		ti = bi[m];
		br[m] = br[k];
		bi[m] = bi[k];
		br[k] = tr;
		bi[k] = ti;
		/* Computing MIN */
		i__3 = *n, i__4 = *lb + k;
		i__2 = MIN(i__3, i__4);
		for (i__ = kp1; i__ <= i__2; ++i__) {
			prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
			prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
			br[i__] += prodr;
			bi[i__] += prodi;
			/* L10: */
		}
		/* L20: */
	}
L25:
	i__1 = nm1;
	for (kb = 1; kb <= i__1; ++kb) {
		km1 = *n - kb;
		k = km1 + 1;
		den = ar[k + k * ar_dim1] * ar[k + k * ar_dim1] + ai[k + k * ai_dim1]
			* ai[k + k * ai_dim1];
		prodr = br[k] * ar[k + k * ar_dim1] + bi[k] * ai[k + k * ai_dim1];
		prodi = bi[k] * ar[k + k * ar_dim1] - br[k] * ai[k + k * ai_dim1];
		br[k] = prodr / den;
		bi[k] = prodi / den;
		tr = -br[k];
		ti = -bi[k];
		i__2 = km1;
		for (i__ = 1; i__ <= i__2; ++i__) {
			prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
			prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
			br[i__] += prodr;
			bi[i__] += prodi;
			/* L30: */
		}
		/* L40: */
	}
L50:
	den = ar[ar_dim1 + 1] * ar[ar_dim1 + 1] + ai[ai_dim1 + 1] * ai[ai_dim1 +
		1];
	prodr = br[1] * ar[ar_dim1 + 1] + bi[1] * ai[ai_dim1 + 1];
	prodi = bi[1] * ar[ar_dim1 + 1] - br[1] * ai[ai_dim1 + 1];
	br[1] = prodr / den;
	bi[1] = prodi / den;
	return 0;
	/* ----------------------- END OF SUBROUTINE SOLHC ----------------------- */
} /* solhc_ */


/* Subroutine */ int decb_(typeLInt *n, typeLInt *ndim, typeRNum *a, typeLInt *
	ml, typeLInt *mu, typeLInt *ip, typeLInt *ier)
{
	/* System generated locals */
	typeLInt a_dim1, a_offset, i__1, i__2, i__3, i__4;
	typeRNum d__1, d__2;

	/* Local variables */
	static typeLInt i__, j, k, m;
	static typeRNum t;
	static typeLInt md, jk, mm, ju, md1, nm1, kp1, mdl, ijk;

	/* ----------------------------------------------------------------------- */
	/*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED */
	/*  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU */
	/*  INPUT.. */
	/*     N       ORDER OF THE ORIGINAL MATRIX A. */
	/*     NDIM    DECLARED DIMENSION OF ARRAY  A. */
	/*     A       CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS */
	/*                OF THE MATRIX ARE STORED IN THE COLUMNS OF  A  AND */
	/*                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS */
	/*                ML+1 THROUGH 2*ML+MU+1 OF  A. */
	/*     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
	/*     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
	/*  OUTPUT.. */
	/*     A       AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND */
	/*                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT. */
	/*     IP      INDEX VECTOR OF PIVOT INDICES. */
	/*     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O . */
	/*     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE */
	/*                SINGULAR AT STAGE K. */
	/*  USE  SOLB  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
	/*  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*.*A(MD,N)  WITH MD=ML+MU+1. */
	/*  IF IP(N)=O, A IS SINGULAR, SOLB WILL DIVIDE BY ZERO. */

	/*  REFERENCE.. */
	/*     THIS IS A MODIFICATION OF */
	/*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
	/*     C.A.C.M. 15 (1972), P. 274. */
	/* ----------------------------------------------------------------------- */
			/* Parameter adjustments */
	--ip;
	a_dim1 = *ndim;
	a_offset = 1 + a_dim1;
	a -= a_offset;

	/* Function Body */
	*ier = 0;
	ip[*n] = 1;
	md = *ml + *mu + 1;
	md1 = md + 1;
	ju = 0;
	if (*ml == 0) {
		goto L70;
	}
	if (*n == 1) {
		goto L70;
	}
	if (*n < *mu + 2) {
		goto L7;
	}
	i__1 = *n;
	for (j = *mu + 2; j <= i__1; ++j) {
		i__2 = *ml;
		for (i__ = 1; i__ <= i__2; ++i__) {
			/* L5: */
			a[i__ + j * a_dim1] = 0;
		}
	}
L7:
	nm1 = *n - 1;
	i__2 = nm1;
	for (k = 1; k <= i__2; ++k) {
		kp1 = k + 1;
		m = md;
		/* Computing MIN */
		i__1 = *ml, i__3 = *n - k;
		mdl = MIN(i__1, i__3) + md;
		i__1 = mdl;
		for (i__ = md1; i__ <= i__1; ++i__) {
			if ((d__1 = a[i__ + k * a_dim1], DABS(d__1)) > (d__2 = a[m + k *
				a_dim1], DABS(d__2))) {
				m = i__;
			}
			/* L10: */
		}
		ip[k] = m + k - md;
		t = a[m + k * a_dim1];
		if (m == md) {
			goto L20;
		}
		ip[*n] = -ip[*n];
		a[m + k * a_dim1] = a[md + k * a_dim1];
		a[md + k * a_dim1] = t;
	L20:
		if (t == 0) {
			goto L80;
		}
		t = (typeRNum) 1. / t;
		i__1 = mdl;
		for (i__ = md1; i__ <= i__1; ++i__) {
			/* L30: */
			a[i__ + k * a_dim1] = -a[i__ + k * a_dim1] * t;
		}
		/* Computing MIN */
		/* Computing MAX */
		i__3 = ju, i__4 = *mu + ip[k];
		i__1 = MAX(i__3, i__4);
		ju = MIN(i__1, *n);
		mm = md;
		if (ju < kp1) {
			goto L55;
		}
		i__1 = ju;
		for (j = kp1; j <= i__1; ++j) {
			--m;
			--mm;
			t = a[m + j * a_dim1];
			if (m == mm) {
				goto L35;
			}
			a[m + j * a_dim1] = a[mm + j * a_dim1];
			a[mm + j * a_dim1] = t;
		L35:
			if (t == 0) {
				goto L45;
			}
			jk = j - k;
			i__3 = mdl;
			for (i__ = md1; i__ <= i__3; ++i__) {
				ijk = i__ - jk;
				/* L40: */
				a[ijk + j * a_dim1] += a[i__ + k * a_dim1] * t;
			}
		L45:
			/* L50: */
			;
		}
	L55:
		/* L60: */
		;
	}
L70:
	k = *n;
	if (a[md + *n * a_dim1] == 0) {
		goto L80;
	}
	return 0;
L80:
	*ier = k;
	ip[*n] = 0;
	return 0;
	/* ----------------------- END OF SUBROUTINE DECB ------------------------ */
} /* decb_ */



/* Subroutine */ int solb_(typeLInt *n, typeLInt *ndim, typeRNum *a, typeLInt *
	ml, typeLInt *mu, typeRNum *b, typeLInt *ip)
{
	/* System generated locals */
	typeLInt a_dim1, a_offset, i__1, i__2, i__3;

	/* Local variables */
	static typeLInt i__, k, m;
	static typeRNum t;
	static typeLInt kb, md, lm, md1, nm1, imd, kmd, mdl, mdm;

	/* ----------------------------------------------------------------------- */
	/*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
	/*  INPUT.. */
	/*    N      ORDER OF MATRIX A. */
	/*    NDIM   DECLARED DIMENSION OF ARRAY  A . */
	/*    A      TRIANGULARIZED MATRIX OBTAINED FROM DECB. */
	/*    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
	/*    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
	/*    B      RIGHT HAND SIDE VECTOR. */
	/*    IP     PIVOT VECTOR OBTAINED FROM DECB. */
	/*  DO NOT USE IF DECB HAS SET IER .NE. 0. */
	/*  OUTPUT.. */
	/*    B      SOLUTION VECTOR, X . */
	/* ----------------------------------------------------------------------- */
			/* Parameter adjustments */
	--ip;
	--b;
	a_dim1 = *ndim;
	a_offset = 1 + a_dim1;
	a -= a_offset;

	/* Function Body */
	md = *ml + *mu + 1;
	md1 = md + 1;
	mdm = md - 1;
	nm1 = *n - 1;
	if (*ml == 0) {
		goto L25;
	}
	if (*n == 1) {
		goto L50;
	}
	i__1 = nm1;
	for (k = 1; k <= i__1; ++k) {
		m = ip[k];
		t = b[m];
		b[m] = b[k];
		b[k] = t;
		/* Computing MIN */
		i__2 = *ml, i__3 = *n - k;
		mdl = MIN(i__2, i__3) + md;
		i__2 = mdl;
		for (i__ = md1; i__ <= i__2; ++i__) {
			imd = i__ + k - md;
			/* L10: */
			b[imd] += a[i__ + k * a_dim1] * t;
		}
		/* L20: */
	}
L25:
	i__1 = nm1;
	for (kb = 1; kb <= i__1; ++kb) {
		k = *n + 1 - kb;
		b[k] /= a[md + k * a_dim1];
		t = -b[k];
		kmd = md - k;
		/* Computing MAX */
		i__2 = 1, i__3 = kmd + 1;
		lm = MAX(i__2, i__3);
		i__2 = mdm;
		for (i__ = lm; i__ <= i__2; ++i__) {
			imd = i__ - kmd;
			/* L30: */
			b[imd] += a[i__ + k * a_dim1] * t;
		}
		/* L40: */
	}
L50:
	b[1] /= a[md + a_dim1];
	return 0;
	/* ----------------------- END OF SUBROUTINE SOLB ------------------------ */
} /* solb_ */


/* Subroutine */ int decbc_(typeLInt *n, typeLInt *ndim, typeRNum *ar,
	typeRNum *ai, typeLInt *ml, typeLInt *mu, typeLInt *ip, typeLInt *ier)
{
	/* System generated locals */
	typeLInt ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2, i__3, i__4;
	typeRNum d__1, d__2, d__3, d__4;

	/* Local variables */
	static typeLInt i__, j, k, m, md, jk, mm;
	static typeRNum ti;
	static typeLInt ju;
	static typeRNum tr;
	static typeLInt md1, nm1, kp1;
	static typeRNum den;
	static typeLInt mdl, ijk;
	static typeRNum prodi, prodr;

	/* ----------------------------------------------------------------------- */
	/*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED COMPLEX */
	/*  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU */
	/*  INPUT.. */
	/*     N       ORDER OF THE ORIGINAL MATRIX A. */
	/*     NDIM    DECLARED DIMENSION OF ARRAY  A. */
	/*     AR, AI     CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS */
	/*                OF THE MATRIX ARE STORED IN THE COLUMNS OF  AR (REAL */
	/*                PART) AND AI (IMAGINARY PART)  AND */
	/*                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS */
	/*                ML+1 THROUGH 2*ML+MU+1 OF  AR AND AI. */
	/*     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
	/*     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
	/*  OUTPUT.. */
	/*     AR, AI  AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND */
	/*                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT. */
	/*     IP      INDEX VECTOR OF PIVOT INDICES. */
	/*     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O . */
	/*     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE */
	/*                SINGULAR AT STAGE K. */
	/*  USE  SOLBC  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
	/*  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*.*A(MD,N)  WITH MD=ML+MU+1. */
	/*  IF IP(N)=O, A IS SINGULAR, SOLBC WILL DIVIDE BY ZERO. */

	/*  REFERENCE.. */
	/*     THIS IS A MODIFICATION OF */
	/*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
	/*     C.A.C.M. 15 (1972), P. 274. */
	/* ----------------------------------------------------------------------- */
			/* Parameter adjustments */
	--ip;
	ai_dim1 = *ndim;
	ai_offset = 1 + ai_dim1;
	ai -= ai_offset;
	ar_dim1 = *ndim;
	ar_offset = 1 + ar_dim1;
	ar -= ar_offset;

	/* Function Body */
	*ier = 0;
	ip[*n] = 1;
	md = *ml + *mu + 1;
	md1 = md + 1;
	ju = 0;
	if (*ml == 0) {
		goto L70;
	}
	if (*n == 1) {
		goto L70;
	}
	if (*n < *mu + 2) {
		goto L7;
	}
	i__1 = *n;
	for (j = *mu + 2; j <= i__1; ++j) {
		i__2 = *ml;
		for (i__ = 1; i__ <= i__2; ++i__) {
			ar[i__ + j * ar_dim1] = 0;
			ai[i__ + j * ai_dim1] = 0;
			/* L5: */
		}
	}
L7:
	nm1 = *n - 1;
	i__2 = nm1;
	for (k = 1; k <= i__2; ++k) {
		kp1 = k + 1;
		m = md;
		/* Computing MIN */
		i__1 = *ml, i__3 = *n - k;
		mdl = MIN(i__1, i__3) + md;
		i__1 = mdl;
		for (i__ = md1; i__ <= i__1; ++i__) {
			if ((d__1 = ar[i__ + k * ar_dim1], DABS(d__1)) + (d__2 = ai[i__ +
				k * ai_dim1], DABS(d__2)) > (d__3 = ar[m + k * ar_dim1],
					DABS(d__3)) + (d__4 = ai[m + k * ai_dim1], DABS(d__4))) {
				m = i__;
			}
			/* L10: */
		}
		ip[k] = m + k - md;
		tr = ar[m + k * ar_dim1];
		ti = ai[m + k * ai_dim1];
		if (m == md) {
			goto L20;
		}
		ip[*n] = -ip[*n];
		ar[m + k * ar_dim1] = ar[md + k * ar_dim1];
		ai[m + k * ai_dim1] = ai[md + k * ai_dim1];
		ar[md + k * ar_dim1] = tr;
		ai[md + k * ai_dim1] = ti;
	L20:
		if (DABS(tr) + DABS(ti) == 0) {
			goto L80;
		}
		den = tr * tr + ti * ti;
		tr /= den;
		ti = -ti / den;
		i__1 = mdl;
		for (i__ = md1; i__ <= i__1; ++i__) {
			prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
			prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
			ar[i__ + k * ar_dim1] = -prodr;
			ai[i__ + k * ai_dim1] = -prodi;
			/* L30: */
		}
		/* Computing MIN */
		/* Computing MAX */
		i__3 = ju, i__4 = *mu + ip[k];
		i__1 = MAX(i__3, i__4);
		ju = MIN(i__1, *n);
		mm = md;
		if (ju < kp1) {
			goto L55;
		}
		i__1 = ju;
		for (j = kp1; j <= i__1; ++j) {
			--m;
			--mm;
			tr = ar[m + j * ar_dim1];
			ti = ai[m + j * ai_dim1];
			if (m == mm) {
				goto L35;
			}
			ar[m + j * ar_dim1] = ar[mm + j * ar_dim1];
			ai[m + j * ai_dim1] = ai[mm + j * ai_dim1];
			ar[mm + j * ar_dim1] = tr;
			ai[mm + j * ai_dim1] = ti;
		L35:
			if (DABS(tr) + DABS(ti) == 0) {
				goto L48;
			}
			jk = j - k;
			if (ti == 0) {
				i__3 = mdl;
				for (i__ = md1; i__ <= i__3; ++i__) {
					ijk = i__ - jk;
					prodr = ar[i__ + k * ar_dim1] * tr;
					prodi = ai[i__ + k * ai_dim1] * tr;
					ar[ijk + j * ar_dim1] += prodr;
					ai[ijk + j * ai_dim1] += prodi;
					/* L40: */
				}
				goto L48;
			}
			if (tr == 0) {
				i__3 = mdl;
				for (i__ = md1; i__ <= i__3; ++i__) {
					ijk = i__ - jk;
					prodr = -ai[i__ + k * ai_dim1] * ti;
					prodi = ar[i__ + k * ar_dim1] * ti;
					ar[ijk + j * ar_dim1] += prodr;
					ai[ijk + j * ai_dim1] += prodi;
					/* L45: */
				}
				goto L48;
			}
			i__3 = mdl;
			for (i__ = md1; i__ <= i__3; ++i__) {
				ijk = i__ - jk;
				prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] *
					ti;
				prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] *
					ti;
				ar[ijk + j * ar_dim1] += prodr;
				ai[ijk + j * ai_dim1] += prodi;
				/* L47: */
			}
		L48:
			/* L50: */
			;
		}
	L55:
		/* L60: */
		;
	}
L70:
	k = *n;
	if ((d__1 = ar[md + *n * ar_dim1], DABS(d__1)) + (d__2 = ai[md + *n *
		ai_dim1], DABS(d__2)) == 0) {
		goto L80;
	}
	return 0;
L80:
	*ier = k;
	ip[*n] = 0;
	return 0;
	/* ----------------------- END OF SUBROUTINE DECBC ------------------------ */
} /* decbc_ */



/* Subroutine */ int solbc_(typeLInt *n, typeLInt *ndim, typeRNum *ar,
	typeRNum *ai, typeLInt *ml, typeLInt *mu, typeRNum *br, typeRNum *
	bi, typeLInt *ip)
{
	/* System generated locals */
	typeLInt ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2, i__3;

	/* Local variables */
	static typeLInt i__, k, m, kb, md, lm;
	static typeRNum ti, tr;
	static typeLInt md1, nm1;
	static typeRNum den;
	static typeLInt imd, kmd, mdl, mdm;
	static typeRNum prodi, prodr;

	/* ----------------------------------------------------------------------- */
	/*  SOLUTION OF LINEAR SYSTEM, A*X = B , */
	/*                  VERSION BANDED AND COMPLEX-DOUBLE PRECISION. */
	/*  INPUT.. */
	/*    N      ORDER OF MATRIX A. */
	/*    NDIM   DECLARED DIMENSION OF ARRAY  A . */
	/*    AR, AI TRIANGULARIZED MATRIX OBTAINED FROM DECB (REAL AND IMAG. PART). */
	/*    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
	/*    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
	/*    BR, BI RIGHT HAND SIDE VECTOR (REAL AND IMAG. PART). */
	/*    IP     PIVOT VECTOR OBTAINED FROM DECBC. */
	/*  DO NOT USE IF DECB HAS SET IER .NE. 0. */
	/*  OUTPUT.. */
	/*    BR, BI SOLUTION VECTOR, X (REAL AND IMAG. PART). */
	/* ----------------------------------------------------------------------- */
			/* Parameter adjustments */
	--ip;
	--bi;
	--br;
	ai_dim1 = *ndim;
	ai_offset = 1 + ai_dim1;
	ai -= ai_offset;
	ar_dim1 = *ndim;
	ar_offset = 1 + ar_dim1;
	ar -= ar_offset;

	/* Function Body */
	md = *ml + *mu + 1;
	md1 = md + 1;
	mdm = md - 1;
	nm1 = *n - 1;
	if (*ml == 0) {
		goto L25;
	}
	if (*n == 1) {
		goto L50;
	}
	i__1 = nm1;
	for (k = 1; k <= i__1; ++k) {
		m = ip[k];
		tr = br[m];
		ti = bi[m];
		br[m] = br[k];
		bi[m] = bi[k];
		br[k] = tr;
		bi[k] = ti;
		/* Computing MIN */
		i__2 = *ml, i__3 = *n - k;
		mdl = MIN(i__2, i__3) + md;
		i__2 = mdl;
		for (i__ = md1; i__ <= i__2; ++i__) {
			imd = i__ + k - md;
			prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
			prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
			br[imd] += prodr;
			bi[imd] += prodi;
			/* L10: */
		}
		/* L20: */
	}
L25:
	i__1 = nm1;
	for (kb = 1; kb <= i__1; ++kb) {
		k = *n + 1 - kb;
		den = ar[md + k * ar_dim1] * ar[md + k * ar_dim1] + ai[md + k *
			ai_dim1] * ai[md + k * ai_dim1];
		prodr = br[k] * ar[md + k * ar_dim1] + bi[k] * ai[md + k * ai_dim1];
		prodi = bi[k] * ar[md + k * ar_dim1] - br[k] * ai[md + k * ai_dim1];
		br[k] = prodr / den;
		bi[k] = prodi / den;
		tr = -br[k];
		ti = -bi[k];
		kmd = md - k;
		/* Computing MAX */
		i__2 = 1, i__3 = kmd + 1;
		lm = MAX(i__2, i__3);
		i__2 = mdm;
		for (i__ = lm; i__ <= i__2; ++i__) {
			imd = i__ - kmd;
			prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
			prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
			br[imd] += prodr;
			bi[imd] += prodi;
			/* L30: */
		}
		/* L40: */
	}
	den = ar[md + ar_dim1] * ar[md + ar_dim1] + ai[md + ai_dim1] * ai[md +
		ai_dim1];
	prodr = br[1] * ar[md + ar_dim1] + bi[1] * ai[md + ai_dim1];
	prodi = bi[1] * ar[md + ar_dim1] - br[1] * ai[md + ai_dim1];
	br[1] = prodr / den;
	bi[1] = prodi / den;
L50:
	return 0;
	/* ----------------------- END OF SUBROUTINE SOLBC ------------------------ */
} /* solbc_ */



/* Subroutine */ int elmhes_(typeLInt *nm, typeLInt *n, typeLInt *low, typeLInt *
	igh, typeRNum *a, typeLInt *int__)
{
	/* System generated locals */
	typeLInt a_dim1, a_offset, i__1, i__2, i__3;
	typeRNum d__1;

	/* Local variables */
	static typeLInt i__, j, m;
	static typeRNum x, y;
	static typeLInt la, mm1, kp1, mp1;



	/*     this subroutine is a translation of the algol procedure elmhes, */
	/*     num. math. 12, 349-368(1968) by martin and wilkinson. */
	/*     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971). */

	/*     given a real general matrix, this subroutine */
	/*     reduces a submatrix situated in rows and columns */
	/*     low through igh to upper hessenberg form by */
	/*     stabilized elementary similarity transformations. */

	/*     on input: */

	/*      nm must be set to the row dimension of two-dimensional */
	/*        array parameters as declared in the calling program */
	/*        dimension statement; */

	/*      n is the order of the matrix; */

	/*      low and igh are integers determined by the balancing */
	/*        subroutine  balanc.      if  balanc  has not been used, */
	/*        set low=1, igh=n; */

	/*      a contains the input matrix. */

	/*     on output: */

	/*      a contains the hessenberg matrix.  the multipliers */
	/*        which were used in the reduction are stored in the */
	/*        remaining triangle under the hessenberg matrix; */

	/*      int contains information on the rows and columns */
	/*        interchanged in the reduction. */
	/*        only elements low through igh are used. */

	/*     questions and comments should be directed to b. s. garbow, */
	/*     applied mathematics division, argonne national laboratory */

	/*     ------------------------------------------------------------------ */

			/* Parameter adjustments */
	a_dim1 = *nm;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--int__;

	/* Function Body */
	la = *igh - 1;
	kp1 = *low + 1;
	if (la < kp1) {
		goto L200;
	}

	i__1 = la;
	for (m = kp1; m <= i__1; ++m) {
		mm1 = m - 1;
		x = 0;
		i__ = m;

		i__2 = *igh;
		for (j = m; j <= i__2; ++j) {
			if ((d__1 = a[j + mm1 * a_dim1], DABS(d__1)) <= DABS(x)) {
				goto L100;
			}
			x = a[j + mm1 * a_dim1];
			i__ = j;
		L100:
			;
		}

		int__[m] = i__;
		if (i__ == m) {
			goto L130;
		}
		/*    :::::::::: interchange rows and columns of a :::::::::: */
		i__2 = *n;
		for (j = mm1; j <= i__2; ++j) {
			y = a[i__ + j * a_dim1];
			a[i__ + j * a_dim1] = a[m + j * a_dim1];
			a[m + j * a_dim1] = y;
			/* L110: */
		}

		i__2 = *igh;
		for (j = 1; j <= i__2; ++j) {
			y = a[j + i__ * a_dim1];
			a[j + i__ * a_dim1] = a[j + m * a_dim1];
			a[j + m * a_dim1] = y;
			/* L120: */
		}
		/*    :::::::::: end interchange :::::::::: */
	L130:
		if (x == 0) {
			goto L180;
		}
		mp1 = m + 1;

		i__2 = *igh;
		for (i__ = mp1; i__ <= i__2; ++i__) {
			y = a[i__ + mm1 * a_dim1];
			if (y == 0) {
				goto L160;
			}
			y /= x;
			a[i__ + mm1 * a_dim1] = y;

			i__3 = *n;
			for (j = m; j <= i__3; ++j) {
				/* L140: */
				a[i__ + j * a_dim1] -= y * a[m + j * a_dim1];
			}

			i__3 = *igh;
			for (j = 1; j <= i__3; ++j) {
				/* L150: */
				a[j + m * a_dim1] += y * a[j + i__ * a_dim1];
			}

		L160:
			;
		}

	L180:
		;
	}

L200:
	return 0;
	/*    :::::::::: last card of elmhes :::::::::: */
} /* elmhes_ */

/* ****************************************** */
/*     VERSION OF SEPTEMBER 18, 1995 */
/* ****************************************** */

/* Subroutine */ int decomr_(typeLInt *n, typeRNum *fjac, typeLInt *ldjac,
	typeRNum *fmas, typeLInt *ldmas, typeLInt *mlmas, typeLInt *mumas,
	typeLInt *m1, typeLInt *m2, typeLInt *nm1, typeRNum *fac1, typeRNum *
	e1, typeLInt *lde1, typeLInt *ip1, typeLInt *ier, typeLInt *ijob, typeLogical
	*calhes, typeLInt *iphes)
{
	/* System generated locals */
	typeLInt fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1,
		e1_offset, i__1, i__2, i__3, i__4, i__5, i__6;

	/* Local variables */
	static typeLInt i__, j, k, j1, ib, mm, jm1;
	extern /* Subroutine */ int dec_(typeLInt *, typeLInt *, typeRNum *,
		typeLInt *, typeLInt *);
	static typeRNum sum;
	extern /* Subroutine */ int decb_(typeLInt *, typeLInt *, typeRNum *,
		typeLInt *, typeLInt *, typeLInt *, typeLInt *), dech_(typeLInt *,
			typeLInt *, typeRNum *, typeLInt *, typeLInt *, typeLInt *),
		elmhes_(typeLInt *, typeLInt *, typeLInt *, typeLInt *, typeRNum *,
			typeLInt *);


	/* Parameter adjustments */
	--iphes;
	fjac_dim1 = *ldjac;
	fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;
	--ip1;
	fmas_dim1 = *ldmas;
	fmas_offset = 1 + fmas_dim1;
	fmas -= fmas_offset;
	e1_dim1 = *lde1;
	e1_offset = 1 + e1_dim1;
	e1 -= e1_offset;

	/* Function Body */
	switch (*ijob) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	case 4:  goto L4;
	case 5:  goto L5;
	case 6:  goto L6;
	case 7:  goto L7;
	case 8:  goto L55;
	case 9:  goto L55;
	case 10:  goto L55;
	case 11:  goto L11;
	case 12:  goto L12;
	case 13:  goto L13;
	case 14:  goto L14;
	case 15:  goto L15;
	}

	/* ----------------------------------------------------------- */

L1:
	/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
			e1[i__ + j * e1_dim1] = -fjac[i__ + j * fjac_dim1];
		}
		e1[j + j * e1_dim1] += *fac1;
	}
	dec_(n, lde1, &e1[e1_offset], &ip1[1], ier);
	return 0;

	/* ----------------------------------------------------------- */

L11:
	/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
	i__1 = *nm1;
	for (j = 1; j <= i__1; ++j) {
		jm1 = j + *m1;
		i__2 = *nm1;
		for (i__ = 1; i__ <= i__2; ++i__) {
			e1[i__ + j * e1_dim1] = -fjac[i__ + jm1 * fjac_dim1];
		}
		e1[j + j * e1_dim1] += *fac1;
	}
L45:
	mm = *m1 / *m2;
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
		i__2 = *nm1;
		for (i__ = 1; i__ <= i__2; ++i__) {
			sum = 0;
			i__3 = mm - 1;
			for (k = 0; k <= i__3; ++k) {
				sum = (sum + fjac[i__ + (j + k * *m2) * fjac_dim1]) / *fac1;
			}
			e1[i__ + j * e1_dim1] -= sum;
		}
	}
	dec_(nm1, lde1, &e1[e1_offset], &ip1[1], ier);
	return 0;

	/* ----------------------------------------------------------- */

L2:
	/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		i__2 = linal_1.mbjac;
		for (i__ = 1; i__ <= i__2; ++i__) {
			e1[i__ + linal_1.mle + j * e1_dim1] = -fjac[i__ + j * fjac_dim1];
		}
		e1[linal_1.mdiag + j * e1_dim1] += *fac1;
	}
	decb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &ip1[1], ier);
	return 0;

	/* ----------------------------------------------------------- */

L12:
	/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
	i__1 = *nm1;
	for (j = 1; j <= i__1; ++j) {
		jm1 = j + *m1;
		i__2 = linal_1.mbjac;
		for (i__ = 1; i__ <= i__2; ++i__) {
			e1[i__ + linal_1.mle + j * e1_dim1] = -fjac[i__ + jm1 * fjac_dim1]
				;
		}
		e1[linal_1.mdiag + j * e1_dim1] += *fac1;
	}
L46:
	mm = *m1 / *m2;
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
		i__2 = linal_1.mbjac;
		for (i__ = 1; i__ <= i__2; ++i__) {
			sum = 0;
			i__3 = mm - 1;
			for (k = 0; k <= i__3; ++k) {
				sum = (sum + fjac[i__ + (j + k * *m2) * fjac_dim1]) / *fac1;
			}
			e1[i__ + linal_1.mle + j * e1_dim1] -= sum;
		}
	}
	decb_(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &ip1[1], ier)
		;
	return 0;

	/* ----------------------------------------------------------- */

L3:
	/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
			e1[i__ + j * e1_dim1] = -fjac[i__ + j * fjac_dim1];
		}
		/* Computing MAX */
		i__2 = 1, i__3 = j - *mumas;
		/* Computing MIN */
		i__5 = *n, i__6 = j + *mlmas;
		i__4 = MIN(i__5, i__6);
		for (i__ = MAX(i__2, i__3); i__ <= i__4; ++i__) {
			e1[i__ + j * e1_dim1] += *fac1 * fmas[i__ - j + linal_1.mbdiag +
				j * fmas_dim1];
		}
	}
	dec_(n, lde1, &e1[e1_offset], &ip1[1], ier);
	return 0;

	/* ----------------------------------------------------------- */

L13:
	/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
	i__1 = *nm1;
	for (j = 1; j <= i__1; ++j) {
		jm1 = j + *m1;
		i__4 = *nm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
			e1[i__ + j * e1_dim1] = -fjac[i__ + jm1 * fjac_dim1];
		}
		/* Computing MAX */
		i__4 = 1, i__2 = j - *mumas;
		/* Computing MIN */
		i__5 = *nm1, i__6 = j + *mlmas;
		i__3 = MIN(i__5, i__6);
		for (i__ = MAX(i__4, i__2); i__ <= i__3; ++i__) {
			e1[i__ + j * e1_dim1] += *fac1 * fmas[i__ - j + linal_1.mbdiag +
				j * fmas_dim1];
		}
	}
	goto L45;

	/* ----------------------------------------------------------- */

L4:
	/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		i__3 = linal_1.mbjac;
		for (i__ = 1; i__ <= i__3; ++i__) {
			e1[i__ + linal_1.mle + j * e1_dim1] = -fjac[i__ + j * fjac_dim1];
		}
		i__3 = linal_1.mbb;
		for (i__ = 1; i__ <= i__3; ++i__) {
			ib = i__ + linal_1.mdiff;
			e1[ib + j * e1_dim1] += *fac1 * fmas[i__ + j * fmas_dim1];
		}
	}
	decb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &ip1[1], ier);
	return 0;

	/* ----------------------------------------------------------- */

L14:
	/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER */
	i__1 = *nm1;
	for (j = 1; j <= i__1; ++j) {
		jm1 = j + *m1;
		i__3 = linal_1.mbjac;
		for (i__ = 1; i__ <= i__3; ++i__) {
			e1[i__ + linal_1.mle + j * e1_dim1] = -fjac[i__ + jm1 * fjac_dim1]
				;
		}
		i__3 = linal_1.mbb;
		for (i__ = 1; i__ <= i__3; ++i__) {
			ib = i__ + linal_1.mdiff;
			e1[ib + j * e1_dim1] += *fac1 * fmas[i__ + j * fmas_dim1];
		}
	}
	goto L46;

	/* ----------------------------------------------------------- */

L5:
	/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		i__3 = *n;
		for (i__ = 1; i__ <= i__3; ++i__) {
			e1[i__ + j * e1_dim1] = fmas[i__ + j * fmas_dim1] * *fac1 - fjac[
				i__ + j * fjac_dim1];
		}
	}
	dec_(n, lde1, &e1[e1_offset], &ip1[1], ier);
	return 0;

	/* ----------------------------------------------------------- */

L15:
	/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
	i__1 = *nm1;
	for (j = 1; j <= i__1; ++j) {
		jm1 = j + *m1;
		i__3 = *nm1;
		for (i__ = 1; i__ <= i__3; ++i__) {
			e1[i__ + j * e1_dim1] = fmas[i__ + j * fmas_dim1] * *fac1 - fjac[
				i__ + jm1 * fjac_dim1];
		}
	}
	goto L45;

	/* ----------------------------------------------------------- */

L6:
	/* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
	/* ---  THIS OPTION IS NOT PROVIDED */
	return 0;

	/* ----------------------------------------------------------- */

L7:
	/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
	if (*calhes) {
		elmhes_(ldjac, n, &c__1, n, &fjac[fjac_offset], &iphes[1]);
	}
	*calhes = FALSE_;
	i__1 = *n - 1;
	for (j = 1; j <= i__1; ++j) {
		j1 = j + 1;
		e1[j1 + j * e1_dim1] = -fjac[j1 + j * fjac_dim1];
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		i__3 = j;
		for (i__ = 1; i__ <= i__3; ++i__) {
			e1[i__ + j * e1_dim1] = -fjac[i__ + j * fjac_dim1];
		}
		e1[j + j * e1_dim1] += *fac1;
	}
	dech_(n, lde1, &e1[e1_offset], &c__1, &ip1[1], ier);
	return 0;

	/* ----------------------------------------------------------- */

L55:
	return 0;
} /* decomr_ */


/*     END OF SUBROUTINE DECOMR */

/* *********************************************************** */

/* Subroutine */ int decomc_(typeLInt *n, typeRNum *fjac, typeLInt *ldjac,
	typeRNum *fmas, typeLInt *ldmas, typeLInt *mlmas, typeLInt *mumas,
	typeLInt *m1, typeLInt *m2, typeLInt *nm1, typeRNum *alphn, typeRNum
	*betan, typeRNum *e2r, typeRNum *e2i, typeLInt *lde1, typeLInt *ip2,
	typeLInt *ier, typeLInt *ijob)
{
	/* System generated locals */
	typeLInt fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e2r_dim1,
		e2r_offset, e2i_dim1, e2i_offset, i__1, i__2, i__3, i__4, i__5,
		i__6;
	typeRNum d__1, d__2;

	/* Local variables */
	static typeLInt i__, j, k, j1;
	static typeRNum bb;
	static typeLInt ib, mm, jm1;
	static typeRNum bet, alp;
	extern /* Subroutine */ int decc_(typeLInt *, typeLInt *, typeRNum *,
		typeRNum *, typeLInt *, typeLInt *);
	static typeRNum ffma, abno;
	static typeLInt imle;
	static typeRNum sumi, sumr, sums;
	extern /* Subroutine */ int decbc_(typeLInt *, typeLInt *, typeRNum *,
		typeRNum *, typeLInt *, typeLInt *, typeLInt *, typeLInt *), dechc_(
			typeLInt *, typeLInt *, typeRNum *, typeRNum *, typeLInt *,
			typeLInt *, typeLInt *);


	/* Parameter adjustments */
	fjac_dim1 = *ldjac;
	fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;
	--ip2;
	fmas_dim1 = *ldmas;
	fmas_offset = 1 + fmas_dim1;
	fmas -= fmas_offset;
	e2i_dim1 = *lde1;
	e2i_offset = 1 + e2i_dim1;
	e2i -= e2i_offset;
	e2r_dim1 = *lde1;
	e2r_offset = 1 + e2r_dim1;
	e2r -= e2r_offset;

	/* Function Body */
	switch (*ijob) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	case 4:  goto L4;
	case 5:  goto L5;
	case 6:  goto L6;
	case 7:  goto L7;
	case 8:  goto L55;
	case 9:  goto L55;
	case 10:  goto L55;
	case 11:  goto L11;
	case 12:  goto L12;
	case 13:  goto L13;
	case 14:  goto L14;
	case 15:  goto L15;
	}

	/* ----------------------------------------------------------- */

L1:
	/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
			e2r[i__ + j * e2r_dim1] = -fjac[i__ + j * fjac_dim1];
			e2i[i__ + j * e2i_dim1] = 0;
		}
		e2r[j + j * e2r_dim1] += *alphn;
		e2i[j + j * e2i_dim1] = *betan;
	}
	decc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &ip2[1], ier);
	return 0;

	/* ----------------------------------------------------------- */

L11:
	/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
	i__1 = *nm1;
	for (j = 1; j <= i__1; ++j) {
		jm1 = j + *m1;
		i__2 = *nm1;
		for (i__ = 1; i__ <= i__2; ++i__) {
			e2r[i__ + j * e2r_dim1] = -fjac[i__ + jm1 * fjac_dim1];
			e2i[i__ + j * e2i_dim1] = 0;
		}
		e2r[j + j * e2r_dim1] += *alphn;
		e2i[j + j * e2i_dim1] = *betan;
	}
L45:
	mm = *m1 / *m2;
	/* Computing 2nd power */
	d__1 = *alphn;
	/* Computing 2nd power */
	d__2 = *betan;
	abno = d__1 * d__1 + d__2 * d__2;
	alp = *alphn / abno;
	bet = *betan / abno;
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
		i__2 = *nm1;
		for (i__ = 1; i__ <= i__2; ++i__) {
			sumr = 0;
			sumi = 0;
			i__3 = mm - 1;
			for (k = 0; k <= i__3; ++k) {
				sums = sumr + fjac[i__ + (j + k * *m2) * fjac_dim1];
				sumr = sums * alp + sumi * bet;
				sumi = sumi * alp - sums * bet;
			}
			e2r[i__ + j * e2r_dim1] -= sumr;
			e2i[i__ + j * e2i_dim1] -= sumi;
		}
	}
	decc_(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &ip2[1], ier);
	return 0;

	/* ----------------------------------------------------------- */

L2:
	/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		i__2 = linal_1.mbjac;
		for (i__ = 1; i__ <= i__2; ++i__) {
			imle = i__ + linal_1.mle;
			e2r[imle + j * e2r_dim1] = -fjac[i__ + j * fjac_dim1];
			e2i[imle + j * e2i_dim1] = 0;
		}
		e2r[linal_1.mdiag + j * e2r_dim1] += *alphn;
		e2i[linal_1.mdiag + j * e2i_dim1] = *betan;
	}
	decbc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
		linal_1.mue, &ip2[1], ier);
	return 0;

	/* ----------------------------------------------------------- */

L12:
	/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
	i__1 = *nm1;
	for (j = 1; j <= i__1; ++j) {
		jm1 = j + *m1;
		i__2 = linal_1.mbjac;
		for (i__ = 1; i__ <= i__2; ++i__) {
			e2r[i__ + linal_1.mle + j * e2r_dim1] = -fjac[i__ + jm1 *
				fjac_dim1];
			e2i[i__ + linal_1.mle + j * e2i_dim1] = 0;
		}
		e2r[linal_1.mdiag + j * e2r_dim1] += *alphn;
		e2i[linal_1.mdiag + j * e2i_dim1] += *betan;
	}
L46:
	mm = *m1 / *m2;
	/* Computing 2nd power */
	d__1 = *alphn;
	/* Computing 2nd power */
	d__2 = *betan;
	abno = d__1 * d__1 + d__2 * d__2;
	alp = *alphn / abno;
	bet = *betan / abno;
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
		i__2 = linal_1.mbjac;
		for (i__ = 1; i__ <= i__2; ++i__) {
			sumr = 0;
			sumi = 0;
			i__3 = mm - 1;
			for (k = 0; k <= i__3; ++k) {
				sums = sumr + fjac[i__ + (j + k * *m2) * fjac_dim1];
				sumr = sums * alp + sumi * bet;
				sumi = sumi * alp - sums * bet;
			}
			imle = i__ + linal_1.mle;
			e2r[imle + j * e2r_dim1] -= sumr;
			e2i[imle + j * e2i_dim1] -= sumi;
		}
	}
	decbc_(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
		linal_1.mue, &ip2[1], ier);
	return 0;

	/* ----------------------------------------------------------- */

L3:
	/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
			e2r[i__ + j * e2r_dim1] = -fjac[i__ + j * fjac_dim1];
			e2i[i__ + j * e2i_dim1] = 0;
		}
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		/* Computing MAX */
		i__2 = 1, i__3 = j - *mumas;
		/* Computing MIN */
		i__5 = *n, i__6 = j + *mlmas;
		i__4 = MIN(i__5, i__6);
		for (i__ = MAX(i__2, i__3); i__ <= i__4; ++i__) {
			bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
			e2r[i__ + j * e2r_dim1] += *alphn * bb;
			e2i[i__ + j * e2i_dim1] = *betan * bb;
		}
	}
	decc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &ip2[1], ier);
	return 0;

	/* ----------------------------------------------------------- */

L13:
	/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
	i__1 = *nm1;
	for (j = 1; j <= i__1; ++j) {
		jm1 = j + *m1;
		i__4 = *nm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
			e2r[i__ + j * e2r_dim1] = -fjac[i__ + jm1 * fjac_dim1];
			e2i[i__ + j * e2i_dim1] = 0;
		}
		/* Computing MAX */
		i__4 = 1, i__2 = j - *mumas;
		/* Computing MIN */
		i__5 = *nm1, i__6 = j + *mlmas;
		i__3 = MIN(i__5, i__6);
		for (i__ = MAX(i__4, i__2); i__ <= i__3; ++i__) {
			ffma = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
			e2r[i__ + j * e2r_dim1] += *alphn * ffma;
			e2i[i__ + j * e2i_dim1] += *betan * ffma;
		}
	}
	goto L45;

	/* ----------------------------------------------------------- */

L4:
	/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		i__3 = linal_1.mbjac;
		for (i__ = 1; i__ <= i__3; ++i__) {
			imle = i__ + linal_1.mle;
			e2r[imle + j * e2r_dim1] = -fjac[i__ + j * fjac_dim1];
			e2i[imle + j * e2i_dim1] = 0;
		}
		/* Computing MAX */
		i__3 = 1, i__4 = *mumas + 2 - j;
		/* Computing MIN */
		i__5 = linal_1.mbb, i__6 = *mumas + 1 - j + *n;
		i__2 = MIN(i__5, i__6);
		for (i__ = MAX(i__3, i__4); i__ <= i__2; ++i__) {
			ib = i__ + linal_1.mdiff;
			bb = fmas[i__ + j * fmas_dim1];
			e2r[ib + j * e2r_dim1] += *alphn * bb;
			e2i[ib + j * e2i_dim1] = *betan * bb;
		}
	}
	decbc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
		linal_1.mue, &ip2[1], ier);
	return 0;

	/* ----------------------------------------------------------- */

L14:
	/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER */
	i__1 = *nm1;
	for (j = 1; j <= i__1; ++j) {
		jm1 = j + *m1;
		i__2 = linal_1.mbjac;
		for (i__ = 1; i__ <= i__2; ++i__) {
			e2r[i__ + linal_1.mle + j * e2r_dim1] = -fjac[i__ + jm1 *
				fjac_dim1];
			e2i[i__ + linal_1.mle + j * e2i_dim1] = 0;
		}
		i__2 = linal_1.mbb;
		for (i__ = 1; i__ <= i__2; ++i__) {
			ib = i__ + linal_1.mdiff;
			ffma = fmas[i__ + j * fmas_dim1];
			e2r[ib + j * e2r_dim1] += *alphn * ffma;
			e2i[ib + j * e2i_dim1] += *betan * ffma;
		}
	}
	goto L46;

	/* ----------------------------------------------------------- */

L5:
	/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
			bb = fmas[i__ + j * fmas_dim1];
			e2r[i__ + j * e2r_dim1] = bb * *alphn - fjac[i__ + j * fjac_dim1];
			e2i[i__ + j * e2i_dim1] = bb * *betan;
		}
	}
	decc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &ip2[1], ier);
	return 0;

	/* ----------------------------------------------------------- */

L15:
	/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
	i__1 = *nm1;
	for (j = 1; j <= i__1; ++j) {
		jm1 = j + *m1;
		i__2 = *nm1;
		for (i__ = 1; i__ <= i__2; ++i__) {
			e2r[i__ + j * e2r_dim1] = *alphn * fmas[i__ + j * fmas_dim1] -
				fjac[i__ + jm1 * fjac_dim1];
			e2i[i__ + j * e2i_dim1] = *betan * fmas[i__ + j * fmas_dim1];
		}
	}
	goto L45;

	/* ----------------------------------------------------------- */

L6:
	/* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
	/* ---  THIS OPTION IS NOT PROVIDED */
	return 0;

	/* ----------------------------------------------------------- */

L7:
	/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
	i__1 = *n - 1;
	for (j = 1; j <= i__1; ++j) {
		j1 = j + 1;
		e2r[j1 + j * e2r_dim1] = -fjac[j1 + j * fjac_dim1];
		e2i[j1 + j * e2i_dim1] = 0;
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) {
			e2i[i__ + j * e2i_dim1] = 0;
			e2r[i__ + j * e2r_dim1] = -fjac[i__ + j * fjac_dim1];
		}
		e2r[j + j * e2r_dim1] += *alphn;
		e2i[j + j * e2i_dim1] = *betan;
	}
	dechc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &c__1, &ip2[1], ier);
	return 0;

	/* ----------------------------------------------------------- */

L55:
	return 0;
} /* decomc_ */


/*     END OF SUBROUTINE DECOMC */

/* *********************************************************** */

/* Subroutine */ int slvrar_(typeLInt *n, typeRNum *fjac, typeLInt *ldjac,
	typeLInt *mljac, typeLInt *mujac, typeRNum *fmas, typeLInt *ldmas,
	typeLInt *mlmas, typeLInt *mumas, typeLInt *m1, typeLInt *m2, typeLInt *
	nm1, typeRNum *fac1, typeRNum *e1, typeLInt *lde1, typeRNum *z1,
	typeRNum *f1, typeLInt *ip1, typeLInt *iphes, typeLInt *ier, typeLInt *
	ijob)
{
	/* System generated locals */
	typeLInt fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1,
		e1_offset, i__1, i__2, i__3, i__4, i__5, i__6;

	/* Local variables */
	static typeLInt i__, j, k;
	static typeRNum s1;
	static typeLInt mm, mp, im1, mp1, jkm;
	extern /* Subroutine */ int sol_(typeLInt *, typeLInt *, typeRNum *,
		typeRNum *, typeLInt *);
	static typeRNum sum1;
	extern /* Subroutine */ int solb_(typeLInt *, typeLInt *, typeRNum *,
		typeLInt *, typeLInt *, typeRNum *, typeLInt *), solh_(typeLInt *,
			typeLInt *, typeRNum *, typeLInt *, typeRNum *, typeLInt *);
	static typeRNum zsafe;


	/* Parameter adjustments */
	--iphes;
	--f1;
	--z1;
	fjac_dim1 = *ldjac;
	fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;
	--ip1;
	fmas_dim1 = *ldmas;
	fmas_offset = 1 + fmas_dim1;
	fmas -= fmas_offset;
	e1_dim1 = *lde1;
	e1_offset = 1 + e1_dim1;
	e1 -= e1_offset;

	/* Function Body */
	switch (*ijob) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	case 4:  goto L4;
	case 5:  goto L5;
	case 6:  goto L6;
	case 7:  goto L7;
	case 8:  goto L55;
	case 9:  goto L55;
	case 10:  goto L55;
	case 11:  goto L11;
	case 12:  goto L12;
	case 13:  goto L13;
	case 14:  goto L13;
	case 15:  goto L15;
	}

	/* ----------------------------------------------------------- */

L1:
	/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		z1[i__] -= f1[i__] * *fac1;
	}
	sol_(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
	return 0;

	/* ----------------------------------------------------------- */

L11:
	/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		z1[i__] -= f1[i__] * *fac1;
	}
L48:
	mm = *m1 / *m2;
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
		sum1 = 0;
		for (k = mm - 1; k >= 0; --k) {
			jkm = j + k * *m2;
			sum1 = (z1[jkm] + sum1) / *fac1;
			i__2 = *nm1;
			for (i__ = 1; i__ <= i__2; ++i__) {
				im1 = i__ + *m1;
				z1[im1] += fjac[i__ + jkm * fjac_dim1] * sum1;
			}
		}
	}
	sol_(nm1, lde1, &e1[e1_offset], &z1[*m1 + 1], &ip1[1]);
L49:
	for (i__ = *m1; i__ >= 1; --i__) {
		z1[i__] = (z1[i__] + z1[*m2 + i__]) / *fac1;
	}
	return 0;

	/* ----------------------------------------------------------- */

L2:
	/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		z1[i__] -= f1[i__] * *fac1;
	}
	solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[1], &ip1[1]
	);
	return 0;

	/* ----------------------------------------------------------- */

L12:
	/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		z1[i__] -= f1[i__] * *fac1;
	}
L45:
	mm = *m1 / *m2;
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
		sum1 = 0;
		for (k = mm - 1; k >= 0; --k) {
			jkm = j + k * *m2;
			sum1 = (z1[jkm] + sum1) / *fac1;
			/* Computing MAX */
			i__2 = 1, i__3 = j - *mujac;
			/* Computing MIN */
			i__5 = *nm1, i__6 = j + *mljac;
			i__4 = MIN(i__5, i__6);
			for (i__ = MAX(i__2, i__3); i__ <= i__4; ++i__) {
				im1 = i__ + *m1;
				z1[im1] += fjac[i__ + *mujac + 1 - j + jkm * fjac_dim1] *
					sum1;
			}
		}
	}
	solb_(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[*m1 + 1],
		&ip1[1]);
	goto L49;

	/* ----------------------------------------------------------- */

L3:
	/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s1 = 0;
		/* Computing MAX */
		i__4 = 1, i__2 = i__ - *mlmas;
		/* Computing MIN */
		i__5 = *n, i__6 = i__ + *mumas;
		i__3 = MIN(i__5, i__6);
		for (j = MAX(i__4, i__2); j <= i__3; ++j) {
			s1 -= fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j];
		}
		z1[i__] += s1 * *fac1;
	}
	sol_(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
	return 0;

	/* ----------------------------------------------------------- */

L13:
	/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
	i__1 = *m1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		z1[i__] -= f1[i__] * *fac1;
	}
	i__1 = *nm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		im1 = i__ + *m1;
		s1 = 0;
		/* Computing MAX */
		i__3 = 1, i__4 = i__ - *mlmas;
		/* Computing MIN */
		i__5 = *nm1, i__6 = i__ + *mumas;
		i__2 = MIN(i__5, i__6);
		for (j = MAX(i__3, i__4); j <= i__2; ++j) {
			s1 -= fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j + *m1]
				;
		}
		z1[im1] += s1 * *fac1;
	}
	if (*ijob == 14) {
		goto L45;
	}
	goto L48;

	/* ----------------------------------------------------------- */

L4:
	/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s1 = 0;
		/* Computing MAX */
		i__2 = 1, i__3 = i__ - *mlmas;
		/* Computing MIN */
		i__5 = *n, i__6 = i__ + *mumas;
		i__4 = MIN(i__5, i__6);
		for (j = MAX(i__2, i__3); j <= i__4; ++j) {
			s1 -= fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j];
		}
		z1[i__] += s1 * *fac1;
	}
	solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[1], &ip1[1]
	);
	return 0;

	/* ----------------------------------------------------------- */

L5:
	/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s1 = 0;
		i__4 = *n;
		for (j = 1; j <= i__4; ++j) {
			s1 -= fmas[i__ + j * fmas_dim1] * f1[j];
		}
		z1[i__] += s1 * *fac1;
	}
	sol_(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
	return 0;

	/* ----------------------------------------------------------- */

L15:
	/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
	i__1 = *m1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		z1[i__] -= f1[i__] * *fac1;
	}
	i__1 = *nm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		im1 = i__ + *m1;
		s1 = 0;
		i__4 = *nm1;
		for (j = 1; j <= i__4; ++j) {
			s1 -= fmas[i__ + j * fmas_dim1] * f1[j + *m1];
		}
		z1[im1] += s1 * *fac1;
	}
	goto L48;

	/* ----------------------------------------------------------- */

L6:
	/* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
	/* ---  THIS OPTION IS NOT PROVIDED */
	return 0;

	/* ----------------------------------------------------------- */

L7:
	/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		z1[i__] -= f1[i__] * *fac1;
	}
	for (mm = *n - 2; mm >= 1; --mm) {
		mp = *n - mm;
		mp1 = mp - 1;
		i__ = iphes[mp];
		if (i__ == mp) {
			goto L746;
		}
		zsafe = z1[mp];
		z1[mp] = z1[i__];
		z1[i__] = zsafe;
	L746:
		i__1 = *n;
		for (i__ = mp + 1; i__ <= i__1; ++i__) {
			z1[i__] -= fjac[i__ + mp1 * fjac_dim1] * z1[mp];
		}
	}
	solh_(n, lde1, &e1[e1_offset], &c__1, &z1[1], &ip1[1]);
	i__1 = *n - 2;
	for (mm = 1; mm <= i__1; ++mm) {
		mp = *n - mm;
		mp1 = mp - 1;
		i__4 = *n;
		for (i__ = mp + 1; i__ <= i__4; ++i__) {
			z1[i__] += fjac[i__ + mp1 * fjac_dim1] * z1[mp];
		}
		i__ = iphes[mp];
		if (i__ == mp) {
			goto L750;
		}
		zsafe = z1[mp];
		z1[mp] = z1[i__];
		z1[i__] = zsafe;
	L750:
		;
	}
	return 0;

	/* ----------------------------------------------------------- */

L55:
	return 0;
} /* slvrar_ */


/*     END OF SUBROUTINE SLVRAR */

/* *********************************************************** */

/* Subroutine */ int slvrai_(typeLInt *n, typeRNum *fjac, typeLInt *ldjac,
	typeLInt *mljac, typeLInt *mujac, typeRNum *fmas, typeLInt *ldmas,
	typeLInt *mlmas, typeLInt *mumas, typeLInt *m1, typeLInt *m2, typeLInt *
	nm1, typeRNum *alphn, typeRNum *betan, typeRNum *e2r,
	typeRNum *e2i, typeLInt *lde1, typeRNum *z2, typeRNum *z3,
	typeRNum *f2, typeRNum *f3, typeRNum *cont, typeLInt *ip2,
	typeLInt *iphes, typeLInt *ier, typeLInt *ijob)
{
	/* System generated locals */
	typeLInt fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e2r_dim1,
		e2r_offset, e2i_dim1, e2i_offset, i__1, i__2, i__3, i__4, i__5,
		i__6;
	typeRNum d__1, d__2;

	/* Local variables */
	static typeLInt i__, j, k;
	static typeRNum s2, s3, bb;
	static typeLInt mm, mp, im1, jm1, mp1;
	static typeRNum z2i, z3i;
	static typeLInt jkm, mpi;
	static typeRNum sum2, sum3, abno;
	extern /* Subroutine */ int solc_(typeLInt *, typeLInt *, typeRNum *,
		typeRNum *, typeRNum *, typeRNum *, typeLInt *);
	static typeLInt iimu;
	static typeRNum sumh, e1imp;
	extern /* Subroutine */ int solbc_(typeLInt *, typeLInt *, typeRNum *,
		typeRNum *, typeLInt *, typeLInt *, typeRNum *, typeRNum *,
		typeLInt *);
	static typeRNum zsafe;
	extern /* Subroutine */ int solhc_(typeLInt *, typeLInt *, typeRNum *,
		typeRNum *, typeLInt *, typeRNum *, typeRNum *, typeLInt *);


	/* Parameter adjustments */
	--iphes;
	--f3;
	--f2;
	--z3;
	--z2;
	fjac_dim1 = *ldjac;
	fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;
	--ip2;
	fmas_dim1 = *ldmas;
	fmas_offset = 1 + fmas_dim1;
	fmas -= fmas_offset;
	e2i_dim1 = *lde1;
	e2i_offset = 1 + e2i_dim1;
	e2i -= e2i_offset;
	e2r_dim1 = *lde1;
	e2r_offset = 1 + e2r_dim1;
	e2r -= e2r_offset;

	/* Function Body */
	switch (*ijob) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	case 4:  goto L4;
	case 5:  goto L5;
	case 6:  goto L6;
	case 7:  goto L7;
	case 8:  goto L55;
	case 9:  goto L55;
	case 10:  goto L55;
	case 11:  goto L11;
	case 12:  goto L12;
	case 13:  goto L13;
	case 14:  goto L13;
	case 15:  goto L15;
	}

	/* ----------------------------------------------------------- */

L1:
	/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s2 = -f2[i__];
		s3 = -f3[i__];
		z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
		z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
	}
	solc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]
	);
	return 0;

	/* ----------------------------------------------------------- */

L11:
	/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s2 = -f2[i__];
		s3 = -f3[i__];
		z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
		z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
	}
L48:
	/* Computing 2nd power */
	d__1 = *alphn;
	/* Computing 2nd power */
	d__2 = *betan;
	abno = d__1 * d__1 + d__2 * d__2;
	mm = *m1 / *m2;
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
		sum2 = 0;
		sum3 = 0;
		for (k = mm - 1; k >= 0; --k) {
			jkm = j + k * *m2;
			sumh = (z2[jkm] + sum2) / abno;
			sum3 = (z3[jkm] + sum3) / abno;
			sum2 = sumh * *alphn + sum3 * *betan;
			sum3 = sum3 * *alphn - sumh * *betan;
			i__2 = *nm1;
			for (i__ = 1; i__ <= i__2; ++i__) {
				im1 = i__ + *m1;
				z2[im1] += fjac[i__ + jkm * fjac_dim1] * sum2;
				z3[im1] += fjac[i__ + jkm * fjac_dim1] * sum3;
			}
		}
	}
	solc_(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[*m1 + 1], &z3[*
		m1 + 1], &ip2[1]);
L49:
	for (i__ = *m1; i__ >= 1; --i__) {
		mpi = *m2 + i__;
		z2i = z2[i__] + z2[mpi];
		z3i = z3[i__] + z3[mpi];
		z3[i__] = (z3i * *alphn - z2i * *betan) / abno;
		z2[i__] = (z2i * *alphn + z3i * *betan) / abno;
	}
	return 0;

	/* ----------------------------------------------------------- */

L2:
	/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s2 = -f2[i__];
		s3 = -f3[i__];
		z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
		z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
	}
	solbc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
		linal_1.mue, &z2[1], &z3[1], &ip2[1]);
	return 0;

	/* ----------------------------------------------------------- */

L12:
	/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s2 = -f2[i__];
		s3 = -f3[i__];
		z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
		z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
	}
L45:
	/* Computing 2nd power */
	d__1 = *alphn;
	/* Computing 2nd power */
	d__2 = *betan;
	abno = d__1 * d__1 + d__2 * d__2;
	mm = *m1 / *m2;
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
		sum2 = 0;
		sum3 = 0;
		for (k = mm - 1; k >= 0; --k) {
			jkm = j + k * *m2;
			sumh = (z2[jkm] + sum2) / abno;
			sum3 = (z3[jkm] + sum3) / abno;
			sum2 = sumh * *alphn + sum3 * *betan;
			sum3 = sum3 * *alphn - sumh * *betan;
			/* Computing MAX */
			i__2 = 1, i__3 = j - *mujac;
			/* Computing MIN */
			i__5 = *nm1, i__6 = j + *mljac;
			i__4 = MIN(i__5, i__6);
			for (i__ = MAX(i__2, i__3); i__ <= i__4; ++i__) {
				im1 = i__ + *m1;
				iimu = i__ + *mujac + 1 - j;
				z2[im1] += fjac[iimu + jkm * fjac_dim1] * sum2;
				z3[im1] += fjac[iimu + jkm * fjac_dim1] * sum3;
			}
		}
	}
	solbc_(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
		linal_1.mue, &z2[*m1 + 1], &z3[*m1 + 1], &ip2[1]);
	goto L49;

	/* ----------------------------------------------------------- */

L3:
	/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s2 = 0;
		s3 = 0;
		/* Computing MAX */
		i__4 = 1, i__2 = i__ - *mlmas;
		/* Computing MIN */
		i__5 = *n, i__6 = i__ + *mumas;
		i__3 = MIN(i__5, i__6);
		for (j = MAX(i__4, i__2); j <= i__3; ++j) {
			bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
			s2 -= bb * f2[j];
			s3 -= bb * f3[j];
		}
		z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
		z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
	}
	solc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]
	);
	return 0;

	/* ----------------------------------------------------------- */

L13:
	/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
	i__1 = *m1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s2 = -f2[i__];
		s3 = -f3[i__];
		z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
		z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
	}
	i__1 = *nm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		im1 = i__ + *m1;
		s2 = 0;
		s3 = 0;
		/* Computing MAX */
		i__3 = 1, i__4 = i__ - *mlmas;
		/* Computing MIN */
		i__5 = *nm1, i__6 = i__ + *mumas;
		i__2 = MIN(i__5, i__6);
		for (j = MAX(i__3, i__4); j <= i__2; ++j) {
			jm1 = j + *m1;
			bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
			s2 -= bb * f2[jm1];
			s3 -= bb * f3[jm1];
		}
		z2[im1] = z2[im1] + s2 * *alphn - s3 * *betan;
		z3[im1] = z3[im1] + s3 * *alphn + s2 * *betan;
	}
	if (*ijob == 14) {
		goto L45;
	}
	goto L48;

	/* ----------------------------------------------------------- */

L4:
	/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s2 = 0;
		s3 = 0;
		/* Computing MAX */
		i__2 = 1, i__3 = i__ - *mlmas;
		/* Computing MIN */
		i__5 = *n, i__6 = i__ + *mumas;
		i__4 = MIN(i__5, i__6);
		for (j = MAX(i__2, i__3); j <= i__4; ++j) {
			bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
			s2 -= bb * f2[j];
			s3 -= bb * f3[j];
		}
		z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
		z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
	}
	solbc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
		linal_1.mue, &z2[1], &z3[1], &ip2[1]);
	return 0;

	/* ----------------------------------------------------------- */

L5:
	/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s2 = 0;
		s3 = 0;
		i__4 = *n;
		for (j = 1; j <= i__4; ++j) {
			bb = fmas[i__ + j * fmas_dim1];
			s2 -= bb * f2[j];
			s3 -= bb * f3[j];
		}
		z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
		z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
	}
	solc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]
	);
	return 0;

	/* ----------------------------------------------------------- */

L15:
	/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
	i__1 = *m1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s2 = -f2[i__];
		s3 = -f3[i__];
		z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
		z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
	}
	i__1 = *nm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		im1 = i__ + *m1;
		s2 = 0;
		s3 = 0;
		i__4 = *nm1;
		for (j = 1; j <= i__4; ++j) {
			jm1 = j + *m1;
			bb = fmas[i__ + j * fmas_dim1];
			s2 -= bb * f2[jm1];
			s3 -= bb * f3[jm1];
		}
		z2[im1] = z2[im1] + s2 * *alphn - s3 * *betan;
		z3[im1] = z3[im1] + s3 * *alphn + s2 * *betan;
	}
	goto L48;

	/* ----------------------------------------------------------- */

L6:
	/* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
	/* ---  THIS OPTION IS NOT PROVIDED */
	return 0;

	/* ----------------------------------------------------------- */

L7:
	/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s2 = -f2[i__];
		s3 = -f3[i__];
		z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
		z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
	}
	for (mm = *n - 2; mm >= 1; --mm) {
		mp = *n - mm;
		mp1 = mp - 1;
		i__ = iphes[mp];
		if (i__ == mp) {
			goto L746;
		}
		zsafe = z2[mp];
		z2[mp] = z2[i__];
		z2[i__] = zsafe;
		zsafe = z3[mp];
		z3[mp] = z3[i__];
		z3[i__] = zsafe;
	L746:
		i__1 = *n;
		for (i__ = mp + 1; i__ <= i__1; ++i__) {
			e1imp = fjac[i__ + mp1 * fjac_dim1];
			z2[i__] -= e1imp * z2[mp];
			z3[i__] -= e1imp * z3[mp];
		}
	}
	solhc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &c__1, &z2[1], &z3[1],
		&ip2[1]);
	i__1 = *n - 2;
	for (mm = 1; mm <= i__1; ++mm) {
		mp = *n - mm;
		mp1 = mp - 1;
		i__4 = *n;
		for (i__ = mp + 1; i__ <= i__4; ++i__) {
			e1imp = fjac[i__ + mp1 * fjac_dim1];
			z2[i__] += e1imp * z2[mp];
			z3[i__] += e1imp * z3[mp];
		}
		i__ = iphes[mp];
		if (i__ == mp) {
			goto L750;
		}
		zsafe = z2[mp];
		z2[mp] = z2[i__];
		z2[i__] = zsafe;
		zsafe = z3[mp];
		z3[mp] = z3[i__];
		z3[i__] = zsafe;
	L750:
		;
	}
	return 0;

	/* ----------------------------------------------------------- */

L55:
	return 0;
} /* slvrai_ */


/*     END OF SUBROUTINE SLVRAI */

/* *********************************************************** */

/* Subroutine */ int slvrad_(typeLInt *n, typeRNum *fjac, typeLInt *ldjac,
	typeLInt *mljac, typeLInt *mujac, typeRNum *fmas, typeLInt *ldmas,
	typeLInt *mlmas, typeLInt *mumas, typeLInt *m1, typeLInt *m2, typeLInt *
	nm1, typeRNum *fac1, typeRNum *alphn, typeRNum *betan,
	typeRNum *e1, typeRNum *e2r, typeRNum *e2i, typeLInt *lde1,
	typeRNum *z1, typeRNum *z2, typeRNum *z3, typeRNum *f1,
	typeRNum *f2, typeRNum *f3, typeRNum *cont, typeLInt *ip1,
	typeLInt *ip2, typeLInt *iphes, typeLInt *ier, typeLInt *ijob)
{
	/* System generated locals */
	typeLInt fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1,
		e1_offset, e2r_dim1, e2r_offset, e2i_dim1, e2i_offset, i__1, i__2,
		i__3, i__4, i__5, i__6;
	typeRNum d__1, d__2;

	/* Local variables */
	static typeLInt i__, j, k;
	static typeRNum s1, s2, s3, bb;
	static typeLInt mm, mp, j1b, j2b, im1, jm1, mp1;
	static typeRNum z2i, z3i;
	static typeLInt jkm, mpi;
	extern /* Subroutine */ int sol_(typeLInt *, typeLInt *, typeRNum *,
		typeRNum *, typeLInt *);
	static typeRNum sum1, sum2, sum3, ffja, abno;
	extern /* Subroutine */ int solb_(typeLInt *, typeLInt *, typeRNum *,
		typeLInt *, typeLInt *, typeRNum *, typeLInt *), solc_(typeLInt *,
			typeLInt *, typeRNum *, typeRNum *, typeRNum *, typeRNum *,
			typeLInt *), solh_(typeLInt *, typeLInt *, typeRNum *, typeLInt *,
				typeRNum *, typeLInt *);
	static typeRNum sumh, e1imp;
	extern /* Subroutine */ int solbc_(typeLInt *, typeLInt *, typeRNum *,
		typeRNum *, typeLInt *, typeLInt *, typeRNum *, typeRNum *,
		typeLInt *);
	static typeRNum zsafe;
	extern /* Subroutine */ int solhc_(typeLInt *, typeLInt *, typeRNum *,
		typeRNum *, typeLInt *, typeRNum *, typeRNum *, typeLInt *);


	/* Parameter adjustments */
	--iphes;
	--f3;
	--f2;
	--f1;
	--z3;
	--z2;
	--z1;
	fjac_dim1 = *ldjac;
	fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;
	--ip2;
	--ip1;
	fmas_dim1 = *ldmas;
	fmas_offset = 1 + fmas_dim1;
	fmas -= fmas_offset;
	e2i_dim1 = *lde1;
	e2i_offset = 1 + e2i_dim1;
	e2i -= e2i_offset;
	e2r_dim1 = *lde1;
	e2r_offset = 1 + e2r_dim1;
	e2r -= e2r_offset;
	e1_dim1 = *lde1;
	e1_offset = 1 + e1_dim1;
	e1 -= e1_offset;

	/* Function Body */
	switch (*ijob) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	case 4:  goto L4;
	case 5:  goto L5;
	case 6:  goto L6;
	case 7:  goto L7;
	case 8:  goto L55;
	case 9:  goto L55;
	case 10:  goto L55;
	case 11:  goto L11;
	case 12:  goto L12;
	case 13:  goto L13;
	case 14:  goto L13;
	case 15:  goto L15;
	}

	/* ----------------------------------------------------------- */

L1:
	/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s2 = -f2[i__];
		s3 = -f3[i__];
		z1[i__] -= f1[i__] * *fac1;
		z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
		z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
	}
	sol_(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
	solc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]
	);
	return 0;

	/* ----------------------------------------------------------- */

L11:
	/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s2 = -f2[i__];
		s3 = -f3[i__];
		z1[i__] -= f1[i__] * *fac1;
		z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
		z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
	}
L48:
	/* Computing 2nd power */
	d__1 = *alphn;
	/* Computing 2nd power */
	d__2 = *betan;
	abno = d__1 * d__1 + d__2 * d__2;
	mm = *m1 / *m2;
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
		sum1 = 0;
		sum2 = 0;
		sum3 = 0;
		for (k = mm - 1; k >= 0; --k) {
			jkm = j + k * *m2;
			sum1 = (z1[jkm] + sum1) / *fac1;
			sumh = (z2[jkm] + sum2) / abno;
			sum3 = (z3[jkm] + sum3) / abno;
			sum2 = sumh * *alphn + sum3 * *betan;
			sum3 = sum3 * *alphn - sumh * *betan;
			i__2 = *nm1;
			for (i__ = 1; i__ <= i__2; ++i__) {
				im1 = i__ + *m1;
				z1[im1] += fjac[i__ + jkm * fjac_dim1] * sum1;
				z2[im1] += fjac[i__ + jkm * fjac_dim1] * sum2;
				z3[im1] += fjac[i__ + jkm * fjac_dim1] * sum3;
			}
		}
	}
	sol_(nm1, lde1, &e1[e1_offset], &z1[*m1 + 1], &ip1[1]);
	solc_(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[*m1 + 1], &z3[*
		m1 + 1], &ip2[1]);
L49:
	for (i__ = *m1; i__ >= 1; --i__) {
		mpi = *m2 + i__;
		z1[i__] = (z1[i__] + z1[mpi]) / *fac1;
		z2i = z2[i__] + z2[mpi];
		z3i = z3[i__] + z3[mpi];
		z3[i__] = (z3i * *alphn - z2i * *betan) / abno;
		z2[i__] = (z2i * *alphn + z3i * *betan) / abno;
	}
	return 0;

	/* ----------------------------------------------------------- */

L2:
	/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s2 = -f2[i__];
		s3 = -f3[i__];
		z1[i__] -= f1[i__] * *fac1;
		z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
		z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
	}
	solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[1], &ip1[1]
	);
	solbc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
		linal_1.mue, &z2[1], &z3[1], &ip2[1]);
	return 0;

	/* ----------------------------------------------------------- */

L12:
	/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s2 = -f2[i__];
		s3 = -f3[i__];
		z1[i__] -= f1[i__] * *fac1;
		z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
		z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
	}
L45:
	/* Computing 2nd power */
	d__1 = *alphn;
	/* Computing 2nd power */
	d__2 = *betan;
	abno = d__1 * d__1 + d__2 * d__2;
	mm = *m1 / *m2;
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
		sum1 = 0;
		sum2 = 0;
		sum3 = 0;
		for (k = mm - 1; k >= 0; --k) {
			jkm = j + k * *m2;
			sum1 = (z1[jkm] + sum1) / *fac1;
			sumh = (z2[jkm] + sum2) / abno;
			sum3 = (z3[jkm] + sum3) / abno;
			sum2 = sumh * *alphn + sum3 * *betan;
			sum3 = sum3 * *alphn - sumh * *betan;
			/* Computing MAX */
			i__2 = 1, i__3 = j - *mujac;
			/* Computing MIN */
			i__5 = *nm1, i__6 = j + *mljac;
			i__4 = MIN(i__5, i__6);
			for (i__ = MAX(i__2, i__3); i__ <= i__4; ++i__) {
				im1 = i__ + *m1;
				ffja = fjac[i__ + *mujac + 1 - j + jkm * fjac_dim1];
				z1[im1] += ffja * sum1;
				z2[im1] += ffja * sum2;
				z3[im1] += ffja * sum3;
			}
		}
	}
	solb_(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[*m1 + 1],
		&ip1[1]);
	solbc_(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
		linal_1.mue, &z2[*m1 + 1], &z3[*m1 + 1], &ip2[1]);
	goto L49;

	/* ----------------------------------------------------------- */

L3:
	/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s1 = 0;
		s2 = 0;
		s3 = 0;
		/* Computing MAX */
		i__4 = 1, i__2 = i__ - *mlmas;
		/* Computing MIN */
		i__5 = *n, i__6 = i__ + *mumas;
		i__3 = MIN(i__5, i__6);
		for (j = MAX(i__4, i__2); j <= i__3; ++j) {
			bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
			s1 -= bb * f1[j];
			s2 -= bb * f2[j];
			s3 -= bb * f3[j];
		}
		z1[i__] += s1 * *fac1;
		z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
		z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
	}
	sol_(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
	solc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]
	);
	return 0;

	/* ----------------------------------------------------------- */

L13:
	/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
	i__1 = *m1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s2 = -f2[i__];
		s3 = -f3[i__];
		z1[i__] -= f1[i__] * *fac1;
		z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
		z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
	}
	i__1 = *nm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		im1 = i__ + *m1;
		s1 = 0;
		s2 = 0;
		s3 = 0;
		/* Computing MAX */
		i__3 = 1, i__4 = i__ - *mlmas;
		j1b = MAX(i__3, i__4);
		/* Computing MIN */
		i__3 = *nm1, i__4 = i__ + *mumas;
		j2b = MIN(i__3, i__4);
		i__3 = j2b;
		for (j = j1b; j <= i__3; ++j) {
			jm1 = j + *m1;
			bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
			s1 -= bb * f1[jm1];
			s2 -= bb * f2[jm1];
			s3 -= bb * f3[jm1];
		}
		z1[im1] += s1 * *fac1;
		z2[im1] = z2[im1] + s2 * *alphn - s3 * *betan;
		z3[im1] = z3[im1] + s3 * *alphn + s2 * *betan;
	}
	if (*ijob == 14) {
		goto L45;
	}
	goto L48;

	/* ----------------------------------------------------------- */

L4:
	/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s1 = 0;
		s2 = 0;
		s3 = 0;
		/* Computing MAX */
		i__3 = 1, i__4 = i__ - *mlmas;
		/* Computing MIN */
		i__5 = *n, i__6 = i__ + *mumas;
		i__2 = MIN(i__5, i__6);
		for (j = MAX(i__3, i__4); j <= i__2; ++j) {
			bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
			s1 -= bb * f1[j];
			s2 -= bb * f2[j];
			s3 -= bb * f3[j];
		}
		z1[i__] += s1 * *fac1;
		z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
		z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
	}
	solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[1], &ip1[1]
	);
	solbc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
		linal_1.mue, &z2[1], &z3[1], &ip2[1]);
	return 0;

	/* ----------------------------------------------------------- */

L5:
	/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s1 = 0;
		s2 = 0;
		s3 = 0;
		i__2 = *n;
		for (j = 1; j <= i__2; ++j) {
			bb = fmas[i__ + j * fmas_dim1];
			s1 -= bb * f1[j];
			s2 -= bb * f2[j];
			s3 -= bb * f3[j];
		}
		z1[i__] += s1 * *fac1;
		z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
		z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
	}
	sol_(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
	solc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]
	);
	return 0;

	/* ----------------------------------------------------------- */

L15:
	/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
	i__1 = *m1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s2 = -f2[i__];
		s3 = -f3[i__];
		z1[i__] -= f1[i__] * *fac1;
		z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
		z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
	}
	i__1 = *nm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		im1 = i__ + *m1;
		s1 = 0;
		s2 = 0;
		s3 = 0;
		i__2 = *nm1;
		for (j = 1; j <= i__2; ++j) {
			jm1 = j + *m1;
			bb = fmas[i__ + j * fmas_dim1];
			s1 -= bb * f1[jm1];
			s2 -= bb * f2[jm1];
			s3 -= bb * f3[jm1];
		}
		z1[im1] += s1 * *fac1;
		z2[im1] = z2[im1] + s2 * *alphn - s3 * *betan;
		z3[im1] = z3[im1] + s3 * *alphn + s2 * *betan;
	}
	goto L48;

	/* ----------------------------------------------------------- */

L6:
	/* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
	/* ---  THIS OPTION IS NOT PROVIDED */
	return 0;

	/* ----------------------------------------------------------- */

L7:
	/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		s2 = -f2[i__];
		s3 = -f3[i__];
		z1[i__] -= f1[i__] * *fac1;
		z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
		z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
	}
	for (mm = *n - 2; mm >= 1; --mm) {
		mp = *n - mm;
		mp1 = mp - 1;
		i__ = iphes[mp];
		if (i__ == mp) {
			goto L746;
		}
		zsafe = z1[mp];
		z1[mp] = z1[i__];
		z1[i__] = zsafe;
		zsafe = z2[mp];
		z2[mp] = z2[i__];
		z2[i__] = zsafe;
		zsafe = z3[mp];
		z3[mp] = z3[i__];
		z3[i__] = zsafe;
	L746:
		i__1 = *n;
		for (i__ = mp + 1; i__ <= i__1; ++i__) {
			e1imp = fjac[i__ + mp1 * fjac_dim1];
			z1[i__] -= e1imp * z1[mp];
			z2[i__] -= e1imp * z2[mp];
			z3[i__] -= e1imp * z3[mp];
		}
	}
	solh_(n, lde1, &e1[e1_offset], &c__1, &z1[1], &ip1[1]);
	solhc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &c__1, &z2[1], &z3[1],
		&ip2[1]);
	i__1 = *n - 2;
	for (mm = 1; mm <= i__1; ++mm) {
		mp = *n - mm;
		mp1 = mp - 1;
		i__2 = *n;
		for (i__ = mp + 1; i__ <= i__2; ++i__) {
			e1imp = fjac[i__ + mp1 * fjac_dim1];
			z1[i__] += e1imp * z1[mp];
			z2[i__] += e1imp * z2[mp];
			z3[i__] += e1imp * z3[mp];
		}
		i__ = iphes[mp];
		if (i__ == mp) {
			goto L750;
		}
		zsafe = z1[mp];
		z1[mp] = z1[i__];
		z1[i__] = zsafe;
		zsafe = z2[mp];
		z2[mp] = z2[i__];
		z2[i__] = zsafe;
		zsafe = z3[mp];
		z3[mp] = z3[i__];
		z3[i__] = zsafe;
	L750:
		;
	}
	return 0;

	/* ----------------------------------------------------------- */

L55:
	return 0;
} /* slvrad_ */


/*     END OF SUBROUTINE SLVRAD */

/* *********************************************************** */

/* Subroutine */ int estrad_(typeLInt *n, typeRNum *fjac, typeLInt *ldjac,
	typeLInt *mljac, typeLInt *mujac, typeRNum *fmas, typeLInt *ldmas,
	typeLInt *mlmas, typeLInt *mumas, typeRNum *h__, typeRNum *dd1,
	typeRNum *dd2, typeRNum *dd3, S_fp fcn, typeLInt *nfcn, typeRNum
	*y0, typeRNum *y, typeLInt *ijob, typeRNum *x, typeLInt *m1,
	typeLInt *m2, typeLInt *nm1, typeRNum *e1, typeLInt *lde1, typeRNum *
	z1, typeRNum *z2, typeRNum *z3, typeRNum *cont, typeRNum *f1,
	typeRNum *f2, typeLInt *ip1, typeLInt *iphes, typeRNum *scal,
	typeRNum *err, typeLogical *first, typeLogical *reject, typeRNum *fac1,
	ctypeRNum *tvec, ctypeRNum *xvec, ctypeRNum *uvec, ctypeRNum *pvec, ctypeRNum *dcdxvec,
	typeGRAMPC *grampc, typeffctPtr pfct)
{
	/* System generated locals */
	typeLInt fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1,
		e1_offset, i__1, i__2, i__3, i__4, i__5, i__6;
	typeRNum d__1;

	/* Local variables */
	static typeLInt i__, j, k, mm, mp, im1;
	extern /* Subroutine */ int sol_(typeLInt *, typeLInt *, typeRNum *,
		typeRNum *, typeLInt *);
	static typeRNum sum, hee1, hee2, hee3, sum1;
	extern /* Subroutine */ int solb_(typeLInt *, typeLInt *, typeRNum *,
		typeLInt *, typeLInt *, typeRNum *, typeLInt *), solh_(typeLInt *,
			typeLInt *, typeRNum *, typeLInt *, typeRNum *, typeLInt *);
	static typeRNum zsafe;

	/* Parameter adjustments */
	--scal;
	--iphes;
	--f2;
	--f1;
	--cont;
	--z3;
	--z2;
	--z1;
	--y;
	--y0;
	fjac_dim1 = *ldjac;
	fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;
	--ip1;
	fmas_dim1 = *ldmas;
	fmas_offset = 1 + fmas_dim1;
	fmas -= fmas_offset;
	e1_dim1 = *lde1;
	e1_offset = 1 + e1_dim1;
	e1 -= e1_offset;

	/* Function Body */
	hee1 = *dd1 / *h__;
	hee2 = *dd2 / *h__;
	hee3 = *dd3 / *h__;
	switch (*ijob) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	case 4:  goto L4;
	case 5:  goto L5;
	case 6:  goto L6;
	case 7:  goto L7;
	case 8:  goto L55;
	case 9:  goto L55;
	case 10:  goto L55;
	case 11:  goto L11;
	case 12:  goto L12;
	case 13:  goto L13;
	case 14:  goto L14;
	case 15:  goto L15;
	}

L1:
	/* ------  B=IDENTITY, JACOBIAN A FULL MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
		cont[i__] = f2[i__] + y0[i__];
	}
	sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
	goto L77;

L11:
	/* ------  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
		cont[i__] = f2[i__] + y0[i__];
	}
L48:
	mm = *m1 / *m2;
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
		sum1 = 0;
		for (k = mm - 1; k >= 0; --k) {
			sum1 = (cont[j + k * *m2] + sum1) / *fac1;
			i__2 = *nm1;
			for (i__ = 1; i__ <= i__2; ++i__) {
				im1 = i__ + *m1;
				cont[im1] += fjac[i__ + (j + k * *m2) * fjac_dim1] * sum1;
			}
		}
	}
	sol_(nm1, lde1, &e1[e1_offset], &cont[*m1 + 1], &ip1[1]);
	for (i__ = *m1; i__ >= 1; --i__) {
		cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
	}
	goto L77;

L2:
	/* ------  B=IDENTITY, JACOBIAN A BANDED MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
		cont[i__] = f2[i__] + y0[i__];
	}
	solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &ip1[
		1]);
	goto L77;

L12:
	/* ------  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
		cont[i__] = f2[i__] + y0[i__];
	}
L45:
	mm = *m1 / *m2;
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
		sum1 = 0;
		for (k = mm - 1; k >= 0; --k) {
			sum1 = (cont[j + k * *m2] + sum1) / *fac1;
			/* Computing MAX */
			i__2 = 1, i__3 = j - *mujac;
			/* Computing MIN */
			i__5 = *nm1, i__6 = j + *mljac;
			i__4 = MIN(i__5, i__6);
			for (i__ = MAX(i__2, i__3); i__ <= i__4; ++i__) {
				im1 = i__ + *m1;
				cont[im1] += fjac[i__ + *mujac + 1 - j + (j + k * *m2) *
					fjac_dim1] * sum1;
			}
		}
	}
	solb_(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[*m1 +
		1], &ip1[1]);
	for (i__ = *m1; i__ >= 1; --i__) {
		cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
	}
	goto L77;

L3:
	/* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		/* Computing MAX */
		i__4 = 1, i__2 = i__ - *mlmas;
		/* Computing MIN */
		i__5 = *n, i__6 = i__ + *mumas;
		i__3 = MIN(i__5, i__6);
		for (j = MAX(i__4, i__2); j <= i__3; ++j) {
			sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j];
		}
		f2[i__] = sum;
		cont[i__] = sum + y0[i__];
	}
	sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
	goto L77;

L13:
	/* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
	i__1 = *m1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
		cont[i__] = f2[i__] + y0[i__];
	}
	i__1 = *n;
	for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
		f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
	}
	i__1 = *nm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		/* Computing MAX */
		i__3 = 1, i__4 = i__ - *mlmas;
		/* Computing MIN */
		i__5 = *nm1, i__6 = i__ + *mumas;
		i__2 = MIN(i__5, i__6);
		for (j = MAX(i__3, i__4); j <= i__2; ++j) {
			sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j + *
				m1];
		}
		im1 = i__ + *m1;
		f2[im1] = sum;
		cont[im1] = sum + y0[im1];
	}
	goto L48;

L4:
	/* ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		/* Computing MAX */
		i__2 = 1, i__3 = i__ - *mlmas;
		/* Computing MIN */
		i__5 = *n, i__6 = i__ + *mumas;
		i__4 = MIN(i__5, i__6);
		for (j = MAX(i__2, i__3); j <= i__4; ++j) {
			sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j];
		}
		f2[i__] = sum;
		cont[i__] = sum + y0[i__];
	}
	solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &ip1[
		1]);
	goto L77;

L14:
	/* ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER */
	i__1 = *m1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
		cont[i__] = f2[i__] + y0[i__];
	}
	i__1 = *n;
	for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
		f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
	}
	i__1 = *nm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		/* Computing MAX */
		i__4 = 1, i__2 = i__ - *mlmas;
		/* Computing MIN */
		i__5 = *nm1, i__6 = i__ + *mumas;
		i__3 = MIN(i__5, i__6);
		for (j = MAX(i__4, i__2); j <= i__3; ++j) {
			sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j + *
				m1];
		}
		im1 = i__ + *m1;
		f2[im1] = sum;
		cont[im1] = sum + y0[im1];
	}
	goto L45;

L5:
	/* ------  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		i__3 = *n;
		for (j = 1; j <= i__3; ++j) {
			sum += fmas[i__ + j * fmas_dim1] * f1[j];
		}
		f2[i__] = sum;
		cont[i__] = sum + y0[i__];
	}
	sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
	goto L77;

L15:
	/* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
	i__1 = *m1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
		cont[i__] = f2[i__] + y0[i__];
	}
	i__1 = *n;
	for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
		f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
	}
	i__1 = *nm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		i__3 = *nm1;
		for (j = 1; j <= i__3; ++j) {
			sum += fmas[i__ + j * fmas_dim1] * f1[j + *m1];
		}
		im1 = i__ + *m1;
		f2[im1] = sum;
		cont[im1] = sum + y0[im1];
	}
	goto L48;

L6:
	/* ------  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
	/* ------  THIS OPTION IS NOT PROVIDED */
	return 0;

L7:
	/* ------  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
		cont[i__] = f2[i__] + y0[i__];
	}
	for (mm = *n - 2; mm >= 1; --mm) {
		mp = *n - mm;
		i__ = iphes[mp];
		if (i__ == mp) {
			goto L310;
		}
		zsafe = cont[mp];
		cont[mp] = cont[i__];
		cont[i__] = zsafe;
	L310:
		i__1 = *n;
		for (i__ = mp + 1; i__ <= i__1; ++i__) {
			cont[i__] -= fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
		}
	}
	solh_(n, lde1, &e1[e1_offset], &c__1, &cont[1], &ip1[1]);
	i__1 = *n - 2;
	for (mm = 1; mm <= i__1; ++mm) {
		mp = *n - mm;
		i__3 = *n;
		for (i__ = mp + 1; i__ <= i__3; ++i__) {
			cont[i__] += fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
		}
		i__ = iphes[mp];
		if (i__ == mp) {
			goto L440;
		}
		zsafe = cont[mp];
		cont[mp] = cont[i__];
		cont[i__] = zsafe;
	L440:
		;
	}

	/* -------------------------------------- */

L77:
	*err = 0;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		/* Computing 2nd power */
		d__1 = cont[i__] / scal[i__];
		*err += d__1 * d__1;
	}
	/* Computing MAX */
	d__1 = SQRT(*err / *n);
	*err = MAX(d__1, (typeRNum) 1e-10);

	if (*err < 1) {
		return 0;
	}
	if (*first || *reject) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			cont[i__] = y[i__] + cont[i__];
		}
		(*fcn)(n, x, &cont[1], &f1[1], tvec, xvec, uvec, pvec, dcdxvec, grampc, pfct);
		++(*nfcn);
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			cont[i__] = f1[i__] + f2[i__];
		}
		switch (*ijob) {
		case 1:  goto L31;
		case 2:  goto L32;
		case 3:  goto L31;
		case 4:  goto L32;
		case 5:  goto L31;
		case 6:  goto L32;
		case 7:  goto L33;
		case 8:  goto L55;
		case 9:  goto L55;
		case 10:  goto L55;
		case 11:  goto L41;
		case 12:  goto L42;
		case 13:  goto L41;
		case 14:  goto L42;
		case 15:  goto L41;
		}
		/* ------ FULL MATRIX OPTION */
	L31:
		sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
		goto L88;
		/* ------ FULL MATRIX OPTION, SECOND ORDER */
	L41:
		i__1 = *m2;
		for (j = 1; j <= i__1; ++j) {
			sum1 = 0;
			for (k = mm - 1; k >= 0; --k) {
				sum1 = (cont[j + k * *m2] + sum1) / *fac1;
				i__3 = *nm1;
				for (i__ = 1; i__ <= i__3; ++i__) {
					im1 = i__ + *m1;
					cont[im1] += fjac[i__ + (j + k * *m2) * fjac_dim1] * sum1;
				}
			}
		}
		sol_(nm1, lde1, &e1[e1_offset], &cont[*m1 + 1], &ip1[1]);
		for (i__ = *m1; i__ >= 1; --i__) {
			cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
		}
		goto L88;
		/* ------ BANDED MATRIX OPTION */
	L32:
		solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &
			ip1[1]);
		goto L88;
		/* ------ BANDED MATRIX OPTION, SECOND ORDER */
	L42:
		i__1 = *m2;
		for (j = 1; j <= i__1; ++j) {
			sum1 = 0;
			for (k = mm - 1; k >= 0; --k) {
				sum1 = (cont[j + k * *m2] + sum1) / *fac1;
				/* Computing MAX */
				i__3 = 1, i__4 = j - *mujac;
				/* Computing MIN */
				i__5 = *nm1, i__6 = j + *mljac;
				i__2 = MIN(i__5, i__6);
				for (i__ = MAX(i__3, i__4); i__ <= i__2; ++i__) {
					im1 = i__ + *m1;
					cont[im1] += fjac[i__ + *mujac + 1 - j + (j + k * *m2) *
						fjac_dim1] * sum1;
				}
			}
		}
		solb_(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[*
			m1 + 1], &ip1[1]);
		for (i__ = *m1; i__ >= 1; --i__) {
			cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
		}
		goto L88;
		/* ------ HESSENBERG MATRIX OPTION */
	L33:
		for (mm = *n - 2; mm >= 1; --mm) {
			mp = *n - mm;
			i__ = iphes[mp];
			if (i__ == mp) {
				goto L510;
			}
			zsafe = cont[mp];
			cont[mp] = cont[i__];
			cont[i__] = zsafe;
		L510:
			i__1 = *n;
			for (i__ = mp + 1; i__ <= i__1; ++i__) {
				cont[i__] -= fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
			}
		}
		solh_(n, lde1, &e1[e1_offset], &c__1, &cont[1], &ip1[1]);
		i__1 = *n - 2;
		for (mm = 1; mm <= i__1; ++mm) {
			mp = *n - mm;
			i__2 = *n;
			for (i__ = mp + 1; i__ <= i__2; ++i__) {
				cont[i__] += fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
			}
			i__ = iphes[mp];
			if (i__ == mp) {
				goto L640;
			}
			zsafe = cont[mp];
			cont[mp] = cont[i__];
			cont[i__] = zsafe;
		L640:
			;
		}
		/* ----------------------------------- */
	L88:
		*err = 0;
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			/* Computing 2nd power */
			d__1 = cont[i__] / scal[i__];
			*err += d__1 * d__1;
		}
		/* Computing MAX */
		d__1 = SQRT(*err / *n);
		*err = MAX(d__1, (typeRNum) 1e-10);
	}
	return 0;
	/* ----------------------------------------------------------- */
L55:
	return 0;
} /* estrad_ */


/*     END OF SUBROUTINE ESTRAD */

/* *********************************************************** */

/* Subroutine */ int estrav_(typeLInt *n, typeRNum *fjac, typeLInt *ldjac,
	typeLInt *mljac, typeLInt *mujac, typeRNum *fmas, typeLInt *ldmas,
	typeLInt *mlmas, typeLInt *mumas, typeRNum *h__, typeRNum *dd, S_fp
	fcn, typeLInt *nfcn, typeRNum *y0, typeRNum *y, typeLInt *ijob,
	typeRNum *x, typeLInt *m1, typeLInt *m2, typeLInt *nm1, typeLInt *ns,
	typeLInt *nns, typeRNum *e1, typeLInt *lde1, typeRNum *zz,
	typeRNum *cont, typeRNum *ff, typeLInt *ip1, typeLInt *iphes,
	typeRNum *scal, typeRNum *err, typeLogical *first, typeLogical *reject,
	typeRNum *fac1, ctypeRNum *tvec, ctypeRNum *xvec, ctypeRNum *uvec, ctypeRNum *pvec, ctypeRNum *dcdxvec, typeGRAMPC *grampc, typeffctPtr pfct)
{
	/* System generated locals */
	typeLInt fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1,
		e1_offset, i__1, i__2, i__3, i__4, i__5, i__6;
	typeRNum d__1;

	/* Local variables */
	static typeLInt i__, j, k, mm, mp, im1;
	extern /* Subroutine */ int sol_(typeLInt *, typeLInt *, typeRNum *,
		typeRNum *, typeLInt *);
	static typeRNum sum, sum1;
	extern /* Subroutine */ int solb_(typeLInt *, typeLInt *, typeRNum *,
		typeLInt *, typeLInt *, typeRNum *, typeLInt *), solh_(typeLInt *,
			typeLInt *, typeRNum *, typeLInt *, typeRNum *, typeLInt *);
	static typeRNum zsafe;

	/* Parameter adjustments */
	--scal;
	--iphes;
	--cont;
	--y;
	--y0;
	fjac_dim1 = *ldjac;
	fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;
	--ip1;
	fmas_dim1 = *ldmas;
	fmas_offset = 1 + fmas_dim1;
	fmas -= fmas_offset;
	--dd;
	--ff;
	--zz;
	e1_dim1 = *lde1;
	e1_offset = 1 + e1_dim1;
	e1 -= e1_offset;

	/* Function Body */
	switch (*ijob) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	case 4:  goto L4;
	case 5:  goto L5;
	case 6:  goto L6;
	case 7:  goto L7;
	case 8:  goto L55;
	case 9:  goto L55;
	case 10:  goto L55;
	case 11:  goto L11;
	case 12:  goto L12;
	case 13:  goto L13;
	case 14:  goto L14;
	case 15:  goto L15;
	}

L1:
	/* ------  B=IDENTITY, JACOBIAN A FULL MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		i__2 = *ns;
		for (k = 1; k <= i__2; ++k) {
			sum += dd[k] * zz[i__ + (k - 1) * *n];
		}
		ff[i__ + *n] = sum / *h__;
		cont[i__] = ff[i__ + *n] + y0[i__];
	}
	sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
	goto L77;

L11:
	/* ------  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		i__2 = *ns;
		for (k = 1; k <= i__2; ++k) {
			sum += dd[k] * zz[i__ + (k - 1) * *n];
		}
		ff[i__ + *n] = sum / *h__;
		cont[i__] = ff[i__ + *n] + y0[i__];
	}
L48:
	mm = *m1 / *m2;
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
		sum1 = 0;
		for (k = mm - 1; k >= 0; --k) {
			sum1 = (cont[j + k * *m2] + sum1) / *fac1;
			i__2 = *nm1;
			for (i__ = 1; i__ <= i__2; ++i__) {
				im1 = i__ + *m1;
				cont[im1] += fjac[i__ + (j + k * *m2) * fjac_dim1] * sum1;
			}
		}
	}
	sol_(nm1, lde1, &e1[e1_offset], &cont[*m1 + 1], &ip1[1]);
	for (i__ = *m1; i__ >= 1; --i__) {
		cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
	}
	goto L77;

L2:
	/* ------  B=IDENTITY, JACOBIAN A BANDED MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		i__2 = *ns;
		for (k = 1; k <= i__2; ++k) {
			sum += dd[k] * zz[i__ + (k - 1) * *n];
		}
		ff[i__ + *n] = sum / *h__;
		cont[i__] = ff[i__ + *n] + y0[i__];
	}
	solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &ip1[
		1]);
	goto L77;

L12:
	/* ------  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		i__2 = *ns;
		for (k = 1; k <= i__2; ++k) {
			sum += dd[k] * zz[i__ + (k - 1) * *n];
		}
		ff[i__ + *n] = sum / *h__;
		cont[i__] = ff[i__ + *n] + y0[i__];
	}
L45:
	mm = *m1 / *m2;
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
		sum1 = 0;
		for (k = mm - 1; k >= 0; --k) {
			sum1 = (cont[j + k * *m2] + sum1) / *fac1;
			/* Computing MAX */
			i__2 = 1, i__3 = j - *mujac;
			/* Computing MIN */
			i__5 = *nm1, i__6 = j + *mljac;
			i__4 = MIN(i__5, i__6);
			for (i__ = MAX(i__2, i__3); i__ <= i__4; ++i__) {
				im1 = i__ + *m1;
				cont[im1] += fjac[i__ + *mujac + 1 - j + (j + k * *m2) *
					fjac_dim1] * sum1;
			}
		}
	}
	solb_(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[*m1 +
		1], &ip1[1]);
	for (i__ = *m1; i__ >= 1; --i__) {
		cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
	}
	goto L77;

L3:
	/* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		i__4 = *ns;
		for (k = 1; k <= i__4; ++k) {
			sum += dd[k] * zz[i__ + (k - 1) * *n];
		}
		ff[i__] = sum / *h__;
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		/* Computing MAX */
		i__4 = 1, i__2 = i__ - *mlmas;
		/* Computing MIN */
		i__5 = *n, i__6 = i__ + *mumas;
		i__3 = MIN(i__5, i__6);
		for (j = MAX(i__4, i__2); j <= i__3; ++j) {
			sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ff[j];
		}
		ff[i__ + *n] = sum;
		cont[i__] = sum + y0[i__];
	}
	sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
	goto L77;

L13:
	/* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
	i__1 = *m1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		i__3 = *ns;
		for (k = 1; k <= i__3; ++k) {
			sum += dd[k] * zz[i__ + (k - 1) * *n];
		}
		ff[i__ + *n] = sum / *h__;
		cont[i__] = ff[i__ + *n] + y0[i__];
	}
	i__1 = *n;
	for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
		sum = 0;
		i__3 = *ns;
		for (k = 1; k <= i__3; ++k) {
			sum += dd[k] * zz[i__ + (k - 1) * *n];
		}
		ff[i__] = sum / *h__;
	}
	i__1 = *nm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		/* Computing MAX */
		i__3 = 1, i__4 = i__ - *mlmas;
		/* Computing MIN */
		i__5 = *nm1, i__6 = i__ + *mumas;
		i__2 = MIN(i__5, i__6);
		for (j = MAX(i__3, i__4); j <= i__2; ++j) {
			sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ff[j + *
				m1];
		}
		im1 = i__ + *m1;
		ff[im1 + *n] = sum;
		cont[im1] = sum + y0[im1];
	}
	goto L48;

L4:
	/* ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		i__2 = *ns;
		for (k = 1; k <= i__2; ++k) {
			sum += dd[k] * zz[i__ + (k - 1) * *n];
		}
		ff[i__] = sum / *h__;
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		/* Computing MAX */
		i__2 = 1, i__3 = i__ - *mlmas;
		/* Computing MIN */
		i__5 = *n, i__6 = i__ + *mumas;
		i__4 = MIN(i__5, i__6);
		for (j = MAX(i__2, i__3); j <= i__4; ++j) {
			sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ff[j];
		}
		ff[i__ + *n] = sum;
		cont[i__] = sum + y0[i__];
	}
	solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &ip1[
		1]);
	goto L77;

L14:
	/* ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER */
	i__1 = *m1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		i__4 = *ns;
		for (k = 1; k <= i__4; ++k) {
			sum += dd[k] * zz[i__ + (k - 1) * *n];
		}
		ff[i__ + *n] = sum / *h__;
		cont[i__] = ff[i__ + *n] + y0[i__];
	}
	i__1 = *n;
	for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
		sum = 0;
		i__4 = *ns;
		for (k = 1; k <= i__4; ++k) {
			sum += dd[k] * zz[i__ + (k - 1) * *n];
		}
		ff[i__] = sum / *h__;
	}
	i__1 = *nm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		/* Computing MAX */
		i__4 = 1, i__2 = i__ - *mlmas;
		/* Computing MIN */
		i__5 = *nm1, i__6 = i__ + *mumas;
		i__3 = MIN(i__5, i__6);
		for (j = MAX(i__4, i__2); j <= i__3; ++j) {
			sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ff[j + *
				m1];
		}
		im1 = i__ + *m1;
		ff[im1 + *n] = sum;
		cont[im1] = sum + y0[im1];
	}
	goto L45;

L5:
	/* ------  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		i__3 = *ns;
		for (k = 1; k <= i__3; ++k) {
			sum += dd[k] * zz[i__ + (k - 1) * *n];
		}
		ff[i__] = sum / *h__;
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		i__3 = *n;
		for (j = 1; j <= i__3; ++j) {
			sum += fmas[i__ + j * fmas_dim1] * ff[j];
		}
		ff[i__ + *n] = sum;
		cont[i__] = sum + y0[i__];
	}
	sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
	goto L77;

L15:
	/* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
	i__1 = *m1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		i__3 = *ns;
		for (k = 1; k <= i__3; ++k) {
			sum += dd[k] * zz[i__ + (k - 1) * *n];
		}
		ff[i__ + *n] = sum / *h__;
		cont[i__] = ff[i__ + *n] + y0[i__];
	}
	i__1 = *n;
	for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
		sum = 0;
		i__3 = *ns;
		for (k = 1; k <= i__3; ++k) {
			sum += dd[k] * zz[i__ + (k - 1) * *n];
		}
		ff[i__] = sum / *h__;
	}
	i__1 = *nm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		i__3 = *nm1;
		for (j = 1; j <= i__3; ++j) {
			sum += fmas[i__ + j * fmas_dim1] * ff[j + *m1];
		}
		im1 = i__ + *m1;
		ff[im1 + *n] = sum;
		cont[im1] = sum + y0[im1];
	}
	goto L48;

L6:
	/* ------  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
	/* ------  THIS OPTION IS NOT PROVIDED */
	return 0;

L7:
	/* ------  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = 0;
		i__3 = *ns;
		for (k = 1; k <= i__3; ++k) {
			sum += dd[k] * zz[i__ + (k - 1) * *n];
		}
		ff[i__ + *n] = sum / *h__;
		cont[i__] = ff[i__ + *n] + y0[i__];
	}
	for (mm = *n - 2; mm >= 1; --mm) {
		mp = *n - mm;
		i__ = iphes[mp];
		if (i__ == mp) {
			goto L310;
		}
		zsafe = cont[mp];
		cont[mp] = cont[i__];
		cont[i__] = zsafe;
	L310:
		i__1 = *n;
		for (i__ = mp + 1; i__ <= i__1; ++i__) {
			cont[i__] -= fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
		}
	}
	solh_(n, lde1, &e1[e1_offset], &c__1, &cont[1], &ip1[1]);
	i__1 = *n - 2;
	for (mm = 1; mm <= i__1; ++mm) {
		mp = *n - mm;
		i__3 = *n;
		for (i__ = mp + 1; i__ <= i__3; ++i__) {
			cont[i__] += fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
		}
		i__ = iphes[mp];
		if (i__ == mp) {
			goto L440;
		}
		zsafe = cont[mp];
		cont[mp] = cont[i__];
		cont[i__] = zsafe;
	L440:
		;
	}

	/* -------------------------------------- */

L77:
	*err = 0;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		/* Computing 2nd power */
		d__1 = cont[i__] / scal[i__];
		*err += d__1 * d__1;
	}
	/* Computing MAX */
	d__1 = SQRT(*err / *n);
	*err = MAX(d__1, (typeRNum) 1e-10);

	if (*err < 1) {
		return 0;
	}
	if (*first || *reject) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			cont[i__] = y[i__] + cont[i__];
		}
		(*fcn)(n, x, &cont[1], &ff[1], tvec, xvec, uvec, pvec, dcdxvec, grampc, pfct);
		++(*nfcn);
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			cont[i__] = ff[i__] + ff[i__ + *n];
		}
		switch (*ijob) {
		case 1:  goto L31;
		case 2:  goto L32;
		case 3:  goto L31;
		case 4:  goto L32;
		case 5:  goto L31;
		case 6:  goto L32;
		case 7:  goto L33;
		case 8:  goto L55;
		case 9:  goto L55;
		case 10:  goto L55;
		case 11:  goto L41;
		case 12:  goto L42;
		case 13:  goto L41;
		case 14:  goto L42;
		case 15:  goto L41;
		}
		/* ------ FULL MATRIX OPTION */
	L31:
		sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
		goto L88;
		/* ------ FULL MATRIX OPTION, SECOND ORDER */
	L41:
		i__1 = *m2;
		for (j = 1; j <= i__1; ++j) {
			sum1 = 0;
			for (k = mm - 1; k >= 0; --k) {
				sum1 = (cont[j + k * *m2] + sum1) / *fac1;
				i__3 = *nm1;
				for (i__ = 1; i__ <= i__3; ++i__) {
					im1 = i__ + *m1;
					cont[im1] += fjac[i__ + (j + k * *m2) * fjac_dim1] * sum1;
				}
			}
		}
		sol_(nm1, lde1, &e1[e1_offset], &cont[*m1 + 1], &ip1[1]);
		for (i__ = *m1; i__ >= 1; --i__) {
			cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
		}
		goto L88;
		/* ------ BANDED MATRIX OPTION */
	L32:
		solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &
			ip1[1]);
		goto L88;
		/* ------ BANDED MATRIX OPTION, SECOND ORDER */
	L42:
		i__1 = *m2;
		for (j = 1; j <= i__1; ++j) {
			sum1 = 0;
			for (k = mm - 1; k >= 0; --k) {
				sum1 = (cont[j + k * *m2] + sum1) / *fac1;
				/* Computing MAX */
				i__3 = 1, i__4 = j - *mujac;
				/* Computing MIN */
				i__5 = *nm1, i__6 = j + *mljac;
				i__2 = MIN(i__5, i__6);
				for (i__ = MAX(i__3, i__4); i__ <= i__2; ++i__) {
					im1 = i__ + *m1;
					cont[im1] += fjac[i__ + *mujac + 1 - j + (j + k * *m2) *
						fjac_dim1] * sum1;
				}
			}
		}
		solb_(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[*
			m1 + 1], &ip1[1]);
		for (i__ = *m1; i__ >= 1; --i__) {
			cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
		}
		goto L88;
		/* ------ HESSENBERG MATRIX OPTION */
	L33:
		for (mm = *n - 2; mm >= 1; --mm) {
			mp = *n - mm;
			i__ = iphes[mp];
			if (i__ == mp) {
				goto L510;
			}
			zsafe = cont[mp];
			cont[mp] = cont[i__];
			cont[i__] = zsafe;
		L510:
			i__1 = *n;
			for (i__ = mp + 1; i__ <= i__1; ++i__) {
				cont[i__] -= fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
			}
		}
		solh_(n, lde1, &e1[e1_offset], &c__1, &cont[1], &ip1[1]);
		i__1 = *n - 2;
		for (mm = 1; mm <= i__1; ++mm) {
			mp = *n - mm;
			i__2 = *n;
			for (i__ = mp + 1; i__ <= i__2; ++i__) {
				cont[i__] += fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
			}
			i__ = iphes[mp];
			if (i__ == mp) {
				goto L640;
			}
			zsafe = cont[mp];
			cont[mp] = cont[i__];
			cont[i__] = zsafe;
		L640:
			;
		}
		/* ----------------------------------- */
	L88:
		*err = 0;
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			/* Computing 2nd power */
			d__1 = cont[i__] / scal[i__];
			*err += d__1 * d__1;
		}
		/* Computing MAX */
		d__1 = SQRT(*err / *n);
		*err = MAX(d__1, (typeRNum) 1e-10);
	}
	return 0;

	/* ----------------------------------------------------------- */

L55:
	return 0;
} /* estrav_ */


/*     END OF SUBROUTINE ESTRAV */

/* *********************************************************** */

/* Subroutine */ int slvrod_(typeLInt *n, typeRNum *fjac, typeLInt *ldjac,
	typeLInt *mljac, typeLInt *mujac, typeRNum *fmas, typeLInt *ldmas,
	typeLInt *mlmas, typeLInt *mumas, typeLInt *m1, typeLInt *m2, typeLInt *
	nm1, typeRNum *fac1, typeRNum *e, typeLInt *lde, typeLInt *ip,
	typeRNum *dy, typeRNum *ak, typeRNum *fx, typeRNum *ynew,
	typeRNum *hd, typeLInt *ijob, typeLogical *stage1)
{
	/* System generated locals */
	typeLInt fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e_dim1, e_offset,
		i__1, i__2, i__3, i__4, i__5, i__6;

	/* Local variables */
	static typeLInt i__, j, k, mm, im1, jkm;
	extern /* Subroutine */ int sol_(typeLInt *, typeLInt *, typeRNum *,
		typeRNum *, typeLInt *);
	static typeRNum sum;
	extern /* Subroutine */ int solb_(typeLInt *, typeLInt *, typeRNum *,
		typeLInt *, typeLInt *, typeRNum *, typeLInt *);


	/* Parameter adjustments */
	--ynew;
	--fx;
	--ak;
	--dy;
	fjac_dim1 = *ldjac;
	fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;
	--ip;
	fmas_dim1 = *ldmas;
	fmas_offset = 1 + fmas_dim1;
	fmas -= fmas_offset;
	e_dim1 = *lde;
	e_offset = 1 + e_dim1;
	e -= e_offset;

	/* Function Body */
	if (*hd == 0) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			ak[i__] = dy[i__];
		}
	}
	else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			ak[i__] = dy[i__] + *hd * fx[i__];
		}
	}

	switch (*ijob) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	case 4:  goto L4;
	case 5:  goto L5;
	case 6:  goto L6;
	case 7:  goto L55;
	case 8:  goto L55;
	case 9:  goto L55;
	case 10:  goto L55;
	case 11:  goto L11;
	case 12:  goto L12;
	case 13:  goto L13;
	case 14:  goto L13;
	case 15:  goto L15;
	}

	/* ----------------------------------------------------------- */

L1:
	/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
	if (*stage1) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			ak[i__] += ynew[i__];
		}
	}
	sol_(n, lde, &e[e_offset], &ak[1], &ip[1]);
	return 0;

	/* ----------------------------------------------------------- */

L11:
	/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
	if (*stage1) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			ak[i__] += ynew[i__];
		}
	}
L48:
	mm = *m1 / *m2;
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
		sum = 0;
		for (k = mm - 1; k >= 0; --k) {
			jkm = j + k * *m2;
			sum = (ak[jkm] + sum) / *fac1;
			i__2 = *nm1;
			for (i__ = 1; i__ <= i__2; ++i__) {
				im1 = i__ + *m1;
				ak[im1] += fjac[i__ + jkm * fjac_dim1] * sum;
			}
		}
	}
	sol_(nm1, lde, &e[e_offset], &ak[*m1 + 1], &ip[1]);
	for (i__ = *m1; i__ >= 1; --i__) {
		ak[i__] = (ak[i__] + ak[*m2 + i__]) / *fac1;
	}
	return 0;

	/* ----------------------------------------------------------- */

L2:
	/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
	if (*stage1) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			ak[i__] += ynew[i__];
		}
	}
	solb_(n, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &ak[1], &ip[1]);
	return 0;

	/* ----------------------------------------------------------- */

L12:
	/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
	if (*stage1) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			ak[i__] += ynew[i__];
		}
	}
L45:
	mm = *m1 / *m2;
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
		sum = 0;
		for (k = mm - 1; k >= 0; --k) {
			jkm = j + k * *m2;
			sum = (ak[jkm] + sum) / *fac1;
			/* Computing MAX */
			i__2 = 1, i__3 = j - *mujac;
			/* Computing MIN */
			i__5 = *nm1, i__6 = j + *mljac;
			i__4 = MIN(i__5, i__6);
			for (i__ = MAX(i__2, i__3); i__ <= i__4; ++i__) {
				im1 = i__ + *m1;
				ak[im1] += fjac[i__ + *mujac + 1 - j + jkm * fjac_dim1] * sum;
			}
		}
	}
	solb_(nm1, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &ak[*m1 + 1], &
		ip[1]);
	for (i__ = *m1; i__ >= 1; --i__) {
		ak[i__] = (ak[i__] + ak[*m2 + i__]) / *fac1;
	}
	return 0;

	/* ----------------------------------------------------------- */

L3:
	/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
	if (*stage1) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			sum = 0;
			/* Computing MAX */
			i__4 = 1, i__2 = i__ - *mlmas;
			/* Computing MIN */
			i__5 = *n, i__6 = i__ + *mumas;
			i__3 = MIN(i__5, i__6);
			for (j = MAX(i__4, i__2); j <= i__3; ++j) {
				sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ynew[
					j];
			}
			ak[i__] += sum;
		}
	}
	sol_(n, lde, &e[e_offset], &ak[1], &ip[1]);
	return 0;

	/* ----------------------------------------------------------- */

L13:
	/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
	if (*stage1) {
		i__1 = *m1;
		for (i__ = 1; i__ <= i__1; ++i__) {
			ak[i__] += ynew[i__];
		}
		i__1 = *nm1;
		for (i__ = 1; i__ <= i__1; ++i__) {
			sum = 0;
			/* Computing MAX */
			i__3 = 1, i__4 = i__ - *mlmas;
			/* Computing MIN */
			i__5 = *nm1, i__6 = i__ + *mumas;
			i__2 = MIN(i__5, i__6);
			for (j = MAX(i__3, i__4); j <= i__2; ++j) {
				sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ynew[
					j + *m1];
			}
			im1 = i__ + *m1;
			ak[im1] += sum;
		}
	}
	if (*ijob == 14) {
		goto L45;
	}
	goto L48;

	/* ----------------------------------------------------------- */

L4:
	/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
	if (*stage1) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			sum = 0;
			/* Computing MAX */
			i__2 = 1, i__3 = i__ - *mlmas;
			/* Computing MIN */
			i__5 = *n, i__6 = i__ + *mumas;
			i__4 = MIN(i__5, i__6);
			for (j = MAX(i__2, i__3); j <= i__4; ++j) {
				sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ynew[
					j];
			}
			ak[i__] += sum;
		}
	}
	solb_(n, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &ak[1], &ip[1]);
	return 0;

	/* ----------------------------------------------------------- */

L5:
	/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
	if (*stage1) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			sum = 0;
			i__4 = *n;
			for (j = 1; j <= i__4; ++j) {
				sum += fmas[i__ + j * fmas_dim1] * ynew[j];
			}
			ak[i__] += sum;
		}
	}
	sol_(n, lde, &e[e_offset], &ak[1], &ip[1]);
	return 0;

	/* ----------------------------------------------------------- */

L15:
	/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
	if (*stage1) {
		i__1 = *m1;
		for (i__ = 1; i__ <= i__1; ++i__) {
			ak[i__] += ynew[i__];
		}
		i__1 = *nm1;
		for (i__ = 1; i__ <= i__1; ++i__) {
			sum = 0;
			i__4 = *nm1;
			for (j = 1; j <= i__4; ++j) {
				sum += fmas[i__ + j * fmas_dim1] * ynew[j + *m1];
			}
			im1 = i__ + *m1;
			ak[im1] += sum;
		}
	}
	goto L48;

	/* ----------------------------------------------------------- */

L6:
	/* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
	/* ---  THIS OPTION IS NOT PROVIDED */
	if (*stage1) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			sum = 0;
			i__4 = *n;
			for (j = 1; j <= i__4; ++j) {
				/* L623: */
				sum += fmas[i__ + j * fmas_dim1] * ynew[j];
			}
			/* L624: */
			ak[i__] += sum;
		}
		solb_(n, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &ak[1], &ip[1]
		);
	}
	return 0;

	/* ----------------------------------------------------------- */

L55:
	return 0;
} /* slvrod_ */


/*     END OF SUBROUTINE SLVROD */


/* *********************************************************** */

/* Subroutine */ int slvseu_(typeLInt *n, typeRNum *fjac, typeLInt *ldjac,
	typeLInt *mljac, typeLInt *mujac, typeRNum *fmas, typeLInt *ldmas,
	typeLInt *mlmas, typeLInt *mumas, typeLInt *m1, typeLInt *m2, typeLInt *
	nm1, typeRNum *fac1, typeRNum *e, typeLInt *lde, typeLInt *ip,
	typeLInt *iphes, typeRNum *del, typeLInt *ijob)
{
	/* System generated locals */
	typeLInt fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e_dim1, e_offset,
		i__1, i__2, i__3, i__4, i__5, i__6;

	/* Local variables */
	static typeLInt i__, j, k, mm, mp, im1, mp1, jkm, mmm;
	extern /* Subroutine */ int sol_(typeLInt *, typeLInt *, typeRNum *,
		typeRNum *, typeLInt *);
	static typeRNum sum;
	extern /* Subroutine */ int solb_(typeLInt *, typeLInt *, typeRNum *,
		typeLInt *, typeLInt *, typeRNum *, typeLInt *), solh_(typeLInt *,
			typeLInt *, typeRNum *, typeLInt *, typeRNum *, typeLInt *);
	static typeRNum zsafe;


	/* Parameter adjustments */
	--del;
	--iphes;
	fjac_dim1 = *ldjac;
	fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;
	--ip;
	fmas_dim1 = *ldmas;
	fmas_offset = 1 + fmas_dim1;
	fmas -= fmas_offset;
	e_dim1 = *lde;
	e_offset = 1 + e_dim1;
	e -= e_offset;

	/* Function Body */
	switch (*ijob) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L1;
	case 4:  goto L2;
	case 5:  goto L1;
	case 6:  goto L55;
	case 7:  goto L7;
	case 8:  goto L55;
	case 9:  goto L55;
	case 10:  goto L55;
	case 11:  goto L11;
	case 12:  goto L12;
	case 13:  goto L11;
	case 14:  goto L12;
	case 15:  goto L11;
	}

	/* ----------------------------------------------------------- */

L1:
	/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
	sol_(n, lde, &e[e_offset], &del[1], &ip[1]);
	return 0;

	/* ----------------------------------------------------------- */

L11:
	/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
	mm = *m1 / *m2;
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
		sum = 0;
		for (k = mm - 1; k >= 0; --k) {
			jkm = j + k * *m2;
			sum = (del[jkm] + sum) / *fac1;
			i__2 = *nm1;
			for (i__ = 1; i__ <= i__2; ++i__) {
				im1 = i__ + *m1;
				del[im1] += fjac[i__ + jkm * fjac_dim1] * sum;
			}
		}
	}
	sol_(nm1, lde, &e[e_offset], &del[*m1 + 1], &ip[1]);
	for (i__ = *m1; i__ >= 1; --i__) {
		del[i__] = (del[i__] + del[*m2 + i__]) / *fac1;
	}
	return 0;

	/* ----------------------------------------------------------- */

L2:
	/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
	solb_(n, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &del[1], &ip[1]);
	return 0;

	/* ----------------------------------------------------------- */

L12:
	/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
	mm = *m1 / *m2;
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
		sum = 0;
		for (k = mm - 1; k >= 0; --k) {
			jkm = j + k * *m2;
			sum = (del[jkm] + sum) / *fac1;
			/* Computing MAX */
			i__2 = 1, i__3 = j - *mujac;
			/* Computing MIN */
			i__5 = *nm1, i__6 = j + *mljac;
			i__4 = MIN(i__5, i__6);
			for (i__ = MAX(i__2, i__3); i__ <= i__4; ++i__) {
				im1 = i__ + *m1;
				del[im1] += fjac[i__ + *mujac + 1 - j + jkm * fjac_dim1] *
					sum;
			}
		}
	}
	solb_(nm1, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &del[*m1 + 1], &
		ip[1]);
	for (i__ = *m1; i__ >= 1; --i__) {
		del[i__] = (del[i__] + del[*m2 + i__]) / *fac1;
	}
	return 0;

	/* ----------------------------------------------------------- */

L7:
	/* ---  HESSENBERG OPTION */
	for (mmm = *n - 2; mmm >= 1; --mmm) {
		mp = *n - mmm;
		mp1 = mp - 1;
		i__ = iphes[mp];
		if (i__ == mp) {
			goto L110;
		}
		zsafe = del[mp];
		del[mp] = del[i__];
		del[i__] = zsafe;
	L110:
		i__1 = *n;
		for (i__ = mp + 1; i__ <= i__1; ++i__) {
			del[i__] -= fjac[i__ + mp1 * fjac_dim1] * del[mp];
		}
	}
	solh_(n, lde, &e[e_offset], &c__1, &del[1], &ip[1]);
	i__1 = *n - 2;
	for (mmm = 1; mmm <= i__1; ++mmm) {
		mp = *n - mmm;
		mp1 = mp - 1;
		i__4 = *n;
		for (i__ = mp + 1; i__ <= i__4; ++i__) {
			del[i__] += fjac[i__ + mp1 * fjac_dim1] * del[mp];
		}
		i__ = iphes[mp];
		if (i__ == mp) {
			goto L240;
		}
		zsafe = del[mp];
		del[mp] = del[i__];
		del[i__] = zsafe;
	L240:
		;
	}
	return 0;

	/* ----------------------------------------------------------- */

L55:
	return 0;
} /* slvseu_ */


/*     END OF SUBROUTINE SLVSEU */

/* Subroutine */ int rodas_(typeLInt *n, U_fp fcn, typeLInt *ifcn,
	typeRNum * x, typeRNum *y, typeRNum *xend, typeRNum *h__, typeRNum *rtol, typeRNum *atol, typeLInt *itol,
	U_fp jac, typeLInt *ijac, typeLInt *mljac, typeLInt *mujac,
	U_fp dfx, typeLInt *idfx,
	U_fp mas, typeLInt *	imas, typeLInt *mlmas, typeLInt *mumas,
	U_fp solout, typeLInt *iout,
	typeRNum *work, typeLInt *lwork, typeLInt *iwork, typeLInt *liwork,
	ctypeRNum *tvec, ctypeRNum *xvec, ctypeRNum *uvec, ctypeRNum *pvec, ctypeRNum *dcdxvec, const typeGRAMPC *grampc, const typeffctPtr pfct, typeLInt *idid)
{
	/* System generated locals */
	typeLInt i__1;

	/* Builtin functions */
	typeLInt s_wsle(cilist *), do_lio(typeLInt *, typeLInt *, char *, ftnlen),
		e_wsle(void);

	/* Local variables */
	static typeLInt i__, m1, m2, nm1, iee, lde;
	static typeRNum fac1, fac2;
	static typeLInt ndec, njac;
	static typeRNum safe;
	static typeLInt ijob, nfcn, ieip;
	static typeLogical pred;
	static typeLInt meth;
	static typeRNum hmax;
	static typeLInt iedy, iefx, nmax, nsol, ieak1, ieak2, ieak3, ieak4, ieak5,
		ieak6, iedy1, iejac, ldjac;
	static typeLogical jband;
	static typeLInt iecon, iemas, ldmas;
	static typeLogical arret;
	static typeLInt nstep, ldmas2, naccpt, nrejct;
	static typeLogical implct;
	static typeLInt ieynew, istore;
	static typeLogical autnms;
	extern /* Subroutine */ int roscor_(typeLInt *, U_fp, typeRNum *,
		typeRNum *, typeRNum *, typeRNum *, typeRNum *,
		typeRNum *, typeRNum *, typeLInt *, U_fp, typeLInt *, typeLInt *,
		typeLInt *, U_fp, typeLInt *, U_fp, typeLInt *, typeLInt *, U_fp,
		typeLInt *, typeLInt *, typeLInt *, typeRNum *, typeLInt *, typeLInt
		*, typeRNum *, typeRNum *, typeRNum *, typeLogical *, typeLogical *,
		typeLogical *, typeLogical *, typeLInt *, typeLInt *, typeLInt *,
		typeRNum *, typeRNum *, typeRNum *, typeRNum *,
		typeRNum *, typeRNum *, typeRNum *, typeRNum *,
		typeRNum *, typeRNum *, typeRNum *, typeRNum *,
		typeRNum *, typeLInt *, typeRNum *, typeLInt *, typeLInt *,
		typeLInt *, typeLInt *, typeLInt *, typeLInt *, typeLInt *, typeLInt *,
		typeLInt *, typeLInt *, ctypeRNum *, ctypeRNum *, ctypeRNum *, ctypeRNum *, ctypeRNum *, const typeGRAMPC *, const typeffctPtr);
	static typeRNum uround;


	/* Fortran I/O blocks */
	static cilist io___282 = { 0, 6, 0, 0, 0 };
	static cilist io___284 = { 0, 6, 0, 0, 0 };
	static cilist io___289 = { 0, 6, 0, 0, 0 };
	static cilist io___291 = { 0, 6, 0, 0, 0 };
	static cilist io___295 = { 0, 6, 0, 0, 0 };
	static cilist io___297 = { 0, 6, 0, 0, 0 };
	static cilist io___298 = { 0, 6, 0, 0, 0 };
	static cilist io___300 = { 0, 6, 0, 0, 0 };
	static cilist io___308 = { 0, 6, 0, 0, 0 };
	static cilist io___325 = { 0, 6, 0, 0, 0 };
	static cilist io___327 = { 0, 6, 0, 0, 0 };


	/* ---------------------------------------------------------- */
	/*     NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC) */
	/*     SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS  MY'=F(X,Y). */
	/*     THIS IS AN EMBEDDED ROSENBROCK METHOD OF ORDER (3)4 */
	/*     (WITH STEP SIZE CONTROL). */
	/*     C.F. SECTIONS IV.7  AND VI.3 */

	/*     AUTHORS: E. HAIRER AND G. WANNER */
	/*              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES */
	/*              CH-1211 GENEVE 24, SWITZERLAND */
	/*              E-MAIL:  Ernst.Hairer@math.unige.ch */
	/*                       Gerhard.Wanner@math.unige.ch */

	/*     THIS CODE IS PART OF THE BOOK: */
	/*         E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL */
	/*         EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS. */
	/*         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS 14, */
	/*         SPRINGER-VERLAG 1991, SECOND EDITION 1996. */

	/*     VERSION OF OCTOBER 28, 1996 */

	/*     INPUT PARAMETERS */
	/*     ---------------- */
	/*     N           DIMENSION OF THE SYSTEM */

	/*     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE */
	/*                 VALUE OF F(X,Y): */
	/*                    SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR) */
	/*                    DOUBLE PRECISION X,Y(N),F(N) */
	/*                    F(1)=...   ETC. */
	/*                 RPAR, IPAR (SEE BELOW) */

	/*     IFCN        GIVES INFORMATION ON FCN: */
	/*                    IFCN=0: F(X,Y) INDEPENDENT OF X (AUTONOMOUS) */
	/*                    IFCN=1: F(X,Y) MAY DEPEND ON X (NON-AUTONOMOUS) */

	/*     X           INITIAL X-VALUE */

	/*     Y(N)        INITIAL VALUES FOR Y */

	/*     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE) */

	/*     H           INITIAL STEP SIZE GUESS; */
	/*                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT, */
	/*                 H=1.D0/(NORM OF F'), USUALLY 1.D-2 OR 1.D-3, IS GOOD. */
	/*                 THIS CHOICE IS NOT VERY IMPORTANT, THE CODE QUICKLY */
	/*                 ADAPTS ITS STEP SIZE (IF H=0.D0, THE CODE PUTS H=1.D-6). */

	/*     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY */
	/*                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N. */

	/*     ITOL        SWITCH FOR RTOL AND ATOL: */
	/*                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS. */
	/*                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF */
	/*                     Y(I) BELOW RTOL*DABS(Y(I))+ATOL */
	/*                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS. */
	/*                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW */
	/*                     RTOL(I)*DABS(Y(I))+ATOL(I). */

	/*     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES */
	/*                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO Y */
	/*                 (THIS ROUTINE IS ONLY CALLED IF IJAC=1; SUPPLY */
	/*                 A DUMMY SUBROUTINE IN THE CASE IJAC=0). */
	/*                 FOR IJAC=1, THIS SUBROUTINE MUST HAVE THE FORM */
	/*                    SUBROUTINE JAC(N,X,Y,DFY,LDFY,RPAR,IPAR) */
	/*                    DOUBLE PRECISION X,Y(N),DFY(LDFY,N) */
	/*                    DFY(1,1)= ... */
	/*                 LDFY, THE COLOMN-LENGTH OF THE ARRAY, IS */
	/*                 FURNISHED BY THE CALLING PROGRAM. */
	/*                 IF (MLJAC.EQ.N) THE JACOBIAN IS SUPPOSED TO */
	/*                    BE FULL AND THE PARTIAL DERIVATIVES ARE */
	/*                    STORED IN DFY AS */
	/*                       DFY(I,J) = PARTIAL F(I) / PARTIAL Y(J) */
	/*                 ELSE, THE JACOBIAN IS TAKEN AS BANDED AND */
	/*                    THE PARTIAL DERIVATIVES ARE STORED */
	/*                    DIAGONAL-WISE AS */
	/*                       DFY(I-J+MUJAC+1,J) = PARTIAL F(I) / PARTIAL Y(J). */

	/*     IJAC        SWITCH FOR THE COMPUTATION OF THE JACOBIAN: */
	/*                    IJAC=0: JACOBIAN IS COMPUTED INTERNALLY BY FINITE */
	/*                       DIFFERENCES, SUBROUTINE "JAC" IS NEVER CALLED. */
	/*                    IJAC=1: JACOBIAN IS SUPPLIED BY SUBROUTINE JAC. */

	/*     MLJAC       SWITCH FOR THE BANDED STRUCTURE OF THE JACOBIAN: */
	/*                    MLJAC=N: JACOBIAN IS A FULL MATRIX. THE LINEAR */
	/*                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION. */
	/*                    0<=MLJAC<N: MLJAC IS THE LOWER BANDWITH OF JACOBIAN */
	/*                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW */
	/*                       THE MAIN DIAGONAL). */

	/*     MUJAC       UPPER BANDWITH OF JACOBIAN  MATRIX (>= NUMBER OF NON- */
	/*                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL). */
	/*                 NEED NOT BE DEFINED IF MLJAC=N. */

	/*     DFX         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES */
	/*                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO X */
	/*                 (THIS ROUTINE IS ONLY CALLED IF IDFX=1 AND IFCN=1; */
	/*                 SUPPLY A DUMMY SUBROUTINE IN THE CASE IDFX=0 OR IFCN=0). */
	/*                 OTHERWISE, THIS SUBROUTINE MUST HAVE THE FORM */
	/*                    SUBROUTINE DFX(N,X,Y,FX,RPAR,IPAR) */
	/*                    DOUBLE PRECISION X,Y(N),FX(N) */
	/*                    FX(1)= ... */

	/*     IDFX        SWITCH FOR THE COMPUTATION OF THE DF/DX: */
	/*                    IDFX=0: DF/DX IS COMPUTED INTERNALLY BY FINITE */
	/*                       DIFFERENCES, SUBROUTINE "DFX" IS NEVER CALLED. */
	/*                    IDFX=1: DF/DX IS SUPPLIED BY SUBROUTINE DFX. */

	/*     ----   MAS,IMAS,MLMAS, AND MUMAS HAVE ANALOG MEANINGS      ----- */
	/*     ----   FOR THE "MASS MATRIX" (THE MATRIX "M" OF SECTION IV.8): - */

	/*     MAS         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE MASS- */
	/*                 MATRIX M. */
	/*                 IF IMAS=0, THIS MATRIX IS ASSUMED TO BE THE IDENTITY */
	/*                 MATRIX AND NEEDS NOT TO BE DEFINED; */
	/*                 SUPPLY A DUMMY SUBROUTINE IN THIS CASE. */
	/*                 IF IMAS=1, THE SUBROUTINE MAS IS OF THE FORM */
	/*                    SUBROUTINE MAS(N,AM,LMAS,RPAR,IPAR) */
	/*                    DOUBLE PRECISION AM(LMAS,N) */
	/*                    AM(1,1)= .... */
	/*                    IF (MLMAS.EQ.N) THE MASS-MATRIX IS STORED */
	/*                    AS FULL MATRIX LIKE */
	/*                         AM(I,J) = M(I,J) */
	/*                    ELSE, THE MATRIX IS TAKEN AS BANDED AND STORED */
	/*                    DIAGONAL-WISE AS */
	/*                         AM(I-J+MUMAS+1,J) = M(I,J). */

	/*     IMAS       GIVES INFORMATION ON THE MASS-MATRIX: */
	/*                    IMAS=0: M IS SUPPOSED TO BE THE IDENTITY */
	/*                       MATRIX, MAS IS NEVER CALLED. */
	/*                    IMAS=1: MASS-MATRIX  IS SUPPLIED. */

	/*     MLMAS       SWITCH FOR THE BANDED STRUCTURE OF THE MASS-MATRIX: */
	/*                    MLMAS=N: THE FULL MATRIX CASE. THE LINEAR */
	/*                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION. */
	/*                    0<=MLMAS<N: MLMAS IS THE LOWER BANDWITH OF THE */
	/*                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW */
	/*                       THE MAIN DIAGONAL). */
	/*                 MLMAS IS SUPPOSED TO BE .LE. MLJAC. */

	/*     MUMAS       UPPER BANDWITH OF MASS-MATRIX (>= NUMBER OF NON- */
	/*                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL). */
	/*                 NEED NOT BE DEFINED IF MLMAS=N. */
	/*                 MUMAS IS SUPPOSED TO BE .LE. MUJAC. */

	/*     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE */
	/*                 NUMERICAL SOLUTION DURING INTEGRATION. */
	/*                 IF IOUT=1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP. */
	/*                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0. */
	/*                 IT MUST HAVE THE FORM */
	/*                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N, */
	/*                                       RPAR,IPAR,IRTRN) */
	/*                    DOUBLE PRECISION X,Y(N),CONT(LRC) */
	/*                    .... */
	/*                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH */
	/*                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS */
	/*                    THE FIRST GRID-POINT). */
	/*                 "XOLD" IS THE PRECEEDING GRID-POINT. */
	/*                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN */
	/*                    IS SET <0, RODAS RETURNS TO THE CALLING PROGRAM. */

	/*          -----  CONTINUOUS OUTPUT: ----- */
	/*                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION */
	/*                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH */
	/*                 THE FUNCTION */
	/*                        >>>   CONTRO(I,S,CONT,LRC)   <<< */
	/*                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH */
	/*                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE */
	/*                 S SHOULD LIE IN THE INTERVAL [XOLD,X]. */

	/*     IOUT        GIVES INFORMATION ON THE SUBROUTINE SOLOUT: */
	/*                    IOUT=0: SUBROUTINE IS NEVER CALLED */
	/*                    IOUT=1: SUBROUTINE IS USED FOR OUTPUT */

	/*     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK". */
	/*                 SERVES AS WORKING SPACE FOR ALL VECTORS AND MATRICES. */
	/*                 "LWORK" MUST BE AT LEAST */
	/*                             N*(LJAC+LMAS+LE1+14)+20 */
	/*                 WHERE */
	/*                    LJAC=N              IF MLJAC=N (FULL JACOBIAN) */
	/*                    LJAC=MLJAC+MUJAC+1  IF MLJAC<N (BANDED JAC) */
	/*                 AND */
	/*                    LMAS=0              IF IMAS=0 */
	/*                    LMAS=N              IF IMAS=1 AND MLMAS=N (FULL) */
	/*                    LMAS=MLMAS+MUMAS+1  IF MLMAS<N (BANDED MASS-M) */
	/*                 AND */
	/*                    LE1=N               IF MLJAC=N (FULL JACOBIAN) */
	/*                    LE1=2*MLJAC+MUJAC+1 IF MLJAC<N (BANDED JAC). */
	/*                 IN THE USUAL CASE WHERE THE JACOBIAN IS FULL AND THE */
	/*                 MASS-MATRIX IS THE INDENTITY (IMAS=0), THE MINIMUM */
	/*                 STORAGE REQUIREMENT IS */
	/*                             LWORK = 2*N*N+14*N+20. */
	/*                 IF IWORK(9)=M1>0 THEN "LWORK" MUST BE AT LEAST */
	/*                          N*(LJAC+14)+(N-M1)*(LMAS+LE1)+20 */
	/*                 WHERE IN THE DEFINITIONS OF LJAC, LMAS AND LE1 THE */
	/*                 NUMBER N CAN BE REPLACED BY N-M1. */

	/*     LWORK       DECLARED LENGTH OF ARRAY "WORK". */

	/*     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK". */
	/*                 "LIWORK" MUST BE AT LEAST N+20. */

	/*     LIWORK      DECLARED LENGTH OF ARRAY "IWORK". */

	/*     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH */
	/*                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING */
	/*                 PROGRAM AND THE FCN, DFX, JAC, MAS, SOLOUT SUBROUTINES. */

	/* ---------------------------------------------------------------------- */

	/*     SOPHISTICATED SETTING OF PARAMETERS */
	/*     ----------------------------------- */
	/*              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK */
	/*              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),..,WORK(4) */
	/*              AS WELL AS IWORK(1),IWORK(2) DIFFERENT FROM ZERO. */
	/*              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES: */

	/*    IWORK(1)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS. */
	/*              THE DEFAULT VALUE (FOR IWORK(1)=0) IS 100000. */

	/*    IWORK(2)  SWITCH FOR THE CHOICE OF THE COEFFICIENTS */
	/*              IF IWORK(2).EQ.1  METHOD (SEE BOOK, PAGE 452) */
	/*              IF IWORK(2).EQ.2  SAME METHOD WITH DIFFERENT PARAMETERS */
	/*              IF IWORK(2).EQ.3  METHOD WITH COEFF. OF GERD STEINEBACH */
	/*              THE DEFAULT VALUE (FOR IWORK(2)=0) IS IWORK(2)=1. */

	/*    IWORK(3)  SWITCH FOR STEP SIZE STRATEGY */
	/*              IF IWORK(3).EQ.1  MOD. PREDICTIVE CONTROLLER (GUSTAFSSON) */
	/*              IF IWORK(3).EQ.2  CLASSICAL APPROACH */
	/*              THE DEFAULT VALUE (FOR IWORK(3)=0) IS IWORK(3)=1. */

	/*       IF THE DIFFERENTIAL SYSTEM HAS THE SPECIAL STRUCTURE THAT */
	/*            Y(I)' = Y(I+M2)   FOR  I=1,...,M1, */
	/*       WITH M1 A MULTIPLE OF M2, A SUBSTANTIAL GAIN IN COMPUTERTIME */
	/*       CAN BE ACHIEVED BY SETTING THE PARAMETERS IWORK(9) AND IWORK(10). */
	/*       E.G., FOR SECOND ORDER SYSTEMS P'=V, V'=G(P,V), WHERE P AND V ARE */
	/*       VECTORS OF DIMENSION N/2, ONE HAS TO PUT M1=M2=N/2. */
	/*       FOR M1>0 SOME OF THE INPUT PARAMETERS HAVE DIFFERENT MEANINGS: */
	/*       - JAC: ONLY THE ELEMENTS OF THE NON-TRIVIAL PART OF THE */
	/*              JACOBIAN HAVE TO BE STORED */
	/*              IF (MLJAC.EQ.N-M1) THE JACOBIAN IS SUPPOSED TO BE FULL */
	/*                 DFY(I,J) = PARTIAL F(I+M1) / PARTIAL Y(J) */
	/*                FOR I=1,N-M1 AND J=1,N. */
	/*              ELSE, THE JACOBIAN IS BANDED ( M1 = M2 * MM ) */
	/*                 DFY(I-J+MUJAC+1,J+K*M2) = PARTIAL F(I+M1) / PARTIAL Y(J+K*M2) */
	/*                FOR I=1,MLJAC+MUJAC+1 AND J=1,M2 AND K=0,MM. */
	/*       - MLJAC: MLJAC=N-M1: IF THE NON-TRIVIAL PART OF THE JACOBIAN IS FULL */
	/*                0<=MLJAC<N-M1: IF THE (MM+1) SUBMATRICES (FOR K=0,MM) */
	/*                     PARTIAL F(I+M1) / PARTIAL Y(J+K*M2),  I,J=1,M2 */
	/*                    ARE BANDED, MLJAC IS THE MAXIMAL LOWER BANDWIDTH */
	/*                    OF THESE MM+1 SUBMATRICES */
	/*       - MUJAC: MAXIMAL UPPER BANDWIDTH OF THESE MM+1 SUBMATRICES */
	/*                NEED NOT BE DEFINED IF MLJAC=N-M1 */
	/*       - MAS: IF IMAS=0 THIS MATRIX IS ASSUMED TO BE THE IDENTITY AND */
	/*              NEED NOT BE DEFINED. SUPPLY A DUMMY SUBROUTINE IN THIS CASE. */
	/*              IT IS ASSUMED THAT ONLY THE ELEMENTS OF RIGHT LOWER BLOCK OF */
	/*              DIMENSION N-M1 DIFFER FROM THAT OF THE IDENTITY MATRIX. */
	/*              IF (MLMAS.EQ.N-M1) THIS SUBMATRIX IS SUPPOSED TO BE FULL */
	/*                 AM(I,J) = M(I+M1,J+M1)     FOR I=1,N-M1 AND J=1,N-M1. */
	/*              ELSE, THE MASS MATRIX IS BANDED */
	/*                 AM(I-J+MUMAS+1,J) = M(I+M1,J+M1) */
	/*       - MLMAS: MLMAS=N-M1: IF THE NON-TRIVIAL PART OF M IS FULL */
	/*                0<=MLMAS<N-M1: LOWER BANDWIDTH OF THE MASS MATRIX */
	/*       - MUMAS: UPPER BANDWIDTH OF THE MASS MATRIX */
	/*                NEED NOT BE DEFINED IF MLMAS=N-M1 */

	/*    IWORK(9)  THE VALUE OF M1.  DEFAULT M1=0. */

	/*    IWORK(10) THE VALUE OF M2.  DEFAULT M2=M1. */

	/*    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16. */

	/*    WORK(2)   MAXIMAL STEP SIZE, DEFAULT XEND-X. */

	/*    WORK(3), WORK(4)   PARAMETERS FOR STEP SIZE SELECTION */
	/*              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION */
	/*                 WORK(3) <= HNEW/HOLD <= WORK(4) */
	/*              DEFAULT VALUES: WORK(3)=0.2D0, WORK(4)=6.D0 */

	/*    WORK(5)   THE SAFETY FACTOR IN STEP SIZE PREDICTION, */
	/*              DEFAULT 0.9D0. */

	/* ----------------------------------------------------------------------- */

	/*     OUTPUT PARAMETERS */
	/*     ----------------- */
	/*     X           X-VALUE WHERE THE SOLUTION IS COMPUTED */
	/*                 (AFTER SUCCESSFUL RETURN X=XEND) */

	/*     Y(N)        SOLUTION AT X */

	/*     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP */

	/*     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN: */
	/*                   IDID= 1  COMPUTATION SUCCESSFUL, */
	/*                   IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT) */
	/*                   IDID=-1  INPUT IS NOT CONSISTENT, */
	/*                   IDID=-2  LARGER NMAX IS NEEDED, */
	/*                   IDID=-3  STEP SIZE BECOMES TOO SMALL, */
	/*                   IDID=-4  MATRIX IS REPEATEDLY SINGULAR. */

	/*   IWORK(14)  NFCN    NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL */
	/*                      EVALUATION OF THE JACOBIAN ARE NOT COUNTED) */
	/*   IWORK(15)  NJAC    NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY */
	/*                      OR NUMERICALLY) */
	/*   IWORK(16)  NSTEP   NUMBER OF COMPUTED STEPS */
	/*   IWORK(17)  NACCPT  NUMBER OF ACCEPTED STEPS */
	/*   IWORK(18)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST), */
	/*                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED) */
	/*   IWORK(19)  NDEC    NUMBER OF LU-DECOMPOSITIONS (N-DIMENSIONAL MATRIX) */
	/*   IWORK(20)  NSOL    NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS */
	/* --------------------------------------------------------- */
	/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
	/*          DECLARATIONS */
	/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
	/* *** *** *** *** *** *** *** */
	/*        SETTING THE PARAMETERS */
	/* *** *** *** *** *** *** *** */
			/* Parameter adjustments */
	--y;
	--rtol;
	--atol;
	--work;
	--iwork;

	/* Function Body */
	nfcn = 0;
	naccpt = 0;
	nrejct = 0;
	nstep = 0;
	njac = 0;
	ndec = 0;
	nsol = 0;
	arret = FALSE_;
	/* -------- NMAX , THE MAXIMAL NUMBER OF STEPS ----- */
	if (iwork[1] == 0) {
		nmax = 100000;
	}
	else {
		nmax = iwork[1];
		if (nmax <= 0) {
			s_wsle(&io___282);
			do_lio(&c__9, &c__1, " WRONG INPUT IWORK(1)=", (ftnlen)22);
			do_lio(&c__3, &c__1, (char *)&iwork[1], (ftnlen)sizeof(typeLInt));
			e_wsle();
			arret = TRUE_;
		}
	}
	/* -------- METH   COEFFICIENTS OF THE METHOD */
	if (iwork[2] == 0) {
		meth = 1;
	}
	else {
		meth = iwork[2];
		if (meth <= 0 || meth >= 4) {
			s_wsle(&io___284);
			do_lio(&c__9, &c__1, " CURIOUS INPUT IWORK(2)=", (ftnlen)24);
			do_lio(&c__3, &c__1, (char *)&iwork[2], (ftnlen)sizeof(typeLInt));
			e_wsle();
			arret = TRUE_;
		}
	}
	/* -------- PRED   STEP SIZE CONTROL */
	if (iwork[3] <= 1) {
		pred = TRUE_;
	}
	else {
		pred = FALSE_;
	}
	/* -------- PARAMETER FOR SECOND ORDER EQUATIONS */
	m1 = iwork[9];
	m2 = iwork[10];
	nm1 = *n - m1;
	if (m1 == 0) {
		m2 = *n;
	}
	if (m2 == 0) {
		m2 = m1;
	}
	if (m1 < 0 || m2 < 0 || m1 + m2 > *n) {
		s_wsle(&io___289);
		do_lio(&c__9, &c__1, " CURIOUS INPUT FOR IWORK(9,10)=", (ftnlen)31);
		do_lio(&c__3, &c__1, (char *)&m1, (ftnlen)sizeof(typeLInt));
		do_lio(&c__3, &c__1, (char *)&m2, (ftnlen)sizeof(typeLInt));
		e_wsle();
		arret = TRUE_;
	}
	/* -------- UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0 */
	if (work[1] == 0) {
		uround = (typeRNum) 1e-16;
	}
	else {
		uround = work[1];
		if (uround < (typeRNum) 1e-16 || uround >= (typeRNum)1) {
			s_wsle(&io___291);
			do_lio(&c__9, &c__1, " COEFFICIENTS HAVE 16 DIGITS, UROUND=", (
				ftnlen)37);
			do_lio(&c__5, &c__1, (char *)&work[1], (ftnlen)sizeof(typeRNum))
				;
			e_wsle();
			arret = TRUE_;
		}
	}
	/* -------- MAXIMAL STEP SIZE */
	if (work[2] == 0) {
		hmax = *xend - *x;
	}
	else {
		hmax = work[2];
	}
	/* -------  FAC1,FAC2     PARAMETERS FOR STEP SIZE SELECTION */
	if (work[3] == 0) {
		fac1 = 5;
	}
	else {
		fac1 = (typeRNum) 1. / work[3];
	}
	if (work[4] == 0) {
		fac2 = (typeRNum).16666666666666666;
	}
	else {
		fac2 = (typeRNum) 1. / work[4];
	}
	if (fac1 < (typeRNum) 1. || fac2 >(typeRNum) 1) {
		s_wsle(&io___295);
		do_lio(&c__9, &c__1, " CURIOUS INPUT WORK(3,4)=", (ftnlen)25);
		do_lio(&c__5, &c__1, (char *)&work[3], (ftnlen)sizeof(typeRNum));
		do_lio(&c__5, &c__1, (char *)&work[4], (ftnlen)sizeof(typeRNum));
		e_wsle();
		arret = TRUE_;
	}
	/* --------- SAFE     SAFETY FACTOR IN STEP SIZE PREDICTION */
	if (work[5] == 0) {
		safe = (typeRNum).9;
	}
	else {
		safe = work[5];
		if (safe <= .001 || safe >= (typeRNum)1) {
			s_wsle(&io___297);
			do_lio(&c__9, &c__1, " CURIOUS INPUT FOR WORK(5)=", (ftnlen)27);
			do_lio(&c__5, &c__1, (char *)&work[5], (ftnlen)sizeof(typeRNum))
				;
			e_wsle();
			arret = TRUE_;
		}
	}
	/* --------- CHECK IF TOLERANCES ARE O.K. */
	if (*itol == 0) {
		if (atol[1] <= 0. || rtol[1] <= uround * 10) {
			s_wsle(&io___298);
			do_lio(&c__9, &c__1, " TOLERANCES ARE TOO SMALL", (ftnlen)25);
			e_wsle();
			arret = TRUE_;
		}
	}
	else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			if (atol[i__] <= 0. || rtol[i__] <= uround * 10) {
				s_wsle(&io___300);
				do_lio(&c__9, &c__1, " TOLERANCES(", (ftnlen)12);
				do_lio(&c__3, &c__1, (char *)&i__, (ftnlen)sizeof(typeLInt));
				do_lio(&c__9, &c__1, ") ARE TOO SMALL", (ftnlen)15);
				e_wsle();
				arret = TRUE_;
			}
		}
	}
	/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
	/*         COMPUTATION OF ARRAY ENTRIES */
	/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
	/* ---- AUTONOMOUS, IMPLICIT, BANDED OR NOT ? */
	autnms = *ifcn == 0;
	implct = *imas != 0;
	jband = *mljac < nm1;
	/* -------- COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS --- */
	/* -- JACOBIAN AND MATRIX E */
	if (jband) {
		ldjac = *mljac + *mujac + 1;
		lde = *mljac + ldjac;
	}
	else {
		*mljac = nm1;
		*mujac = nm1;
		ldjac = nm1;
		lde = nm1;
	}
	/* -- MASS MATRIX */
	if (implct) {
		if (*mlmas != nm1) {
			ldmas = *mlmas + *mumas + 1;
			if (jband) {
				ijob = 4;
			}
			else {
				ijob = 3;
			}
		}
		else {
			ldmas = nm1;
			ijob = 5;
		}
		/* ------ BANDWITH OF "MAS" NOT LARGER THAN BANDWITH OF "JAC" */
		if (*mlmas > *mljac || *mumas > *mujac) {
			s_wsle(&io___308);
			do_lio(&c__9, &c__1, "BANDWITH OF \"MAS\" NOT LARGER THAN BANDWI"
				"TH OF \"JAC\"", (ftnlen)51);
			e_wsle();
			arret = TRUE_;
		}
	}
	else {
		ldmas = 0;
		if (jband) {
			ijob = 2;
		}
		else {
			ijob = 1;
		}
	}
	ldmas2 = MAX(1, ldmas);
	/* ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK ----- */
	ieynew = 21;
	iedy1 = ieynew + *n;
	iedy = iedy1 + *n;
	ieak1 = iedy + *n;
	ieak2 = ieak1 + *n;
	ieak3 = ieak2 + *n;
	ieak4 = ieak3 + *n;
	ieak5 = ieak4 + *n;
	ieak6 = ieak5 + *n;
	iefx = ieak6 + *n;
	iecon = iefx + *n;
	iejac = iecon + (*n << 2);
	iemas = iejac + *n * ldjac;
	iee = iemas + nm1 * ldmas;

	/* ------ TOTAL STORAGE REQUIREMENT ----------- */
	istore = iee + nm1 * lde - 1;
	if (istore > *lwork) {
		s_wsle(&io___325);
		do_lio(&c__9, &c__1, " INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=", (
			ftnlen)43);
		do_lio(&c__3, &c__1, (char *)&istore, (ftnlen)sizeof(typeLInt));
		e_wsle();
		arret = TRUE_;
	}
	/* ------- ENTRY POINTS FOR INTEGER WORKSPACE ----- */
	ieip = 21;
	istore = ieip + nm1 - 1;
	if (istore > *liwork) {
		s_wsle(&io___327);
		do_lio(&c__9, &c__1, " INSUFF. STORAGE FOR IWORK, MIN. LIWORK=", (
			ftnlen)40);
		do_lio(&c__3, &c__1, (char *)&istore, (ftnlen)sizeof(typeLInt));
		e_wsle();
		arret = TRUE_;
	}
	/* ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1 */
	if (arret) {
		*idid = -1;
		return 0;
	}


	/* -------- CALL TO CORE INTEGRATOR ------------ */
	roscor_(n, (U_fp)fcn, x, &y[1], xend, &hmax, h__, &rtol[1], &atol[1],
		itol, (U_fp)jac, ijac, mljac, mujac, (U_fp)dfx, idfx, (U_fp)mas,
		mlmas, mumas, (U_fp)solout, iout, idid, &nmax, &uround, &meth, &
		ijob, &fac1, &fac2, &safe, &autnms, &implct, &jband, &pred, &
		ldjac, &lde, &ldmas2, &work[ieynew], &work[iedy1], &work[iedy], &
		work[ieak1], &work[ieak2], &work[ieak3], &work[ieak4], &work[
			ieak5], &work[ieak6], &work[iefx], &work[iejac], &work[iee], &
				work[iemas], &iwork[ieip], &work[iecon], &m1, &m2, &nm1, &nfcn, &
				njac, &nstep, &naccpt, &nrejct, &ndec, &nsol, tvec, xvec, uvec, pvec, dcdxvec, grampc, pfct);

	iwork[14] = nfcn;
	iwork[15] = njac;
	iwork[16] = nstep;
	iwork[17] = naccpt;
	iwork[18] = nrejct;
	iwork[19] = ndec;
	iwork[20] = nsol;
	/* ----------- RETURN ----------- */
	return 0;
} /* rodas_ */




/*  ----- ... AND HERE IS THE CORE INTEGRATOR  ---------- */

/* Subroutine */ int roscor_(typeLInt *n, S_fp fcn, typeRNum *x, typeRNum *
	y, typeRNum *xend, typeRNum *hmax, typeRNum *h__, typeRNum *
	rtol, typeRNum *atol, typeLInt *itol, S_fp jac, typeLInt *ijac,
	typeLInt *mljac, typeLInt *mujac, S_fp dfx, typeLInt *idfx, S_fp mas,
	typeLInt *mlmas, typeLInt *mumas, S_fp solout, typeLInt *iout, typeLInt *
	idid, typeLInt *nmax, typeRNum *uround, typeLInt *meth, typeLInt *ijob,
	typeRNum *fac1, typeRNum *fac2, typeRNum *safe, typeLogical *
	autnms, typeLogical *implct, typeLogical *banded, typeLogical *pred, typeLInt *
	ldjac, typeLInt *lde, typeLInt *ldmas, typeRNum *ynew, typeRNum *
	dy1, typeRNum *dy, typeRNum *ak1, typeRNum *ak2, typeRNum *
	ak3, typeRNum *ak4, typeRNum *ak5, typeRNum *ak6, typeRNum *
	fx, typeRNum *fjac, typeRNum *e, typeRNum *fmas, typeLInt *ip,
	typeRNum *cont, typeLInt *m1, typeLInt *m2, typeLInt *nm1, typeLInt *
	nfcn, typeLInt *njac, typeLInt *nstep, typeLInt *naccpt, typeLInt *nrejct,
	typeLInt *ndec, typeLInt *nsol, ctypeRNum *tvec, ctypeRNum *xvec, ctypeRNum *uvec, ctypeRNum *pvec, ctypeRNum *dcdxvec, const typeGRAMPC *grampc, const typeffctPtr pfct)
{
	/* Format strings */
	static char fmt_979[] = "(\002 EXIT OF RODAS AT X=\002,e18.4)";

	/* System generated locals */
	typeLInt fjac_dim1, fjac_offset, e_dim1, e_offset, fmas_dim1, fmas_offset,
		i__1, i__2, i__3, i__4;
	typeRNum d__1, d__2, d__3, d__4;

	/* Builtin functions */
	typeRNum d_sign(typeRNum *, typeRNum *), pow_dd(
		typeRNum *, typeRNum *);
	typeLInt s_wsfe(cilist *), do_fio(typeLInt *, char *, ftnlen), e_wsfe(void),
		s_wsle(cilist *), do_lio(typeLInt *, typeLInt *, char *, ftnlen),
		e_wsle(void);

	/* Local variables */
	static typeLInt i__, j, k, l;
	static typeRNum c2, c3, c4, d1, d2, d3, d4;
	static typeLInt j1;
	static typeRNum a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, c21,
		c31, c32, c41, c42, c43, c51, c52, c53, c54, c61, c62, c63, c64,
		c65, d21, d22, d23, d24, d25, d31, d32, d33, d34, d35;
	static typeLInt md, mm;
	static typeRNum sk, hd1, hd2, hd3, hd4;
	static typeLInt nn2, nn3;
	static typeRNum fac, hc21, hc31, hc32, hc41, hc42, hc43, hc51, hc52,
		hc53, hc54, hc61, hc62;
	static typeLInt ier, lrc;
	static typeRNum hc63, hc64, hc65, err, hacc;
	static typeLInt lbeg, lend;
	static typeRNum delt, hnew;
	static typeLogical last;
	static typeRNum hopt, gamma;
	extern /* Subroutine */ int rocoe_(typeLInt *, typeRNum *, typeRNum *,
		typeRNum *, typeRNum *, typeRNum *, typeRNum *,
		typeRNum *, typeRNum *, typeRNum *, typeRNum *,
		typeRNum *, typeRNum *, typeRNum *, typeRNum *,
		typeRNum *, typeRNum *, typeRNum *, typeRNum *,
		typeRNum *, typeRNum *, typeRNum *, typeRNum *,
		typeRNum *, typeRNum *, typeRNum *, typeRNum *,
		typeRNum *, typeRNum *, typeRNum *, typeRNum *,
		typeRNum *, typeRNum *, typeRNum *, typeRNum *,
		typeRNum *, typeRNum *, typeRNum *, typeRNum *,
		typeRNum *, typeRNum *, typeRNum *, typeRNum *,
		typeRNum *);
	static typeRNum ysafe, hmaxn;
	static typeLInt nsing;
	static typeRNum xdelt;
	static typeLInt irtrn;
	static typeRNum erracc;
	static typeLInt mujacj;
	extern /* Subroutine */ int decomr_(typeLInt *, typeRNum *, typeLInt *,
		typeRNum *, typeLInt *, typeLInt *, typeLInt *, typeLInt *, typeLInt
		*, typeLInt *, typeRNum *, typeRNum *, typeLInt *, typeLInt *,
		typeLInt *, typeLInt *, typeLogical *, typeLInt *);
	static typeRNum facgus;
	static typeLogical reject;
	static typeLInt mujacp;
	static typeRNum posneg;
	extern /* Subroutine */ int slvrod_(typeLInt *, typeRNum *, typeLInt *,
		typeLInt *, typeLInt *, typeRNum *, typeLInt *, typeLInt *, typeLInt
		*, typeLInt *, typeLInt *, typeLInt *, typeRNum *, typeRNum *,
		typeLInt *, typeLInt *, typeRNum *, typeRNum *, typeRNum *,
		typeRNum *, typeRNum *, typeLInt *, typeLogical *);

	/* Fortran I/O blocks */
	static cilist io___422 = { 0, 6, 0, fmt_979, 0 };
	static cilist io___423 = { 0, 6, 0, 0, 0 };
	static cilist io___424 = { 0, 6, 0, fmt_979, 0 };
	static cilist io___425 = { 0, 6, 0, 0, 0 };
	static cilist io___426 = { 0, 6, 0, fmt_979, 0 };
	static cilist io___427 = { 0, 6, 0, 0, 0 };
	static cilist io___428 = { 0, 6, 0, fmt_979, 0 };

	/* ---------------------------------------------------------- */
	/*     CORE INTEGRATOR FOR RODAS */
	/*     PARAMETERS SAME AS IN RODAS WITH WORKSPACE ADDED */
	/* ---------------------------------------------------------- */
	/*         DECLARATIONS */
	/* ---------------------------------------------------------- */
	/* *** *** *** *** *** *** *** */
	/*  INITIALISATIONS */
	/* *** *** *** *** *** *** *** */
			/* Parameter adjustments */
	--cont;
	--fx;
	--ak6;
	--ak5;
	--ak4;
	--ak3;
	--ak2;
	--ak1;
	--dy;
	--dy1;
	--ynew;
	--y;
	--rtol;
	--atol;
	fjac_dim1 = *ldjac;
	fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;
	--ip;
	fmas_dim1 = *ldmas;
	fmas_offset = 1 + fmas_dim1;
	fmas -= fmas_offset;
	e_dim1 = *lde;
	e_offset = 1 + e_dim1;
	e -= e_offset;

	/* Function Body */
	conros_1.nn = *n;
	nn2 = *n << 1;
	nn3 = *n * 3;
	lrc = *n << 2;
	/* ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ---------- */
	if (*implct) {
		(*mas)(nm1, &fmas[fmas_offset], ldmas, grampc, pfct);
	}
	/* ------ SET THE PARAMETERS OF THE METHOD ----- */
	rocoe_(meth, &a21, &a31, &a32, &a41, &a42, &a43, &a51, &a52, &a53, &a54, &
		c21, &c31, &c32, &c41, &c42, &c43, &c51, &c52, &c53, &c54, &c61, &
		c62, &c63, &c64, &c65, &gamma, &c2, &c3, &c4, &d1, &d2, &d3, &d4,
		&d21, &d22, &d23, &d24, &d25, &d31, &d32, &d33, &d34, &d35);
	/* --- INITIAL PREPARATIONS */
	if (*m1 > 0) {
		*ijob += 10;
	}
	d__1 = *xend - *x;
	posneg = d_sign(&c_b361, &d__1);
	/* Computing MIN */
	d__2 = DABS(*hmax), d__3 = (d__1 = *xend - *x, DABS(d__1));
	hmaxn = MIN(d__2, d__3);
	if (DABS(*h__) <= *uround * ((typeRNum)10)) {
		*h__ = (typeRNum)1e-6;
	}
	/* Computing MIN */
	d__1 = DABS(*h__);
	*h__ = MIN(d__1, hmaxn);
	*h__ = d_sign(h__, &posneg);
	reject = FALSE_;
	last = FALSE_;
	nsing = 0;
	irtrn = 1;
	if (*autnms) {
		hd1 = 0;
		hd2 = 0;
		hd3 = 0;
		hd4 = 0;
	}
	/* -------- PREPARE BAND-WIDTHS -------- */
	linal_1.mbdiag = *mumas + 1;
	if (*banded) {
		linal_1.mle = *mljac;
		linal_1.mue = *mujac;
		linal_1.mbjac = *mljac + *mujac + 1;
		linal_1.mbb = *mlmas + *mumas + 1;
		linal_1.mdiag = linal_1.mle + linal_1.mue + 1;
		linal_1.mdiff = linal_1.mle + linal_1.mue - *mumas;
	}
	if (*iout != 0) {
		conros_1.xold = *x;
		irtrn = 1;
		conros_1.hout = *h__;
		i__1 = *naccpt + 1;
		(*solout)(&i__1, &conros_1.xold, x, &conros_1.hout, &y[1], &cont[1],
			&lrc, n, grampc, pfct, &irtrn);
		if (irtrn < 0) {
			goto L179;
		}
	}
	/* --- BASIC INTEGRATION STEP */
L1:
	if (*nstep > *nmax) {
		goto L178;
	}
	if (DABS(*h__) * .1 <= DABS(*x) * *uround) {
		goto L177;
	}
	if (last) {
		*h__ = hopt;
		*idid = 1;
		return 0;
	}
	hopt = *h__;
	if ((*x + *h__ * 1.0001 - *xend) * posneg >= 0) {
		*h__ = *xend - *x;
		last = TRUE_;
	}
	/* *** *** *** *** *** *** *** */
	/*  COMPUTATION OF THE JACOBIAN */
	/* *** *** *** *** *** *** *** */
	(*fcn)(n, x, &y[1], &dy1[1], tvec, xvec, uvec, pvec, dcdxvec, grampc, pfct);
	++(*nfcn);
	++(*njac);
	if (*ijac == 0) {
		/* --- COMPUTE JACOBIAN MATRIX NUMERICALLY */
		if (*banded) {
			/* --- JACOBIAN IS BANDED */
			mujacp = *mujac + 1;
			md = MIN(linal_1.mbjac, *n);
			i__1 = *m1 / *m2 + 1;
			for (mm = 1; mm <= i__1; ++mm) {
				i__2 = md;
				for (k = 1; k <= i__2; ++k) {
					j = k + (mm - 1) * *m2;
				L12:
					ak2[j] = y[j];
					/* Computing MAX */
					d__2 = (typeRNum)1e-5, d__3 = (d__1 = y[j], DABS(d__1));
					ak3[j] = SQRT(*uround * MAX(d__2, d__3));
					y[j] += ak3[j];
					j += md;
					if (j <= mm * *m2) {
						goto L12;
					}
					(*fcn)(n, x, &y[1], &ak1[1], tvec, xvec, uvec, pvec, dcdxvec, grampc, pfct);
					j = k + (mm - 1) * *m2;
					j1 = k;
					/* Computing MAX */
					i__3 = 1, i__4 = j1 - *mujac;
					lbeg = MAX(i__3, i__4) + *m1;
				L14:
					/* Computing MIN */
					i__3 = *m2, i__4 = j1 + *mljac;
					lend = MIN(i__3, i__4) + *m1;
					y[j] = ak2[j];
					mujacj = mujacp - j1 - *m1;
					i__3 = lend;
					for (l = lbeg; l <= i__3; ++l) {
						fjac[l + mujacj + j * fjac_dim1] = (ak1[l] - dy1[l]) /
							ak3[j];
					}
					j += md;
					j1 += md;
					lbeg = lend + 1;
					if (j <= mm * *m2) {
						goto L14;
					}
				}
			}
		}
		else {
			/* --- JACOBIAN IS FULL */
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
				ysafe = y[i__];
				/* Computing MAX */
				d__1 = (typeRNum)1e-5, d__2 = DABS(ysafe);
				delt = SQRT(*uround * MAX(d__1, d__2));
				y[i__] = ysafe + delt;
				(*fcn)(n, x, &y[1], &ak1[1], tvec, xvec, uvec, pvec, dcdxvec, grampc, pfct);
				i__2 = *n;
				for (j = *m1 + 1; j <= i__2; ++j) {
					fjac[j - *m1 + i__ * fjac_dim1] = (ak1[j] - dy1[j]) /
						delt;
				}
				y[i__] = ysafe;

			}
		}
	}
	else {
		/* --- COMPUTE JACOBIAN MATRIX ANALYTICALLY */
		(*jac)(n, x, &y[1], &fjac[fjac_offset], ldjac, tvec, xvec, uvec, pvec, dcdxvec, grampc, pfct);
	}
	if (!(*autnms)) {
		if (*idfx == 0) {
			/* --- COMPUTE NUMERICALLY THE DERIVATIVE WITH RESPECT TO X */
			/* Computing MAX */
			d__1 = (typeRNum)1e-5, d__2 = DABS(*x);
			delt = SQRT(*uround * MAX(d__1, d__2));
			xdelt = *x + delt;
			(*fcn)(n, &xdelt, &y[1], &ak1[1], tvec, xvec, uvec, pvec, dcdxvec, grampc, pfct);
			i__1 = *n;
			for (j = 1; j <= i__1; ++j) {
				fx[j] = (ak1[j] - dy1[j]) / delt;
			}
		}
		else {
			/* --- COMPUTE ANALYTICALLY THE DERIVATIVE WITH RESPECT TO X */
			(*dfx)(n, x, &y[1], &fx[1], tvec, xvec, uvec, pvec, dcdxvec, grampc, pfct);
		}
	}
L2:
	/* *** *** *** *** *** *** *** */
	/*  COMPUTE THE STAGES */
	/* *** *** *** *** *** *** *** */
	fac = (typeRNum) 1. / (*h__ * gamma);
	decomr_(n, &fjac[fjac_offset], ldjac, &fmas[fmas_offset], ldmas, mlmas,
		mumas, m1, m2, nm1, &fac, &e[e_offset], lde, &ip[1], &ier, ijob,
		implct, &ip[1]);
	if (ier != 0) {
		goto L80;
	}
	++(*ndec);
	/* --- PREPARE FOR THE COMPUTATION OF THE 6 STAGES */
	hc21 = c21 / *h__;
	hc31 = c31 / *h__;
	hc32 = c32 / *h__;
	hc41 = c41 / *h__;
	hc42 = c42 / *h__;
	hc43 = c43 / *h__;
	hc51 = c51 / *h__;
	hc52 = c52 / *h__;
	hc53 = c53 / *h__;
	hc54 = c54 / *h__;
	hc61 = c61 / *h__;
	hc62 = c62 / *h__;
	hc63 = c63 / *h__;
	hc64 = c64 / *h__;
	hc65 = c65 / *h__;
	if (!(*autnms)) {
		hd1 = *h__ * d1;
		hd2 = *h__ * d2;
		hd3 = *h__ * d3;
		hd4 = *h__ * d4;
	}
	/* --- THE STAGES */
	slvrod_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
		ldmas, mlmas, mumas, m1, m2, nm1, &fac, &e[e_offset], lde, &ip[1],
		&dy1[1], &ak1[1], &fx[1], &ynew[1], &hd1, ijob, &c_false);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		ynew[i__] = y[i__] + a21 * ak1[i__];
	}
	d__1 = *x + c2 * *h__;
	(*fcn)(n, &d__1, &ynew[1], &dy[1], tvec, xvec, uvec, pvec, dcdxvec, grampc, pfct);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		ynew[i__] = hc21 * ak1[i__];
	}
	slvrod_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
		ldmas, mlmas, mumas, m1, m2, nm1, &fac, &e[e_offset], lde, &ip[1],
		&dy[1], &ak2[1], &fx[1], &ynew[1], &hd2, ijob, &c_true);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		ynew[i__] = y[i__] + a31 * ak1[i__] + a32 * ak2[i__];
	}
	d__1 = *x + c3 * *h__;
	(*fcn)(n, &d__1, &ynew[1], &dy[1], tvec, xvec, uvec, pvec, dcdxvec, grampc, pfct);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		ynew[i__] = hc31 * ak1[i__] + hc32 * ak2[i__];
	}
	slvrod_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
		ldmas, mlmas, mumas, m1, m2, nm1, &fac, &e[e_offset], lde, &ip[1],
		&dy[1], &ak3[1], &fx[1], &ynew[1], &hd3, ijob, &c_true);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		ynew[i__] = y[i__] + a41 * ak1[i__] + a42 * ak2[i__] + a43 * ak3[i__];
	}
	d__1 = *x + c4 * *h__;
	(*fcn)(n, &d__1, &ynew[1], &dy[1], tvec, xvec, uvec, pvec, dcdxvec, grampc, pfct);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		ynew[i__] = hc41 * ak1[i__] + hc42 * ak2[i__] + hc43 * ak3[i__];
	}
	slvrod_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
		ldmas, mlmas, mumas, m1, m2, nm1, &fac, &e[e_offset], lde, &ip[1],
		&dy[1], &ak4[1], &fx[1], &ynew[1], &hd4, ijob, &c_true);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		ynew[i__] = y[i__] + a51 * ak1[i__] + a52 * ak2[i__] + a53 * ak3[i__]
			+ a54 * ak4[i__];
	}
	d__1 = *x + *h__;
	(*fcn)(n, &d__1, &ynew[1], &dy[1], tvec, xvec, uvec, pvec, dcdxvec, grampc, pfct);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		ak6[i__] = hc52 * ak2[i__] + hc54 * ak4[i__] + hc51 * ak1[i__] + hc53
			* ak3[i__];
	}
	slvrod_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
		ldmas, mlmas, mumas, m1, m2, nm1, &fac, &e[e_offset], lde, &ip[1],
		&dy[1], &ak5[1], &fx[1], &ak6[1], &c_b374, ijob, &c_true);
	/* ------------ EMBEDDED SOLUTION --------------- */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		ynew[i__] += ak5[i__];
	}
	d__1 = *x + *h__;
	(*fcn)(n, &d__1, &ynew[1], &dy[1], tvec, xvec, uvec, pvec, dcdxvec, grampc, pfct);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		cont[i__] = hc61 * ak1[i__] + hc62 * ak2[i__] + hc65 * ak5[i__] +
			hc64 * ak4[i__] + hc63 * ak3[i__];
	}
	slvrod_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
		ldmas, mlmas, mumas, m1, m2, nm1, &fac, &e[e_offset], lde, &ip[1],
		&dy[1], &ak6[1], &fx[1], &cont[1], &c_b374, ijob, &c_true);
	/* ------------ NEW SOLUTION --------------- */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		ynew[i__] += ak6[i__];
	}
	*nsol += 6;
	*nfcn += 5;
	/* ------------ DENSE OUTPUT ---------- */
	if (*iout != 0) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			cont[i__] = y[i__];
			cont[i__ + nn2] = d21 * ak1[i__] + d22 * ak2[i__] + d23 * ak3[i__]
				+ d24 * ak4[i__] + d25 * ak5[i__];
			cont[i__ + nn3] = d31 * ak1[i__] + d32 * ak2[i__] + d33 * ak3[i__]
				+ d34 * ak4[i__] + d35 * ak5[i__];
		}
	}
	/* *** *** *** *** *** *** *** */
	/*  ERROR ESTIMATION */
	/* *** *** *** *** *** *** *** */
	++(*nstep);
	/* ------------ COMPUTE ERROR ESTIMATION ---------------- */
	err = 0;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		if (*itol == 0) {
			/* Computing MAX */
			d__3 = (d__1 = y[i__], DABS(d__1)), d__4 = (d__2 = ynew[i__], DABS(
				d__2));
			sk = atol[1] + rtol[1] * MAX(d__3, d__4);
		}
		else {
			/* Computing MAX */
			d__3 = (d__1 = y[i__], DABS(d__1)), d__4 = (d__2 = ynew[i__], DABS(
				d__2));
			sk = atol[i__] + rtol[i__] * MAX(d__3, d__4);
		}
		/* Computing 2nd power */
		d__1 = ak6[i__] / sk;
		err += d__1 * d__1;
	}
	err = SQRT(err / *n);
	/* --- COMPUTATION OF HNEW */
	/* --- WE REQUIRE .2<=HNEW/H<=6. */
	/* Computing MAX */
	/* Computing MIN */
	d__3 = *fac1, d__4 = pow_dd(&err, &c_b378) / *safe;
	d__1 = *fac2, d__2 = MIN(d__3, d__4);
	fac = MAX(d__1, d__2);
	hnew = *h__ / fac;
	/* *** *** *** *** *** *** *** */
	/*  IS THE ERROR SMALL ENOUGH ? */
	/* *** *** *** *** *** *** *** */
	if (err <= (typeRNum)1) {
		/* --- STEP IS ACCEPTED */
		++(*naccpt);
		if (*pred) {
			/*       --- PREDICTIVE CONTROLLER OF GUSTAFSSON */
			if (*naccpt > 1) {
				/* Computing 2nd power */
				d__2 = err;
				d__1 = d__2 * d__2 / erracc;
				facgus = hacc / *h__ * pow_dd(&d__1, &c_b378) / *safe;
				/* Computing MAX */
				d__1 = *fac2, d__2 = MIN(*fac1, facgus);
				facgus = MAX(d__1, d__2);
				fac = MAX(fac, facgus);
				hnew = *h__ / fac;
			}
			hacc = *h__;
			erracc = MAX((typeRNum).01, err);
		}
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
			y[i__] = ynew[i__];
		}
		conros_1.xold = *x;
		*x += *h__;
		if (*iout != 0) {
			i__1 = *n;
			for (i__ = 1; i__ <= i__1; ++i__) {
				cont[conros_1.nn + i__] = y[i__];
			}
			irtrn = 1;
			conros_1.hout = *h__;
			i__1 = *naccpt + 1;
			(*solout)(&i__1, &conros_1.xold, x, &conros_1.hout, &y[1],
				&cont[1], &lrc, n, grampc, pfct, &irtrn);
			if (irtrn < 0) {
				goto L179;
			}
		}
		if (DABS(hnew) > hmaxn) {
			hnew = posneg * hmaxn;
		}
		if (reject) {
			/* Computing MIN */
			d__1 = DABS(hnew), d__2 = DABS(*h__);
			hnew = posneg * MIN(d__1, d__2);
		}
		reject = FALSE_;
		*h__ = hnew;
		goto L1;
	}
	else {
		/* --- STEP IS REJECTED */
		reject = TRUE_;
		last = FALSE_;
		*h__ = hnew;
		if (*naccpt >= 1) {
			++(*nrejct);
		}
		goto L2;
	}
	/* --- SINGULAR MATRIX */
L80:
	++nsing;
	if (nsing >= 5) {
		goto L176;
	}
	*h__ *= .5;
	reject = TRUE_;
	last = FALSE_;
	goto L2;
	/* --- FAIL EXIT */
L176:
	s_wsfe(&io___422);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(typeRNum));
	e_wsfe();
	s_wsle(&io___423);
	do_lio(&c__9, &c__1, " MATRIX IS REPEATEDLY SINGULAR, IER=", (ftnlen)36);
	do_lio(&c__3, &c__1, (char *)&ier, (ftnlen)sizeof(typeLInt));
	e_wsle();
	*idid = -4;
	return 0;
L177:
	s_wsfe(&io___424);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(typeRNum));
	e_wsfe();
	s_wsle(&io___425);
	do_lio(&c__9, &c__1, " STEP SIZE T0O SMALL, H=", (ftnlen)24);
	do_lio(&c__5, &c__1, (char *)&(*h__), (ftnlen)sizeof(typeRNum));
	e_wsle();
	*idid = -3;
	return 0;
L178:
	s_wsfe(&io___426);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(typeRNum));
	e_wsfe();
	s_wsle(&io___427);
	do_lio(&c__9, &c__1, " MORE THAN NMAX =", (ftnlen)17);
	do_lio(&c__3, &c__1, (char *)&(*nmax), (ftnlen)sizeof(typeLInt));
	do_lio(&c__9, &c__1, "STEPS ARE NEEDED", (ftnlen)16);
	e_wsle();
	*idid = -2;
	return 0;
	/* --- EXIT CAUSED BY SOLOUT */
L179:
	s_wsfe(&io___428);
	do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(typeRNum));
	e_wsfe();
	*idid = 2;
	return 0;
} /* roscor_ */


typeRNum contro(int i, int *N, typeRNum x, typeRNum *xold, typeRNum *h,
	typeRNum *cont, int *lrc)
{
	/* ---------------------------------------------------------- */
	/*     THIS FUNCTION CAN BE USED FOR CONTINUOUS OUTPUT IN CONNECTION */
	/*     WITH THE OUTPUT-SUBROUTINE FOR RODAS. IT PROVIDES AN */
	/*     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION AT X. */
	/* ---------------------------------------------------------- */
			/* System generated locals */
	typeRNum ret_val;

	/* Local variables */
	static typeRNum s;

	/* Function Body */
	s = (x - xold[0]) / h[0];

	ret_val = cont[i] * (1 - s) + s * (cont[i + N[0]] + (1 - s) *
		(cont[i + N[0] * 2] + s * cont[i + N[0] * 3]));
	return ret_val;
} /* contro_ */


/* Subroutine */ int rocoe_(typeLInt *meth, typeRNum *a21, typeRNum *a31,
	typeRNum *a32, typeRNum *a41, typeRNum *a42, typeRNum *a43,
	typeRNum *a51, typeRNum *a52, typeRNum *a53, typeRNum *a54,
	typeRNum *c21, typeRNum *c31, typeRNum *c32, typeRNum *c41,
	typeRNum *c42, typeRNum *c43, typeRNum *c51, typeRNum *c52,
	typeRNum *c53, typeRNum *c54, typeRNum *c61, typeRNum *c62,
	typeRNum *c63, typeRNum *c64, typeRNum *c65, typeRNum *gamma,
	typeRNum *c2, typeRNum *c3, typeRNum *c4, typeRNum *d1,
	typeRNum *d2, typeRNum *d3, typeRNum *d4, typeRNum *d21,
	typeRNum *d22, typeRNum *d23, typeRNum *d24, typeRNum *d25,
	typeRNum *d31, typeRNum *d32, typeRNum *d33, typeRNum *d34,
	typeRNum *d35)
{
	switch (*meth) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	}
L1:
	*c2 = (typeRNum) .386;
	*c3 = (typeRNum) .21;
	*c4 = (typeRNum) .63;
	*d1 = (typeRNum) .25;
	*d2 = (typeRNum)-.1043;
	*d3 = (typeRNum) .1035;
	*d4 = (typeRNum)-.03620000000000023;
	*a21 = (typeRNum) 1.544;
	*a31 = (typeRNum) .9466785280815826;
	*a32 = (typeRNum) .2557011698983284;
	*a41 = (typeRNum) 3.314825187068521;
	*a42 = (typeRNum) 2.896124015972201;
	*a43 = (typeRNum) .9986419139977817;
	*a51 = (typeRNum) 1.221224509226641;
	*a52 = (typeRNum) 6.019134481288629;
	*a53 = (typeRNum) 12.53708332932087;
	*a54 = (typeRNum)-.687886036105895;
	*c21 = (typeRNum)-5.6688;
	*c31 = (typeRNum)-2.430093356833875;
	*c32 = (typeRNum)-.2063599157091915;
	*c41 = (typeRNum)-.1073529058151375;
	*c42 = (typeRNum)-9.594562251023355;
	*c43 = (typeRNum)-20.47028614809616;
	*c51 = (typeRNum) 7.496443313967647;
	*c52 = (typeRNum)-10.24680431464352;
	*c53 = (typeRNum)-33.99990352819905;
	*c54 = (typeRNum) 11.7089089320616;
	*c61 = (typeRNum) 8.083246795921522;
	*c62 = (typeRNum)-7.981132988064893;
	*c63 = (typeRNum)-31.52159432874371;
	*c64 = (typeRNum) 16.31930543123136;
	*c65 = (typeRNum)-6.058818238834054;
	*gamma = (typeRNum) .25;
	*d21 = (typeRNum) 10.12623508344586;
	*d22 = (typeRNum)-7.487995877610167;
	*d23 = (typeRNum)-34.80091861555747;
	*d24 = (typeRNum)-7.992771707568823;
	*d25 = (typeRNum) 1.025137723295662;
	*d31 = (typeRNum)-.6762803392801253;
	*d32 = (typeRNum) 6.087714651680015;
	*d33 = (typeRNum) 16.43084320892478;
	*d34 = (typeRNum) 24.76722511418386;
	*d35 = (typeRNum)-6.594389125716872;
	return 0;

L2:
	*c2 = (typeRNum) .3507221;
	*c3 = (typeRNum) .2557041;
	*c4 = (typeRNum) .681779;
	*d1 = (typeRNum) .25;
	*d2 = (typeRNum)-.06902209999999998;
	*d3 = (typeRNum)-9.671999999999459e-4;
	*d4 = (typeRNum)-.08797900000000025;
	*a21 = (typeRNum) 1.4028884;
	*a31 = (typeRNum) .6581212688557198;
	*a32 = (typeRNum)-1.320936088384301;
	*a41 = (typeRNum) 7.131197445744498;
	*a42 = (typeRNum) 16.02964143958207;
	*a43 = (typeRNum)-5.561572550509766;
	*a51 = (typeRNum) 22.73885722420363;
	*a52 = (typeRNum) 67.38147284535289;
	*a53 = (typeRNum)-31.2187749303856;
	*a54 = (typeRNum) .7285641833203814;
	*c21 = (typeRNum)-5.1043536;
	*c31 = (typeRNum)-2.899967805418783;
	*c32 = (typeRNum) 4.040399359702244;
	*c41 = (typeRNum)-32.64449927841361;
	*c42 = (typeRNum)-99.35311008728094;
	*c43 = (typeRNum) 49.99119122405989;
	*c51 = (typeRNum)-76.46023087151691;
	*c52 = (typeRNum)-278.5942120829058;
	*c53 = (typeRNum) 153.9294840910643;
	*c54 = (typeRNum) 10.97101866258358;
	*c61 = (typeRNum)-76.29701586804983;
	*c62 = (typeRNum)-294.2795630511232;
	*c63 = (typeRNum) 162.0029695867566;
	*c64 = (typeRNum) 23.6516690309527;
	*c65 = (typeRNum)-7.652977706771382;
	*gamma = (typeRNum) .25;
	*d21 = (typeRNum)-38.71940424117216;
	*d22 = (typeRNum)-135.8025833007622;
	*d23 = (typeRNum) 64.51068857505875;
	*d24 = (typeRNum)-4.192663174613162;
	*d25 = (typeRNum)-2.53193205033506;
	*d31 = (typeRNum)-14.99268484949843;
	*d32 = (typeRNum)-76.30242396627033;
	*d33 = (typeRNum) 58.65928432851416;
	*d34 = (typeRNum) 16.61359034616402;
	*d35 = (typeRNum)-.6758691794084156;
	return 0;

	/* Coefficients for RODAS with order 4 for linear parabolic problems */
	/* Gerd Steinebach (1993) */
L3:
	*gamma = (typeRNum) .25;
	*c2 = *gamma * ((typeRNum)3);
	*c3 = (typeRNum) .21;
	*c4 = (typeRNum) .63;
	*d1 = (typeRNum) .25;
	*d2 = (typeRNum)-.5;
	*d3 = (typeRNum)-.023504;
	*d4 = (typeRNum)-.0362;
	*a21 = (typeRNum)3;
	*a31 = (typeRNum) 1.831036793486759;
	*a32 = (typeRNum) .4955183967433795;
	*a41 = (typeRNum) 2.304376582692669;
	*a42 = (typeRNum)-.05249275245743001;
	*a43 = (typeRNum)-1.176798761832782;
	*a51 = (typeRNum)-7.170454962423024;
	*a52 = (typeRNum)-4.741636671481785;
	*a53 = (typeRNum)-16.31002631330971;
	*a54 = (typeRNum)-1.062004044111401;
	*c21 = (typeRNum)-12;
	*c31 = (typeRNum)-8.791795173947035;
	*c32 = (typeRNum)-2.207865586973518;
	*c41 = (typeRNum) 10.81793056857153;
	*c42 = (typeRNum) 6.780270611428266;
	*c43 = (typeRNum) 19.5348594464241;
	*c51 = (typeRNum) 34.19095006749676;
	*c52 = (typeRNum) 15.49671153725963;
	*c53 = (typeRNum) 54.7476087596413;
	*c54 = (typeRNum) 14.16005392148534;
	*c61 = (typeRNum) 34.62605830930532;
	*c62 = (typeRNum) 15.30084976114473;
	*c63 = (typeRNum) 56.99955578662667;
	*c64 = (typeRNum) 18.40807009793095;
	*c65 = (typeRNum)-5.714285714285717;

	*d21 = (typeRNum) 25.09876703708589;
	*d22 = (typeRNum) 11.62013104361867;
	*d23 = (typeRNum) 28.49148307714626;
	*d24 = (typeRNum)-5.664021568594133;
	*d25 = (typeRNum)0;
	*d31 = (typeRNum) 1.638054557396973;
	*d32 = (typeRNum)-.7373619806678748;
	*d33 = (typeRNum) 8.47791821923899;
	*d34 = (typeRNum) 15.9925314877952;
	*d35 = (typeRNum)-1.882352941176471;
	return 0;
} /* rocoe_ */

#endif
