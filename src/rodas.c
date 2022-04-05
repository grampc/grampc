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


#include "rodas.h"
#include "f2cmod.h"
#include "rodas_decsol_f2c.h"



void intsysRodas(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum *tvec, ctypeRNum *xvec, ctypeRNum *uvec, ctypeRNum *pvec, const typeGRAMPC *grampc, const typeffctPtr pfct)
{
	typeRNum *dcdxvec = grampc->rws->dcdx;

	int i, j;

	int IFCN = grampc->opt->FlagsRodas[0];  /* 0 --> right hand side independent of time t */
	int IDFX = grampc->opt->FlagsRodas[1];  /* 0 --> DF/DX is numerically computed */


	int IOUT = 1;  /* 1(0) --> subroutine 'solout' is called after each succesive step of integration (is never called) */
	int ITOL = 0;  /* rel. & abs. tol. are scalars */

	int IDID = 0;  /* output flag */
	               /*  1 --> computation successful */
	               /*  2 --> computation successful (interrupted by solout) */
	               /* -1 --> input is not consistent */
	               /* -2 --> larger nmax is needed */
	               /* -3 --> step size becomes too small */
	               /* -4 --> matrix is repeatedly singular */

	typeRNum t0, t1;                                /* initial and end time */
	typeRNum tolr = grampc->opt->IntegratorRelTol;	/* relative tolerance */
	typeRNum tola = grampc->opt->IntegratorAbsTol;	/* absolut tolerance */
	typeRNum HINIT = (typeRNum)1e-2;                /* initial stepsize */

	int IJAC;  /* 1(0) -> analytical (numerical) jacobian (partial derivatives of right hand side w.r.t. state) */
	int MLJAC; /* no. of lower diagonals of jacobian */
	int MUJAC; /* no. of upper diagonals of jacobian */

	int IMAS;  /* 1 --> mass matrix is supplied */
	int MLMAS; /* no. of lower diagonals of mass matrix */
	int MUMAS; /* no. of upper diagonals of mass matrix */

	grampc->rws->iparRodas[0] = pInt;
	grampc->rws->iparRodas[1] = Nint;
	grampc->rws->iparRodas[2] = 0;

	grampc->rws->iworkRodas[0] = grampc->opt->IntegratorMaxSteps; /* MAXIMAL NUMBER OF ALLOWED STEPS (default is 100.000) */
	grampc->rws->iworkRodas[1] = 2;  /* SWITCH FOR THE CHOICE OF THE COEFFICIENTS (default: 0/1) */
	grampc->rws->iworkRodas[2] = 0;  /* SWITCH FOR STEP SIZE STRATEGY (default: 0/1) */

	/* grampc->rws->workRodas[0] = 1e-16;*/	/* UROUND, THE ROUNDING UNIT (default: e-16) */
	/* grampc->rws->workRodas[1] = 1e-3;*/	/* MAXIMAL STEP SIZE (default: TEND-T) */

	/* grampc->rws->workRodas[2] = 0.2;	*/ /* workRodas[2] & workRodas[2] --> PARAMETERs FOR STEP SIZE SELECTION (default: 0,2 & 6.0) */
	/* grampc->rws->workRodas[3] = 6.0;	*/ /* the new step size is chosen subject to the restriction workRodas[2] <= HNEW/HOLD <= workRodas[3] */



	if (pInt == FWINT) {
		t0 = tvec[0];
		t1 = t0 + grampc->param->Thor / (typeRNum)(grampc->opt->Nhor - 1)*(Nint - 1);

		IJAC = grampc->opt->FlagsRodas[2];   /* 1(0) -> analytical (numerical) jacobian (partial derivatives of right hand side w.r.t. state) */
		MLJAC = grampc->opt->FlagsRodas[4];  /* no. of lower diagonals of jacobian */
		MUJAC = grampc->opt->FlagsRodas[5];  /* no. of upper diagonals of jacobian */

		IMAS = grampc->opt->FlagsRodas[3];   /* 1 --> mass matrix is supplied */
		MLMAS = grampc->opt->FlagsRodas[6];  /* no. of lower diagonals of mass matrix */
		MUMAS = grampc->opt->FlagsRodas[7];  /* no. of upper diagonals of mass matrix */
	}
	else {
		t0 = tvec[-(Nint - 1)];
		t1 = t0 + grampc->param->Thor / (typeRNum)(grampc->opt->Nhor - 1)*(Nint - 1);
		/* calling RODAS ... */
		tvec = tvec - (Nint - 1);
		xvec = xvec - (Nint - 1)*grampc->param->Nx;
		uvec = uvec - (Nint - 1)*grampc->param->Nu;

		/*jacobian of right hand side of adjoint system equal to transpose(dfdx)*/
		IJAC = grampc->opt->FlagsRodas[2];   /* 1(0) -> analytical (numerical) jacobian (partial derivatives of right hand side w.r.t. state) */
		MLJAC = grampc->opt->FlagsRodas[5];  /* no. of lower diagonals of jacobian */
		MUJAC = grampc->opt->FlagsRodas[4];  /* no. of upper diagonals of jacobian */

		IMAS = grampc->opt->FlagsRodas[3];   /* 1 --> mass matrix is supplied */
		MLMAS = grampc->opt->FlagsRodas[7];  /* no. of lower diagonals of mass matrix */
		MUMAS = grampc->opt->FlagsRodas[6];  /* no. of upper diagonals of mass matrix */
	}

	/* calling RODAS ... */
	rodas_(&grampc->param->Nx,
		(U_fp)ffctRodas, &IFCN,
		&t0, y, &t1, &HINIT, &tolr, &tola, &ITOL,
		(U_fp)dfdxRodas, &IJAC, &MLJAC, &MUJAC,
		(U_fp)dfdtRodas, &IDFX,
		(U_fp)MfctRodas, &IMAS, &MLMAS, &MUMAS,
		(U_fp)solout, &IOUT,
		&grampc->rws->workRodas[0], &grampc->rws->lworkRodas, &grampc->rws->iworkRodas[0], &grampc->rws->liworkRodas,
		tvec, xvec, uvec, pvec, dcdxvec, grampc, pfct, &IDID);

	if (IDID == -1) { grampc->sol->status |= STATUS_INTEGRATOR_INPUT_NOT_CONSISTENT; }
	if (IDID == -2) { grampc->sol->status |= STATUS_INTEGRATOR_MAXSTEPS; }
	if (IDID == -3) { grampc->sol->status |= STATUS_INTEGRATOR_STEPS_TOO_SMALL; }
	if (IDID == -4) { grampc->sol->status |= STATUS_INTEGRATOR_MATRIX_IS_SINGULAR; }

	if (pInt == FWINT) {
		for (i = 0; i < grampc->param->Nx*Nint; i++) { y[i] = grampc->rws->rparRodas[i]; }
	}
	else
	{
		/* solution (REVERSE) */
		y = y - (Nint - 1)*grampc->param->Nx;
		for (i = 0; i < Nint; i++) {
			for (j = 0; j < grampc->param->Nx; j++) { y[i*grampc->param->Nx + j] = grampc->rws->rparRodas[((Nint - 1) - i)*grampc->param->Nx + j]; }
		}
	}



}

void ffctRodas(int *N, typeRNum *t, typeRNum *x, typeRNum *rhs, ctypeRNum *tvec, ctypeRNum *xvec, ctypeRNum *uvec, ctypeRNum *pvec, ctypeRNum *dcdxvec, const typeGRAMPC *grampc, const typeffctPtr pfct)
{
	typeRNum *uakt = grampc->rws->rwsGeneral + LWadjsys; /* size NU */
	typeRNum *xakt = uakt + grampc->param->Nu;           /* size NX */
	typeRNum *dcdxakt = xakt + grampc->param->Nx;        /* size NX */

	typeInt i;
	typeRNum takt;
	ctypeInt pInt = grampc->rws->iparRodas[0];
	ctypeInt Nres = grampc->rws->iparRodas[1];

	if (pInt == BWINT) {
		takt = grampc->param->Thor - t[0];
	}
	else {
		takt = t[0];
	}

	interplin(uakt, tvec, uvec, takt, grampc->param->Nu, Nres, 1);

	if (pInt == BWINT) {
		interplin(xakt, tvec, xvec, takt, grampc->param->Nx, Nres, 1);
		interplin(dcdxakt, tvec, dcdxvec, takt, grampc->param->Nx, Nres, 1);

		(*pfct)(rhs, x, &takt, xakt, uakt, pvec, dcdxakt, grampc);
		for (i = 0; i < grampc->param->Nx; i++) {
			rhs[i] = -rhs[i];
		}
	}
	else {
		(*pfct)(rhs, x, &takt, xakt, uakt, pvec, dcdxakt, grampc);
	}
}

void dfdxRodas(int *N, typeRNum *t, typeRNum *x, typeRNum *dfdx_val, int *ldfdy, ctypeRNum *tvec, ctypeRNum *xvec, ctypeRNum *uvec, ctypeRNum *pvec, ctypeRNum *dcdxvec, const typeGRAMPC *grampc, const typeffctPtr pfct)
{
	typeRNum takt;
	typeRNum *uakt = grampc->rws->rwsGeneral + LWadjsys;  /* size NU */
	typeRNum *xakt = uakt + grampc->param->Nu;            /* size NX */

	ctypeInt pInt = grampc->rws->iparRodas[0];
	ctypeInt Nres = grampc->rws->iparRodas[1];

	if (pInt == BWINT) {
		takt = grampc->param->Thor - t[0];
	}
	else {
		takt = t[0];
	}

	interplin(uakt, tvec, uvec, takt, grampc->param->Nu, Nres, 1);

	if (pInt == FWINT) {
		dfdx(dfdx_val, takt, x, uakt, pvec, grampc->userparam);
	}
	else
	{
		interplin(xakt, tvec, xvec, takt, grampc->param->Nx, Nres, 1);
		dfdxtrans(dfdx_val, takt, xakt, uakt, pvec, grampc->userparam);
	}
}

void dfdtRodas(int *N, typeRNum *t, typeRNum *x, typeRNum *dfdt_val, int *ldfdy, ctypeRNum *tvec, ctypeRNum *xvec, ctypeRNum *uvec, ctypeRNum *pvec, ctypeRNum *dcdxvec, const typeGRAMPC *grampc, const typeffctPtr pfct)
{
	typeRNum takt;
	typeRNum *uakt = grampc->rws->rwsGeneral + LWadjsys;  /* size NU */
	typeRNum *xakt = uakt + grampc->param->Nu;            /* size NX */

	ctypeInt pInt = grampc->rws->iparRodas[0];
	ctypeInt Nres = grampc->rws->iparRodas[1];

	if (pInt == BWINT) {
		takt = grampc->param->Thor - t[0];
	}
	else {
		takt = t[0];
	}

	interplin(uakt, tvec, uvec, takt, grampc->param->Nu, Nres, 1);

	if (pInt == FWINT) {
		dfdt(dfdt_val, takt, x, uakt, pvec, grampc->userparam);
	}
	else
	{
		interplin(xakt, tvec, xvec, takt, grampc->param->Nx, Nres, 1);
		dHdxdt(dfdt_val, takt, xakt, uakt, x, pvec, grampc->userparam);
	}
}

void MfctRodas(int *N, typeRNum *out, int *LMAS, const typeGRAMPC *grampc, const typeffctPtr pfct)
{
	ctypeInt pInt = grampc->rws->iparRodas[0];

	if (pInt == FWINT) {
		Mfct(out, grampc->userparam);
	}
	else
	{
		Mtrans(out, grampc->userparam);
	}


}

void solout(int *nr, typeRNum *xold, typeRNum *x, typeRNum *h, typeRNum *y, typeRNum *cont, int *lrc, int *n, const typeGRAMPC *grampc, const typeffctPtr pfct, int *irtrn)
{
	/* Local variables */
	typeInt i;
	typeRNum dt = grampc->param->Thor / (typeRNum)(grampc->opt->Nhor - 1);

	if (*nr == 1)
	{
		for (i = 0; i < n[0]; i++) { grampc->rws->rparRodas[grampc->rws->iparRodas[2] + i] = y[i]; }
		grampc->rws->iparRodas[2] += 1;
	}
	else
	{
		while (*x >= grampc->rws->iparRodas[2] * dt)
		{
			for (i = 0; i < n[0]; i++) { grampc->rws->rparRodas[grampc->rws->iparRodas[2] * n[0] + i] = contro(i, n, grampc->rws->iparRodas[2] * dt, xold, h, cont, lrc); }
			grampc->rws->iparRodas[2] += 1;
		}
	}
} /* solout_ */