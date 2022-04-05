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



#include "ruku45.h"


void intsysRuKu45(typeRNum *y, ctypeInt pInt, ctypeInt Nint, ctypeRNum *t, ctypeRNum *x,
	ctypeRNum *u, ctypeRNum *p_, const typeGRAMPC *grampc, const typeffctPtr pfct)
{
	typeInt i;
	typeInt Nres = Nint;
	typeInt nstep = 0;
	typeRNum h, hnext;
	typeRNum *s1 = grampc->rws->rwsGeneral + LWadjsys; /* size:  18*NX + Nu */
	typeRNum *ys1 = s1 + grampc->param->Nx;
	typeRNum *s2 = ys1 + grampc->param->Nx;
	typeRNum *ys2 = s2 + grampc->param->Nx;
	typeRNum *s3 = ys2 + grampc->param->Nx;
	typeRNum *ys3 = s3 + grampc->param->Nx;
	typeRNum *s4 = ys3 + grampc->param->Nx;
	typeRNum *ys4 = s4 + grampc->param->Nx;
	typeRNum *s5 = ys4 + grampc->param->Nx;
	typeRNum *ys5 = s5 + grampc->param->Nx;
	typeRNum *s6 = ys5 + grampc->param->Nx;
	typeRNum *ys6 = s6 + grampc->param->Nx;
	typeRNum *s7 = ys6 + grampc->param->Nx;
	typeRNum *ys7 = s7 + grampc->param->Nx;

	ctypeRNum *tnow = t;
	ctypeRNum *xnow = x;
	ctypeRNum *unow = u;
	typeRNum *ynow = y;
	ctypeRNum tt = t[pInt*(Nint - 1)];
	ctypeRNum *tvec;
	ctypeRNum *uvec;
	ctypeRNum *zvec;

	/* only used for integrating the adjoint system */
	typeRNum *dcdxvec = grampc->rws->dcdx;
	typeRNum *dcdxnow = grampc->rws->dcdx + grampc->param->Nx * (grampc->opt->Nhor - 1);

	/* intermediate steps */
	typeRNum t1;
	typeRNum t2;
	typeRNum *y1 = ys7 + grampc->param->Nx;
	typeRNum *y2 = y1 + grampc->param->Nx;
	typeRNum *uakt = y2 + grampc->param->Nx;
	typeRNum *xakt = uakt + grampc->param->Nu;
	typeRNum *dcdxakt = xakt + grampc->param->Nx;
	typeRNum tdummy;

	/* Dormand-Prince-formula parameters */
	typeRNum a21 = (typeRNum)1.0 / (typeRNum)5.0,
		a31 = (typeRNum)3.0 / (typeRNum)40.0,
		a32 = (typeRNum)9.0 / (typeRNum)40.0,
		a41 = (typeRNum)44.0 / (typeRNum)45.0,
		a42 = (typeRNum)-56.0 / (typeRNum)15.0,
		a43 = (typeRNum)32.0 / (typeRNum)9.0,
		a51 = (typeRNum)19372.0 / (typeRNum)6561.0,
		a52 = (typeRNum)-25360.0 / (typeRNum)2187.0,
		a53 = (typeRNum)64448.0 / (typeRNum)6561.0,
		a54 = (typeRNum)-212.0 / (typeRNum)729.0,
		a61 = (typeRNum)9017.0 / (typeRNum)3168.0,
		a62 = (typeRNum)-355.0 / (typeRNum)33.0,
		a63 = (typeRNum)46732.0 / (typeRNum)5247.0,
		a64 = (typeRNum)49.0 / (typeRNum)176.0,
		a65 = (typeRNum)-5103.0 / (typeRNum)18656.0;
	typeRNum b1 = (typeRNum)35.0 / (typeRNum)384.0,
		b2 = (typeRNum)0.0,
		b3 = (typeRNum)500.0 / (typeRNum)1113.0,
		b4 = (typeRNum)125.0 / (typeRNum)192.0,
		b5 = (typeRNum)-2187.0 / (typeRNum)6784.0,
		b6 = (typeRNum)11.0 / (typeRNum)84.0,
		b7 = (typeRNum)0.0;
	typeRNum a71 = (typeRNum)35.0 / (typeRNum)384.0,
		a72 = (typeRNum)0.0,
		a73 = (typeRNum)500.0 / (typeRNum)1113.0,
		a74 = (typeRNum)125.0 / (typeRNum)192.0,
		a75 = (typeRNum)-2187.0 / (typeRNum)6784.0,
		a76 = (typeRNum)11.0 / (typeRNum)84.0;
	typeRNum b1s = (typeRNum)5179.0 / (typeRNum)57600.0,
		b2s = (typeRNum)0.0,
		b3s = (typeRNum)7571.0 / (typeRNum)16695.0,
		b4s = (typeRNum)393.0 / (typeRNum)640.0,
		b5s = (typeRNum)-92097.0 / (typeRNum)339200.0,
		b6s = (typeRNum)187.0 / (typeRNum)2100.0,
		b7s = (typeRNum)1.0 / (typeRNum)40.0;
	typeRNum db1 = b1 - b1s,
		db2 = b2 - b2s,
		db3 = b3 - b3s,
		db4 = b4 - b4s,
		db5 = b5 - b5s,
		db6 = b6 - b6s,
		db7 = b7 - b7s;
	typeRNum c2 = (typeRNum)1.0 / (typeRNum)5.0,
		c3 = (typeRNum)3.0 / (typeRNum)10.0,
		c4 = (typeRNum)4.0 / (typeRNum)5.0,
		c5 = (typeRNum)8.0 / (typeRNum)9.0,
		c6 = (typeRNum)1.0,
		c7 = (typeRNum)1.0;
	typeRNum d1 = (typeRNum)-12715105075.0 / (typeRNum)11282082432.0,
		d3 = (typeRNum)87487479700.0 / (typeRNum)32700410799.0,
		d4 = (typeRNum)-10690763975.0 / (typeRNum)1880347072.0,
		d5 = (typeRNum)701980252875.0 / (typeRNum)199316789632.0,
		d6 = (typeRNum)-1453857185.0 / (typeRNum)822651844.0,
		d7 = (typeRNum)69997945.0 / (typeRNum)29380423.0;

	typeInt reject = 0;

	typeRNum hmin = grampc->opt->IntegratorMinStepSize;						/*step size control accounts for minimal step size hmin*/

	typeRNum minscale = (typeRNum)1.0 / (typeRNum)5.0;
	typeRNum maxscale = (typeRNum)10.0;
	typeRNum errold = (typeRNum)0.0001;
	typeRNum bb = (typeRNum)0.08;
	typeRNum k = (typeRNum)5.0;
	typeRNum aa = (typeRNum)1 / k - bb * (typeRNum)0.75;
	typeRNum safe = (typeRNum)0.9;
	typeInt ende = 0;
	typeRNum scale;
	typeRNum sk;
	typeRNum err;
	typeRNum errPow;

	typeRNum theta;
	typeRNum interp1, interp2, interp3, interp4, interp5;

	if (pInt == FWINT) {
		tvec = t;
		uvec = u;
		zvec = x;
	}
	else {
		tvec = t - (Nint - 1);
		uvec = u - grampc->param->Nu*(Nint - 1);
		zvec = x - grampc->param->Nx*(Nint - 1);
	}

	/* initialization */
	t1 = tnow[0];
	for (i = 0; i < grampc->param->Nx; i++) {
		y1[i] = ynow[i];
	}
	for (i = 0; i < grampc->param->Nu; i++) {
		uakt[i] = unow[i];
	}
	if (y != x) {
		for (i = 0; i < grampc->param->Nx; i++) {
			xakt[i] = xnow[i];
		}
	}
	h = tnow[pInt] - tnow[0];

	/* Stage 1 */
	(*pfct)(s1, ynow, tnow, xnow, unow, p_, dcdxnow, grampc);

	while (pInt*(tt - tnow[pInt]) >= 0 && ende == 0) {
		/* START WHILE LOOP *******************************************************/

		if (nstep <= grampc->opt->IntegratorMaxSteps) {
			/* Stage 2 */
			for (i = 0; i < grampc->param->Nx; i++) {
				ys1[i] = y1[i] + h * (a21*s1[i]);
			}
			tdummy = t1 + c2 * h;
			interplin(uakt, tvec, uvec, tdummy, grampc->param->Nu, Nres, pInt);
			if (y != x) {
				interplin(xakt, tvec, zvec, tdummy, grampc->param->Nx, Nres, pInt);
			}
			interplin(dcdxakt, tvec, dcdxvec, tdummy, grampc->param->Nx, Nres, pInt);
			(*pfct)(s2, ys1, &tdummy, xakt, uakt, p_, dcdxakt, grampc);

			/* Stage 3 */
			for (i = 0; i < grampc->param->Nx; i++) {
				ys2[i] = y1[i] + h * (a31*s1[i] + a32 * s2[i]);
			}
			tdummy = t1 + c3 * h;
			interplin(uakt, tvec, uvec, tdummy, grampc->param->Nu, Nres, pInt);
			if (y != x) {
				interplin(xakt, tvec, zvec, tdummy, grampc->param->Nx, Nres, pInt);
			}
			interplin(dcdxakt, tvec, dcdxvec, tdummy, grampc->param->Nx, Nres, pInt);
			(*pfct)(s3, ys2, &tdummy, xakt, uakt, p_, dcdxakt, grampc);

			/* Stage 4 */
			for (i = 0; i < grampc->param->Nx; i++) {
				ys3[i] = y1[i] + h * (a41*s1[i] + a42 * s2[i] + a43 * s3[i]);
			}
			tdummy = t1 + c4 * h;
			interplin(uakt, tvec, uvec, tdummy, grampc->param->Nu, Nres, pInt);
			if (y != x) {
				interplin(xakt, tvec, zvec, tdummy, grampc->param->Nx, Nres, pInt);
			}
			interplin(dcdxakt, tvec, dcdxvec, tdummy, grampc->param->Nx, Nres, pInt);
			(*pfct)(s4, ys3, &tdummy, xakt, uakt, p_, dcdxakt, grampc);

			/* Stage 5 */
			for (i = 0; i < grampc->param->Nx; i++) {
				ys4[i] = y1[i] + h * (a51*s1[i] + a52 * s2[i] + a53 * s3[i] + a54 * s4[i]);
			}
			tdummy = t1 + c5 * h;
			interplin(uakt, tvec, uvec, tdummy, grampc->param->Nu, Nres, pInt);
			if (y != x) {
				interplin(xakt, tvec, zvec, tdummy, grampc->param->Nx, Nres, pInt);
			}
			interplin(dcdxakt, tvec, dcdxvec, tdummy, grampc->param->Nx, Nres, pInt);
			(*pfct)(s5, ys4, &tdummy, xakt, uakt, p_, dcdxakt, grampc);

			/* Stage 6 */
			for (i = 0; i < grampc->param->Nx; i++) {
				ys5[i] = y1[i] + h * (a61*s1[i] + a62 * s2[i] + a63 * s3[i] + a64 * s4[i] + a65 * s5[i]);
			}
			tdummy = t1 + c6 * h;
			interplin(uakt, tvec, uvec, tdummy, grampc->param->Nu, Nres, pInt);
			if (y != x) {
				interplin(xakt, tvec, zvec, tdummy, grampc->param->Nx, Nres, pInt);
			}
			interplin(dcdxakt, tvec, dcdxvec, tdummy, grampc->param->Nx, Nres, pInt);
			(*pfct)(s6, ys5, &tdummy, xakt, uakt, p_, dcdxakt, grampc);

			/* Stage 7 */
			for (i = 0; i < grampc->param->Nx; i++) {
				ys6[i] = y1[i] + h * (a71*s1[i] + a72 * s2[i] + a73 * s3[i] + a74 * s4[i] + a75 * s5[i] + a76 * s6[i]);
			}
			tdummy = t1 + c7 * h;
			interplin(uakt, tvec, uvec, tdummy, grampc->param->Nu, Nres, pInt);
			if (y != x) {
				interplin(xakt, tvec, zvec, tdummy, grampc->param->Nx, Nres, pInt);
			}
			interplin(dcdxakt, tvec, dcdxvec, tdummy, grampc->param->Nx, Nres, pInt);
			(*pfct)(s7, ys6, &tdummy, xakt, uakt, p_, dcdxakt, grampc);

			/* Error */
			t2 = t1 + h;
			err = 0.0;

			for (i = 0; i < grampc->param->Nx; i++) {
				y2[i] = y1[i] + h * (b1*s1[i] + b2 * s2[i] + b3 * s3[i] + b4 * s4[i] + b5 * s5[i] + b6 * s6[i]);
				sk = grampc->opt->IntegratorAbsTol + grampc->opt->IntegratorRelTol*MAX(ABS(y2[i]), ABS(y1[i]));
				errPow = h * (db1*s1[i] + db2 * s2[i] + db3 * s3[i] + db4 * s4[i] + db5 * s5[i] + db6 * s6[i] + db7 * s7[i]) / sk;
				err += errPow * errPow;
			}
			err = SQRT(err / (grampc->param->Nx));

			/* step size control */
			if (err <= 1.0) {
				if (err == 0.0) {
					scale = maxscale;
				}
				else {
					scale = safe * POW(err, -aa)*POW(errold, bb);
					if (scale < minscale) {
						scale = minscale;
					}
					if (scale > maxscale) {
						scale = maxscale;
					}
				}
				if (reject == 1) {
					hnext = h * MIN(scale, 1);
				}
				else {
					hnext = h * scale;
				}
				errold = MAX(err, (typeRNum)0.0001);
				reject = 0;
			} /*if (err>1.0)*/
			else {
				scale = MAX(safe*POW(err, -aa), minscale);
				h = h * scale;
				reject = 1;

				if (ABS(h) < hmin) {
					h = (typeRNum)pInt*hmin;
					hnext = h;
					reject = 0;
					grampc->sol->status |= STATUS_INTEGRATOR_H_MIN;
				}
			}
			if (reject == 0) {
				while (pInt*(t2 - tnow[pInt]) >= 0 && ende == 0) {
					theta = (tnow[pInt] - t1) / (t2 - t1);
					for (i = 0; i < grampc->param->Nx; i++) {
						interp1 = y1[i];
						interp2 = y2[i] - y1[i];
						interp3 = h * s1[i] - interp2;
						interp4 = interp2 - h * s7[i] - interp3;
						interp5 = h * (d1*s1[i] + d3 * s3[i] + d4 * s4[i] + d5 * s5[i] + d6 * s6[i] + d7 * s7[i]);
						ynow[i + pInt * grampc->param->Nx] = interp1 + theta * (interp2 + (1 - theta)*(interp3 + theta * (interp4 + (1 - theta)*interp5)));
					}
					if (tnow[pInt] == tt) {
						ende = 1;
					}
					else {
						Nres -= 1;
						tnow += pInt;
						xnow += pInt * grampc->param->Nx;
						unow += pInt * grampc->param->Nu;
						ynow += pInt * grampc->param->Nx;
						if (pInt == FWINT) {
							tvec += 1;
							zvec += grampc->param->Nx;
							uvec += grampc->param->Nu;
						}
					}
				}
				/* Rewriting of local time */
				t1 = t2;

				/* Ensuring that the time remains within the integration interval. */
				if (pInt* (t1+hnext) > tt) {
					h = tt - t1;
				}
				else {
					h = hnext;
				}

				/* Rewriting of local state variables */
				for (i = 0; i < grampc->param->Nx; i++) {
					y1[i] = y2[i];
				}
				for (i = 0; i < grampc->param->Nx; i++) {
					s1[i] = s7[i]; /* Using the previous last stage for the next step */
				}
			}
		}
		else
		{
			ende = 1;
			grampc->sol->status |= STATUS_INTEGRATOR_MAXSTEPS;
		}
		nstep++;
		/* END WHILE LOOP *********************************************************/
	}
}
