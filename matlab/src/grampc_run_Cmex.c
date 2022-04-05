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

#include "mex.h"
#include "grampc.h"
#include "grampc_conversion_Cmex.h"

#ifndef N_TIMER
#define N_TIMER 1
#endif

/* Runtime measurement */
#ifdef N_TIMER
#ifdef _WIN32
#include <Windows.h>
#elif defined(__linux__) || defined(__APPLE__)
#include <sys/time.h> 
#endif
#endif


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	typeGRAMPC* grampc;
	typeRNum * comptime;

	mxArray *mxuserparam = NULL;

	/* variables for runtime measurement */
#ifdef N_TIMER
	typeInt i;
	typeGRAMPC* grampcArray[N_TIMER];
#ifdef _WIN32
	LARGE_INTEGER StartingTime, EndingTime, ElapsedMicroseconds;
	LARGE_INTEGER Frequency;
#elif defined(__linux__) || defined(__APPLE__)
	struct timeval StartingTime, EndingTime;
#endif
#endif

	/* check proper number of input arguments */
	if (nrhs < 1) {
		mexErrMsgTxt("Not enough input arguments.");
	}
	if (nrhs > 1) {
		mexErrMsgTxt("Too many input arguments.");
	}
	/* first argument: param structure */
	if (!mxIsStruct(prhs[0])) {
		mexErrMsgTxt("First argument (grampc structure) must be a struct.");
	}
	/* check proper number of output arguments */
	if (nlhs > 7) {
		mexErrMsgTxt("Too many output arguments.");
	}

	/* Allocate grampc struct and copy mxData */
	mxuserparam = mxDuplicateArray(mxGetField(prhs[0], 0, "userparam"));
	mx2typeGRAMPC(&grampc, prhs[0], mxuserparam);

	/* Allocate memory for the outputs */
	plhs[1] = mxCreateNumericMatrix(1, 1, mxtypeRNum_CLASS, mxREAL);
	comptime = (typeRNum*)mxGetPr(plhs[1]);
	*comptime = -1;

	/* runtime analysis */

#ifdef N_TIMER
	/* copy GRAMPC struct for runtime analysis */
	for (i = 0; i < N_TIMER; i++) {
		mx2typeGRAMPC((grampcArray + i), prhs[0], mxuserparam);
	}

	/* initialize variables for runtime check */
#ifdef _WIN32
	QueryPerformanceFrequency(&Frequency);
	QueryPerformanceCounter(&StartingTime);
#elif defined(__linux__) || defined(__APPLE__)
	gettimeofday(&StartingTime, NULL);
#endif

	/* run GRAMPC while measuring time */
	grampc_run(grampc);

	for (i = 0; i < N_TIMER-1; i++) {
		grampc_run(grampcArray[i]);
	}

	/* runtime evaluation */
#ifdef _WIN32
	QueryPerformanceCounter(&EndingTime);
	ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
	ElapsedMicroseconds.QuadPart *= (100000000 / N_TIMER); /* 100ns */
	ElapsedMicroseconds.QuadPart /= (Frequency.QuadPart);
	*comptime = ((typeRNum)(ElapsedMicroseconds.QuadPart) / 100000); /* time in ms */
#elif defined(__linux__) || defined(__APPLE__)
	gettimeofday(&EndingTime, NULL);
	*comptime = (EndingTime.tv_sec - StartingTime.tv_sec) * 1000000.0;      /* sec to µs */
	*comptime += (EndingTime.tv_usec - StartingTime.tv_usec) * 1.0;         /* µs to µs */
	*comptime /= 1000.0;																										/* µs to ms */
	*comptime /= N_TIMER; /* Time of one MPC step */
#endif

	/* free copies of GRAMPC struct */
	for (i = 0; i < N_TIMER; i++) {
		grampc_free((grampcArray + i));
	}
#endif

	/* create mx structure */
	typeGRAMPC2mx(plhs, grampc, mxuserparam);

	/* free allocated memory */
	grampc_free(&grampc);
}
