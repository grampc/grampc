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


#ifndef TIMING_H
#define TIMING_H

#include "grampc_macro.h"

#ifdef _WIN32

    #define WIN32_LEAN_AND_MEAN
    #include <Windows.h>

    typedef LARGE_INTEGER typeTime;

#elif defined(__linux__) || defined(__APPLE__)

    #include <sys/time.h> 

    typedef struct timeval typeTime;

#else

    #include <time.h>

    typedef clock_t typeTime;

#endif

/* Get current timestamp */
void timer_now(typeTime* time);

/* Get elapsed time between two timestamps in milliseconds  */
typeRNum timer_diff_ms(const typeTime* start, const typeTime* end);

#endif