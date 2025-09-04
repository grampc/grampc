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


#include "timing.h"

#ifdef _WIN32

    void timer_now(typeTime* time)
    {
        QueryPerformanceCounter(time);
    }

    typeRNum timer_diff_ms(const typeTime* start, const typeTime* end)
    {
        LARGE_INTEGER frequency;
        QueryPerformanceFrequency(&frequency);
        return (end->QuadPart - start->QuadPart) * ((typeRNum)1000) / frequency.QuadPart;
    }

#elif defined(__linux__) || defined(__APPLE__)

    void timer_now(typeTime* time)
    {
        gettimeofday(time, 0);
    }

    typeRNum timer_diff_ms(const typeTime* start, const typeTime* end)
    {
        return (end->tv_sec - start->tv_sec) * ((typeRNum)1000)
            + (end->tv_usec - start->tv_usec) / ((typeRNum)1000);
    }

#else

    void timer_now(typeTime* time)
    {
        *time = clock();
    }

    typeRNum timer_diff_ms(const typeTime* start, const typeTime* end)
    {
        return (*end - *start) * ((typeRNum)1000) / CLOCKS_PER_SEC;
    }

#endif