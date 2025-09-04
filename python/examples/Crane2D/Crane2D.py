"""
This file is part of GRAMPC - (https://github.com/grampc/grampc)

GRAMPC -- A software framework for embedded nonlinear model predictive
control using a gradient-based augmented Lagrangian approach

Copyright 2014-2025 by Knut Graichen, Andreas Voelz, Thore Wietzke,
Tobias Englert (<v2.3), Felix Mesmer (<v2.3), Soenke Rhein (<v2.3),
Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
All rights reserved.

GRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
"""
import numpy as np
import matplotlib.pyplot as plt
from pygrampc import Grampc, GrampcResults
from scipy.integrate import solve_ivp

# Compare the speed and implementation between Python and C++
# Before using the C++ problem description, one has to install it.
# The easiest way is via pip with the supplied pyproject.toml and CMakeLists.txt file

from crane_problem import Crane2D
# from crane_problem_py import Crane2D

if __name__ == "__main__":
    Tsim = 12.5
    plotSteps = 150
    plotPause = False
    options = "Crane2D.json"

    Q = np.array([1.0, 2.0, 2.0, 1.0, 1.0, 4.0])
    R = np.array([0.05, 0.05])

    # initialize problem and GRAMPC
    problem = Crane2D(Q, R, 0.2, 1.25, 0.3)
    grampc = Grampc(problem, options, plot_prediction=False)

    # estimate penaltyMin and set option
    grampc.estim_penmin(True)
    grampc.print_opts()
    grampc.print_params()

    # construct solution structure
    vec = GrampcResults(grampc, Tsim, plot_results=True, plot_statistics=True)

    dt = grampc.param.dt

    for i, t in enumerate(vec.t):
        vec.CPUtime[i] = grampc.run()
        vec.update(grampc, i)

        if i + 1 > len(vec.t):
            break

        # simulate system
        sol = solve_ivp(grampc.ffct, [t, t+dt], grampc.param.x0,
                        args=(grampc.sol.unext, grampc.sol.pnext, grampc.param))

        # set current time and state
        grampc.set_param({"x0": sol.y[:, -1],
                          "t0": t + dt})

        # plots of the grampc predictions
        if i % plotSteps == 0:
            grampc.plot()
            vec.plot()
            if plotPause:
                input("Press Enter to continue...")

    vec.plot()
    plt.show()
