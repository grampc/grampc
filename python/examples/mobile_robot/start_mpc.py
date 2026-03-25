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
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from pygrampc import Grampc, GrampcResults
from grampc_mobile_robot import MobileRobotProblem

Parameters = {
    "x0": [0, 0, 0],
    "umax": [2, 2],
    "umin": [-2,-2],
    "Thor": 1.0,
    "dt": 0.01
}

Options = {
    "Nhor": 30,
    "MaxGradIter": 3,
    "MaxMultIter": 1,
    "InequalityConstraints": "on"
}

if __name__ == "__main__":
    plotSteps = 10

    # problem-specific parameters
    pSys = []
    pCost = [1, 5, 0.1, 0.125, 0.0125]
    pCon = [-0.9, 0.9, -0.9, 0.9]

    # create problem description and solver
    problem = MobileRobotProblem(pSys, pCost, pCon)
    grampc = Grampc(problem)

    # set parameters and options for GRAMPC
    grampc.set_param(Parameters)
    grampc.set_opt(Options)
    grampc.estim_penmin(True)
    grampc.check_gradients(1e-6, 1e-8)

    # MPC main loop
    Tsim = 5.5
    vec = GrampcResults(grampc, Tsim)

    # prepare plots (interactive mode)
    plt.ion()
    fig, axs = plt.subplots(1, 2)
    axs[0].add_patch(plt.Circle((0, 0), 1, fill=False, linestyle='--'))
    axs[0].add_patch(plt.Rectangle((-0.9, -0.9), 1.8, 1.8, fill=False, linestyle=':'))
    ph_vec_pos = axs[0].plot(vec.x[:, 0], vec.x[:, 1])[0]
    ph_pred_pos = axs[0].plot(grampc.rws.x[0, :], grampc.rws.x[1, :])[0]
    axs[0].axis('equal')
    axs[0].set_xlabel('Position x')
    axs[0].set_ylabel('Position y')

    axs[1].plot([0, Tsim+grampc.param.Thor], [0.9, 0.9], 'k--')
    axs[1].plot([0, Tsim+grampc.param.Thor], [-0.9, -0.9], 'k--')
    ph_vec_x = axs[1].plot(vec.t, vec.x)
    axs[1].set_prop_cycle(None)
    ph_pred_x = axs[1].plot(grampc.rws.t, grampc.rws.x.T, '.')
    axs[1].set_prop_cycle(None)
    axs[1].set_xlim(0, Tsim+grampc.param.Thor)
    axs[1].set_ylim(-1.1, 1.1)
    axs[1].set_xlabel('Time t')
    axs[1].set_ylabel('States x')

    for i, t in enumerate(vec.t):
        # solve problem
        vec.CPUtime[i] = grampc.run()
        vec.update(grampc, i)
        vec.x[i, :] = grampc.param.x0 # store current state instead of predicted state

        # simulate system
        sol = solve_ivp(grampc.ffct, [t, t + grampc.param.dt], grampc.param.x0,
                        args=(grampc.sol.unext, grampc.param.p0, grampc.param))
        next_state = sol.y[:grampc.param.Nx, -1]

	    # update current time and state
        grampc.set_param({"t0": t + grampc.param.dt, "x0": next_state})
        
        # update plots
        if i % plotSteps == 0:
            ph_vec_pos.set_data(vec.x[:, 0], vec.x[:, 1])
            ph_pred_pos.set_data(grampc.rws.x[0, :], grampc.rws.x[1, :])
            for j in range(grampc.param.Nx):
                ph_vec_x[j].set_data(vec.t, vec.x[:, j])
                ph_pred_x[j].set_data(t + grampc.rws.t, grampc.rws.x[j, :])
            plt.pause(grampc.param.dt)

print(f'Computation time: mean {np.mean(vec.CPUtime):.3f} ms and max {np.max(vec.CPUtime):.3f} ms')

plt.ioff()
plt.show()