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
from pygrampc import ProblemDescription

class Crane2D(ProblemDescription):
    def __init__(self, Q: np.array, R: np.array, ScaleConstraint, MaxConstraintHeight, MaxAngularDeflection):
        ProblemDescription.__init__(self, Nx=6, Nu=2, Np=0, Ng=0, Nh=3, NgT=0, NhT=0)

        self.Q = Q
        self.R = R
        self.ScaleConstraint = ScaleConstraint
        self.MaxConstraintHeight = MaxConstraintHeight
        self.MaxAngularDeflection = MaxAngularDeflection

    def ffct(self, out, t, x, u, p, param):
        out[0] = x[1]
        out[1] = u[0]
        out[2] = x[3]
        out[3] = u[1]
        out[4] = x[5]
        out[5] = -((9.81 * np.sin(x[4]) + np.cos(x[4]) * u[0] + 2 * x[3] * x[5]) / x[2])

    def dfdx_vec(self, out, t, x, u, p, vec, param):
        sinX = np.sin(x[4])
        cosX = np.cos(x[4])
        g = 9.81

        out[0] = 0
        out[1] = vec[0]
        out[2] = (g * sinX + cosX * u[0] + 2 * x[3] * x[5]) * vec[5] / x[2]**2
        out[3] = vec[2] - (2 * x[5] * vec[5]) / x[2]
        out[4] = -((g * cosX - sinX * u[0]) * vec[5] / x[2])
        out[5] = vec[4] - (2 * x[3] * vec[5]) / x[2]

    def dfdu_vec(self, out, t, x, u, p, vec, param):
        out[0] = vec[1] - (np.cos(x[4]) * vec[5]) / x[2]
        out[1] = vec[3]

    def lfct(self, out, t, x, u, p, param):
        out[0] = self.Q @ (x - param.xdes)**2 + self.R @ (u - param.udes)**2

    def dldx(self, out, t, x, u, p, param):
        out[:] = 2 * self.Q * (x - param.xdes)

    def dldu(self, out, t, x, u, p, param):
        out[:] = 2 * self.R * (u - param.udes)

    def hfct(self, out, t, x, u, p, param):
        Position = x[0] + np.sin(x[4]) * x[2]

        out[0] = np.cos(x[4]) * x[2] - self.ScaleConstraint * Position**2 - self.MaxConstraintHeight
        out[1] = x[5] - self.MaxAngularDeflection
        out[2] = -x[5] - self.MaxAngularDeflection

    def dhdx_vec(self, out, t, x, u, p, vec, param):
        tmp = self.ScaleConstraint * (x[0] + np.sin(x[4]) * x[2])

        out[0] = -2 * tmp * vec[0]
        out[1] = 0
        out[2] = (np.sin(x[4]) * tmp + np.cos(x[4])) * vec[0]
        out[3] = 0
        out[4] = (-2 * np.cos(x[4]) * x[2] * tmp - np.sin(x[4]) * x[2]) * vec[0]
        out[5] = 0 + vec[1] - vec[2]

    def dhdu_vec(self, out, t, x, u, p, vec, param):
        out[0] = 0