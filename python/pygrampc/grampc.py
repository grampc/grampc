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
import json
from math import nan, inf
from typing import NamedTuple, Optional
import matplotlib.pyplot as plt
import numpy as np

from ._core import GrampcBinding, ProblemDescription


class ValidatorOptions(NamedTuple):
    allowedValues: tuple
    Class: str
    dimension: int = 0


class Grampc(GrampcBinding):
    """
    An interface for the GRAMPC solver. Handles the initialization and destruction of the GRAMPC solver.

    Attributes:
        problem: The problem description for this instance.
        fig: handle for the prediction plot figure.
        param: The GRAMPC parameters struct.
        opt: The GRAMPC options struct.
        rws: The GRAMPC real workspace struct.
        sol: The GRAMPC solution struct.
    """
    _parameters_real = {
        "dt": ValidatorOptions((0.0, inf), "OpenInterval"),
        "Thor": ValidatorOptions((0.0, inf), "OpenInterval"),
        "Tmax": ValidatorOptions((0.0, inf), "OpenInterval"),
        "Tmin": ValidatorOptions((0.0, inf), "OpenInterval"),
        "t0": ValidatorOptions((-inf, inf), "OpenInterval")}

    _options_real = {
        "IntegratorRelTol": ValidatorOptions((0.0, inf), "OpenInterval"),
        "IntegratorAbsTol": ValidatorOptions((0.0, inf), "OpenInterval"),
        "IntegratorMinStepSize": ValidatorOptions((0.0, inf), "OpenInterval"),
        "LineSearchMax": ValidatorOptions((0.0, inf), "OpenInterval"),
        "LineSearchMin": ValidatorOptions((0.0, inf), "OpenInterval"),
        "LineSearchInit": ValidatorOptions((0.0, inf), "OpenInterval"),
        "LineSearchAdaptAbsTol": ValidatorOptions((0.0, inf), "LeftClosedInterval"),
        "LineSearchAdaptFactor": ValidatorOptions((1.0, inf), "OpenInterval"),
        "LineSearchIntervalTol": ValidatorOptions((0.0, 0.5), "OpenInterval"),
        "LineSearchIntervalFactor": ValidatorOptions((0.0, 1.0), "OpenInterval"),
        "OptimParamLineSearchFactor": ValidatorOptions((0.0, inf), "OpenInterval"),
        "OptimTimeLineSearchFactor": ValidatorOptions((0.0, inf), "OpenInterval"),
        "TScale": ValidatorOptions((0.0, inf), "OpenInterval"),
        "TOffset": ValidatorOptions((-inf, inf), "OpenInterval"),
        "JScale": ValidatorOptions((0.0, inf), "OpenInterval"),
        "MultiplierMax": ValidatorOptions((0.0, inf), "OpenInterval"),
        "MultiplierDampingFactor": ValidatorOptions((0.0, 1.0), "ClosedInterval"),
        "PenaltyMax": ValidatorOptions((0.0, inf), "OpenInterval"),
        "PenaltyMin": ValidatorOptions((0.0, inf), "OpenInterval"),
        "PenaltyIncreaseFactor": ValidatorOptions((1.0, inf), "LeftClosedInterval"),
        "PenaltyDecreaseFactor": ValidatorOptions((0.0, 1.0), "ClosedInterval"),
        "PenaltyIncreaseThreshold": ValidatorOptions((0.0, inf), "LeftClosedInterval"),
        "AugLagUpdateGradientRelTol": ValidatorOptions((0.0, 1.0), "ClosedInterval"),
        "ConvergenceGradientRelTol": ValidatorOptions((0.0, 1.0), "ClosedInterval")}

    _options_int = {
        "Nhor": ValidatorOptions((2, inf), "LeftClosedInterval"),
        "MaxGradIter": ValidatorOptions((1, inf), "LeftClosedInterval"),
        "MaxMultIter": ValidatorOptions((1, inf), "LeftClosedInterval"),
        "IntegratorMaxSteps": ValidatorOptions((1, inf), "LeftClosedInterval")}

    _options_str = {
        "ShiftControl": ValidatorOptions(("on", "off"), "Enum"),
        "IntegralCost": ValidatorOptions(("on", "off"), "Enum"),
        "TerminalCost": ValidatorOptions(("on", "off"), "Enum"),
        "IntegratorCost": ValidatorOptions(("trapezoidal", "simpson"), "Enum"),
        "Integrator": ValidatorOptions(("erk1", "erk2", "erk3", "erk4", "ruku45", "rodas", "discrete"), "Enum"),
        "LineSearchType": ValidatorOptions(("adaptive", "explicit1", "explicit2"), "Enum"),
        "LineSearchExpAutoFallback": ValidatorOptions(("on", "off"), "Enum"),
        "OptimControl": ValidatorOptions(("on", "off"), "Enum"),
        "OptimParam": ValidatorOptions(("on", "off"), "Enum"),
        "OptimTime": ValidatorOptions(("on", "off"), "Enum"),
        "ScaleProblem": ValidatorOptions(("on", "off"), "Enum"),
        "EqualityConstraints": ValidatorOptions(("on", "off"), "Enum"),
        "InequalityConstraints": ValidatorOptions(("on", "off"), "Enum"),
        "TerminalEqualityConstraints": ValidatorOptions(("on", "off"), "Enum"),
        "TerminalInequalityConstraints": ValidatorOptions(("on", "off"), "Enum"),
        "ConstraintsHandling": ValidatorOptions(("auglag", "extpen"), "Enum"),
        "ConvergenceCheck": ValidatorOptions(("on", "off"), "Enum")}

    def __init__(self, problem: ProblemDescription, 
                       options_path: Optional[str] = None, 
                       plot_prediction: Optional[bool] = False):
        """
        Initilizes the interface and the underlying GRAMPC object.

        Args:
            problem: The concrete problem description object.
            options_path: The path to the options.json file.
            plot_prediction: If True enables the prediction plot utility.
        """
        super().__init__(problem)

        self._parameters_real_vec = {
            "x0": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nx),
            "xdes": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nx),
            "u0": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nu),
            "udes": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nu),
            "umax": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nu),
            "umin": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nu),
            "p0": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Np),
            "pmax": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Np),
            "pmin": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Np)
        }

        self._options_real_vec = {
            "xScale": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nx),
            "xOffset": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nx),
            "uScale": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nu),
            "uOffset": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nu),
            "pScale": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Np),
            "pOffset": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Np),
            "cScale": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nc),
            "ConstraintsAbsTol": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nc)
        }

        self._options_int_vec = {"FlagsRodas": ValidatorOptions((0, self.param.Nx), "Rodas", 8)}

        if options_path is not None:
            with open(options_path) as f:
                grampc_options = json.load(f)
            self.set_param(grampc_options["Parameters"])
            self.set_opt(grampc_options["Options"])

        if plot_prediction:
            self.fig, _ = plt.subplots(2, 3, sharex=True, figsize=(10, 5))

        else:
            self.fig = None

    @staticmethod
    def _len(value):
        "Returns the length of a list or numpy array. For scalar inputs the function returns 0"
        try:
            return len(value)
        except TypeError:
            return 0

    def _validator(self, key: str, value, options: ValidatorOptions) -> bool:
        """
        Helper function to validate the inputs for the GRAMPC parameters and options.

        Args:
            key: The option/parameter key to check.
            value: The value of the option/parameter.
            options: NamedTuple with the option class and allowed values.

        Returns:
            bool: True if value is valid.
        
        Raises:
            ValueError: Raised if value is not valid.
            ValueError: Raised if dimension does not match.
        """
        valid_input = False
        if options.Class == "Enum":
            valid_input = value in options.allowedValues
        elif options.Class == "OpenInterval":
            if self._len(value) == options.dimension:
                if isinstance(value, (list, np.ndarray)):
                    valid_input = all(
                        [options.allowedValues[0] < item < options.allowedValues[1] for item in value])
                else:
                    valid_input = options.allowedValues[0] < value < options.allowedValues[1]
            else:
                raise ValueError(
                    f"{key}: dimension does not match problem description. Expected {options.dimension}, "
                    f"got {self._len(value)}"
                )

        elif options.Class == "LeftClosedInterval":
            if self._len(value) == options.dimension:
                if isinstance(value, (list, np.ndarray)):
                    valid_input = all(
                        [options.allowedValues[0] <= item < options.allowedValues[1] for item in value])
                else:
                    valid_input = options.allowedValues[0] <= value < options.allowedValues[1]
            else:
                raise ValueError(
                    f"{key}: dimension does not match problem description. Expected {options.dimension}, "
                    f"got {self._len(value)}"
                )

        elif options.Class == "ClosedInterval":
            if self._len(value) == options.dimension:
                if isinstance(value, (list, np.ndarray)):
                    valid_input = all(
                        [options.allowedValues[0] <= item <= options.allowedValues[1] for item in value])
                else:
                    valid_input = options.allowedValues[0] <= value <= options.allowedValues[1]
            else:
                raise ValueError(
                    f"{key}: dimension does not match problem description. Expected {options.dimension}, "
                    f"got {self._len(value)}"
                )

        elif options.Class == "Rodas":
            if self._len(value) == 8:
                # Check if IFCN, IDFX, IJAC and IMAS are either 1 or 0
                valid_input = all([value[i] == 1 or value[i] == 0 for i in range(4)])
                # check if MLJAC, MUJAX, MLMAS and MUMAS are between 0 and Nx
                valid_input = valid_input and all(
                    [options.allowedValues[0] <= value[i + 4] <= options.allowedValues[1] for i in range(4)])
            else:
                raise ValueError(
                    f"{key}: dimension does not match problem description. Expected 8, got {self._len(value)}"
                )

        if not valid_input:
            raise ValueError(f"{key}: unallowed value detected. Expected {options.allowedValues}, got {value}")
        return valid_input

    def set_param(self, parameters: dict):
        """
        Sets the GRAMPC parameters with value validation.

        Args:
            parameters: Mapping with parameter name and values.
        
        Raises: 
            KeyError: if a parameter key is not a valid GRAMPC parameter.
            ValueError: Raised if value is not valid.
            ValueError: Raised if dimension does not match.
        """
        for key, value in parameters.items():
            if key in self._parameters_real:
                if self._validator(key, value, self._parameters_real[key]):
                    self._set_param_real(key, value)
            elif key in self._parameters_real_vec:
                if self._validator(key, value, self._parameters_real_vec[key]):
                    self._set_param_real_vec(key, value)
            else:
                raise KeyError(f"{key} not found as GRAMPC parameter")

    def set_opt(self, options: dict):
        """
        Sets the GRAMPC options with value validation.

        Args:
            options: Mapping with parameter name and values.
        
        Raises: 
            KeyError: if a parameter key is not a valid GRAMPC parameter.
            ValueError: Raised if value is not valid.
            ValueError: Raised if dimension does not match.
        """
        for key, value in options.items():
            if key in self._options_str:
                if self._validator(key, value, self._options_str[key]):
                    self._set_opt_str(key, value)
            elif key in self._options_int:
                if self._validator(key, value, self._options_int[key]):
                    self._set_opt_int(key, value)
            elif key in self._options_int_vec:
                self._set_opt_int_vec(key, value)
            elif key in self._options_real:
                if self._validator(key, value, self._options_real[key]):
                    self._set_opt_real(key, value)
            elif key in self._options_real_vec:
                if self._validator(key, value, self._options_real_vec[key]):
                    self._set_opt_real_vec(key, value)
            else:
                raise KeyError(f"{key} not found as GRAMPC option")

    def plot(self):
        """
        Draws the current prediction plot with the values of the rws struct.
        """
        if self.fig is not None:

            self.fig.suptitle("Prediction")
            axs = self.fig.axes
            for ax in axs:
                ax.clear()
                ax.grid()

            axs[0].plot(self.rws.t[:, None], self.rws.x.T)
            axs[0].set_title("Predicted states")

            axs[1].plot(self.rws.t[:, None], self.rws.adj.T)
            axs[1].set_title("Predicted adjoint states")

            axs[2].plot(self.rws.t[:, None], self.rws.u.T)
            axs[2].set_title("Predicted controls")

            if self.param.Nc == self.param.NgT + self.param.NhT:
                axs[3].scatter(np.full((self.param.Nc, 1), self.rws.t[-1]), self.rws.cfct[:, -1][:, None])
                axs[4].scatter(np.full((self.param.Nc, 1), self.rws.t[-1]), self.rws.mult[:, -1][:, None])
                axs[5].scatter(np.full((self.param.Nc, 1), self.rws.t[-1]), self.rws.pen[:, -1][:, None])
            else:
                axs[3].plot(self.rws.t[:, None], self.rws.cfct.T)
                axs[4].plot(self.rws.t[:, None], self.rws.mult.T)
                axs[5].plot(self.rws.t[:, None], self.rws.pen.T)
            axs[3].set_title("Predicted constraints")
            axs[4].set_title("Predicted Lagrange multipliers")
            axs[5].set_title("Predicted penalty parameters")

            axs[0].set_xlim(left=self.rws.t[0])

            plt.draw()
            plt.pause(0.0001)


class GrampcResults:
    """
    Dataclass for storing the results of the GRAMPC solver.

    Attributes:
        t: Time array.
        x: State array.
        u: Control array.
        p: Parameter array.
        adj: Adjoint states array.
        J: Costs and augmented costs.
        T: Prediction Horizon.
        iter: Number of gradient iterations per multiplier iteration.
        Nconstr: Norm of the constraint violation of the prediction horizon.
        Npen: Norm of penalty parameters over horizon.
        CPUtime: CPU wall clock time for every time step.
        linesearch_step_size: Step size of the line search.
        constr: Evaluated equality and inequality constraints.
        mult: Lagrange Multipliers for the constraints.
        pen: Penalty parameters for the constraints.
        constrT: Evaluated Terminal equality and inequality constraints.
        multT: Lagrange Multipliers for the terminal constraints.
        penT: Penalty parameters for the terminal constraints.
    """

    class index:
        def __init__(self, grampc):
            self.g = grampc.param.Ng
            self.h = grampc.param.Nh + self.g
            self.gT = grampc.param.NgT + self.h
            self.hT = grampc.param.NhT + self.gT

    def __init__(self, grampc: Grampc, 
                       Tsim: float, 
                       plot_statistics: Optional[bool] = False, 
                       plot_results: Optional[bool] = False):
        """
        Initilizes the result data container for GRAMPC. 
        Creates numpy arrays with length (Tsim/dt + 1) and fills them with NaN

        Args:
            grampc: Reference to the GRAMPC interface.
            Tsim: End time of your simulation.
            plot_statistics: If True enables the statistics plot utility.
            plot_results: If True enables the results plot utility.
        """

        num_data_points = int(Tsim / grampc.param.dt) + 1
        self._id = self.index(grampc)
        self.t = np.linspace(grampc.param.t0, Tsim, num_data_points)
        self.x = np.full((num_data_points, grampc.param.Nx), nan)
        self.u = np.full((num_data_points, grampc.param.Nu), nan)
        self.p = np.full((num_data_points, grampc.param.Np), nan)
        self.adj = np.full((num_data_points, grampc.param.Nx), nan)
        self.J = np.full((num_data_points, 2), nan)
        self.T = np.full((num_data_points, 1), nan)
        self.iter = np.full((num_data_points, grampc.opt.MaxMultIter), nan)

        self.CPUtime = np.full((num_data_points, 1), nan)
        self.linesearch_step_size = np.full((num_data_points, 1), nan)
        self._constraints = (grampc.param.Ng + grampc.param.Nh) > 0
        self._terminal_constraints = (grampc.param.NgT + grampc.param.NhT) > 0
        self._OptimTime = grampc.opt.OptimTime == 1
        self._OptimParam = grampc.opt.OptimParam == 1

        if self._constraints or self._terminal_constraints:
            self.Nconstr = np.full((num_data_points, 1), nan)
            self.Npen = np.full((num_data_points, 1), nan)

        if self._constraints:
            self.constr = np.full((num_data_points, grampc.param.Ng + grampc.param.Nh), nan)
            self.mult = np.full((num_data_points, grampc.param.Ng + grampc.param.Nh), nan)
            self.pen = np.full((num_data_points, grampc.param.Ng + grampc.param.Nh), nan)

        if self._terminal_constraints:
            self.constrT = np.full((num_data_points, grampc.param.NgT + grampc.param.NhT), nan)
            self.multT = np.full((num_data_points, grampc.param.NgT + grampc.param.NhT), nan)
            self.penT = np.full((num_data_points, grampc.param.NgT + grampc.param.NhT), nan)

        if plot_statistics:
            Nrows = 1
            Ncols = 0
            if grampc.opt.OptimTime or grampc.opt.OptimParam:
                Nrows = 2
                Ncols = grampc.opt.OptimTime + grampc.opt.OptimParam

            if self._constraints or self._terminal_constraints:
                if self._terminal_constraints:
                    Nrows = 3
                else:
                    Nrows = 2
                Ncols += 2

            self._stat_Nrows = Nrows
            self._stat_Ncols = Ncols
            self._fig_statistics = plt.figure()
        else:
            self._fig_statistics = None

        if plot_results:
            self._fig_results, _ = plt.subplots(1 + self._constraints, 3, sharex=True)
        else:
            self._fig_results = None

    def update(self, grampc: Grampc, index: int):
        """
        Inserts the current datapoints into the numpy arrays.
        Evaluates the constraint functions of GRAMPC.
        
        Args:
            grampc: Reference to the GRAMPC interface.
            index: Current time index. 
        """
        self.x[index, :] = grampc.sol.xnext
        self.u[index, :] = grampc.sol.unext
        self.adj[index, :] = grampc.rws.adj[:, 0] / grampc.opt.xScale
        self.J[index, :] = grampc.sol.J
        self.T[index] = grampc.sol.Tnext
        if self._constraints or self._terminal_constraints:
            self.Nconstr[index] = grampc.sol.cfct
            self.Npen[index] = grampc.sol.pen
        self.iter[index, :] = grampc.sol.iter

        if grampc.opt.LineSearchType == 0:
            self.linesearch_step_size[index] = grampc.rws.lsAdapt[3 + grampc.opt.MaxGradIter * 8]
        elif grampc.opt.LineSearchType >= 1:
            self.linesearch_step_size[index] = grampc.rws.lsExplicit[2]

        # retrieve the equality and inequality constraints
        if self._constraints:
            if grampc.param.Ng:
                self.constr[index, 0:self._id.g] = grampc.gfct(0.0, self.x[index, :], self.u[index, :], self.p[index, :], grampc.param)

            if grampc.param.Nh:
                self.constr[index, self._id.g:self._id.h] = grampc.hfct(self.t[index], self.x[index, :], self.u[index, :], self.p[index, :], grampc.param)

            self.mult[index, :] = grampc.rws.mult[0:self._id.h, 0]
            self.pen[index, :] = grampc.rws.pen[0:self._id.h, 0]

        if self._terminal_constraints:
            if grampc.param.NgT:
                self.constrT[index, 0:grampc.param.NgT] = grampc.gTfct(self.t[index], self.x[index, :], self.p[index, :], grampc.param)

            if grampc.param.NhT:
                self.constrT[index, grampc.param.NgT:grampc.param.Ng+grampc.param.NhT] = grampc.hTfct(self.t[index], self.x[index, :], self.p[index, :], grampc.param)

            self.multT[index, :] = grampc.rws.mult[self._id.h:self._id.hT, -1]
            self.penT[index, :] = grampc.rws.pen[self._id.h:self._id.hT, -1]

    def plot(self):
        """
        Utility function to plot the statistics and results stored in this class.
        """
        if self._fig_statistics is not None:
            self._fig_statistics.clear()
            self._fig_statistics.suptitle("Statistics")

            ax = self._fig_statistics.add_subplot(self._stat_Nrows, 3, 1)
            ax.plot(self.t, self.J)
            ax.set_title("Costs")
            ax.legend(["Original costs", "Augmented costs"])

            ax = self._fig_statistics.add_subplot(self._stat_Nrows, 3, 2, sharex=ax)
            ax.plot(self.t, self.CPUtime)
            ax.set_title("Computation time in ms")

            ax = self._fig_statistics.add_subplot(self._stat_Nrows, 3, 3, sharex=ax)
            ax.plot(self.t, self.linesearch_step_size)
            ax.set_title("Line search step size")
            ax.set_yscale("log")

            column = 1
            if self._OptimTime:
                ax = self._fig_statistics.add_subplot(self._stat_Nrows, self._stat_Ncols, self._stat_Ncols+column, sharex=ax)
                ax.plot(self.t, self.T)
                ax.set_title("Prediction Horizon")
                column += 1

            if self._OptimParam:
                ax = self._fig_statistics.add_subplot(self._stat_Nrows, self._stat_Ncols, self._stat_Ncols+column, sharex=ax)
                ax.plot(self.t, self.p)
                ax.set_title("Parameter")
                column += 1

            if self._constraints or self._terminal_constraints:
                ax = self._fig_statistics.add_subplot(self._stat_Nrows, self._stat_Ncols, self._stat_Ncols+column, sharex=ax)
                ax.plot(self.t, self.Nconstr)
                ax.set_title("Norm of constraint violation over horizon")
                column += 1

                ax = self._fig_statistics.add_subplot(self._stat_Nrows, self._stat_Ncols, self._stat_Ncols+column, sharex=ax)
                ax.plot(self.t, self.Npen)
                ax.set_title("Norm of penalty parameters over horizon")

            if self._terminal_constraints:
                ax = self._fig_statistics.add_subplot(self._stat_Nrows, 3, 7, sharex=ax)
                ax.plot(self.t, self.constrT)
                ax.set_title("Terminal Constraints")

                ax = self._fig_statistics.add_subplot(self._stat_Nrows, 3, 8, sharex=ax)
                ax.plot(self.t, self.multT)
                ax.set_title("Terminal Lagrangian multipliers")

                ax = self._fig_statistics.add_subplot(self._stat_Nrows, 3, 9, sharex=ax)
                ax.plot(self.t, self.penT)
                ax.set_title("Terminal penalty parameters")
            
            ax.set_xlim((self.t[0], self.t[-1]))

            axs = self._fig_statistics.axes
            for ax in axs:
                ax.grid()

        if self._fig_results is not None:
            self._fig_results.suptitle("Results")

            result_axs = self._fig_results.axes

            for ax in result_axs:
                ax.clear()
                ax.grid()

            result_axs[0].plot(self.t, self.x)
            result_axs[0].set_title("States")

            result_axs[1].plot(self.t, self.adj)
            result_axs[1].set_title("Adjoint states")

            result_axs[2].plot(self.t, self.u)
            result_axs[2].set_title("Controls")

            if self._constraints:
                result_axs[3].plot(self.t, self.constr)
                result_axs[3].set_title("Constraints")

                result_axs[4].plot(self.t, self.mult)
                result_axs[4].set_title("Lagrange multipliers")

                result_axs[5].plot(self.t, self.pen)
                result_axs[5].set_title("Penalty parameters")
            
            result_axs[0].set_xlim(self.t[0], self.t[-1])

        plt.draw()
        plt.pause(0.0001)
