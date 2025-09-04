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

class GrampcParam:
    """
    GRAMPC parameters struct, read only
    """
    def __init__(self) -> None: ...
    Nc: int
    Ng: int
    NgT: int
    Nh: int
    NhT: int
    Np: int
    Nu: int
    Nx: int
    Thor: float
    Tmax: float
    Tmin: float
    dt: float
    p0: np.ndarray
    pmax: np.ndarray
    pmin: np.ndarray
    t0: float
    u0: np.ndarray
    udes: np.ndarray
    umax: np.ndarray
    umin: np.ndarray
    x0: np.ndarray
    xdes: np.ndarray

class GrampcOpt:
    """
    GRAMPC options struct, read only
    """
    def __init__(self) -> None: ...
    AugLagUpdateGradientRelTol: float
    ConstraintsAbsTol: np.ndarray
    ConstraintsHandling: int
    ConvergenceCheck: int
    ConvergenceGradientRelTol: float
    EqualityConstraints: int
    FlagsRodas: list[int]
    InequalityConstraints: int
    IntegralCost: int
    Integrator: int
    IntegratorAbsTol: float
    IntegratorCost: int
    IntegratorMaxSteps: int
    IntegratorMinStepSize: float
    IntegratorRelTol: float
    JScale: float
    LineSearchAdaptAbsTol: float
    LineSearchAdaptFactor: float
    LineSearchExpAutoFallback: int
    LineSearchInit: float
    LineSearchIntervalFactor: float
    LineSearchIntervalTol: float
    LineSearchMax: float
    LineSearchMin: float
    LineSearchType: int
    MaxGradIter: int
    MaxMultIter: int
    MultiplierDampingFactor: float
    MultiplierMax: float
    Nhor: int
    OptimControl: int
    OptimParam: int
    OptimParamLineSearchFactor: float
    OptimTime: int
    OptimTimeLineSearchFactor: float
    PenaltyDecreaseFactor: float
    PenaltyIncreaseFactor: float
    PenaltyIncreaseThreshold: float
    PenaltyMax: float
    PenaltyMin: float
    ScaleProblem: int
    ShiftControl: int
    TOffset: float
    TScale: float
    TerminalCost: int
    TerminalEqualityConstraints: int
    TerminalInequalityConstraints: int
    TimeDiscretization: int
    cScale: np.ndarray
    pOffset: np.ndarray
    pScale: np.ndarray
    uOffset: np.ndarray
    uScale: np.ndarray
    xOffset: np.ndarray
    xScale: np.ndarray

class GrampcRWS:
    """
    GRAMPC private real workspace struct with all internal variables, read only
    """
    def __init__(self) -> None: ...
    T: float
    Tprev: float
    adj: np.ndarray
    cfct: np.ndarray
    cfctAbsTol: np.ndarray
    cfctprev: np.ndarray
    dcdp: np.ndarray
    dcdt: float
    dcdu: np.ndarray
    dcdx: np.ndarray
    gradT: float
    gradp: np.ndarray
    gradpprev: np.ndarray
    gradu: np.ndarray
    graduprev: np.ndarray
    iparRodas: list[int]
    iworkRodas: list[int]
    liworkRodas: int
    lrwsGeneral: int
    lsAdapt: np.ndarray
    lsExplicit: np.ndarray
    lworkRodas: int
    mult: np.ndarray
    p: np.ndarray
    pen: np.ndarray
    pls: np.ndarray
    pprev: np.ndarray
    rparRodas: np.ndarray
    rwsGeneral: np.ndarray
    rwsScale: np.ndarray
    t: np.ndarray
    tls: np.ndarray
    u: np.ndarray
    uls: np.ndarray
    uprev: np.ndarray
    workRodas: np.ndarray
    x: np.ndarray

class GrampcSol:
    """
    GRAMPC solution struct, read only
    """
    def __init__(self) -> None: ...
    J: np.ndarray
    Tnext: float
    cfct: float
    iter: list[int]
    pen: float
    pnext: np.ndarray
    status: int
    unext: np.ndarray
    xnext: np.ndarray

class GrampcBinding:

    problem: ProblemDescription
    opt: GrampcOpt
    param: GrampcParam
    rws: GrampcRWS
    sol: GrampcSol

    def __init__(self, problem: ProblemDescription) -> None: ...

    def run(self) -> float:
        """
        Calls grampc_run which executes one MPC step. Updates the values inside rws and sol.

        Returns:
            float: CPU wall clock time of one function call in milliseconds.

        Raises:
            RuntimeError: if dt or Thor don't have valid values
        """

    def get_config_from_file(self, file_name: str) -> None:
        """
        Reads a config file with the config reader from GRAMPC
        
        Args: 
            file_name (str): file name of the config file.
        """
        
    def estim_penmin(self, run_grampc: bool) -> None:
        """
        Estimates the minimal penalty parameter value. 
        
        Args: 
            run_grampc (bool): Specifies if grampc_run() shall be called.
        
        Raises:
            RuntimeError: if dt or Thor don't have valid values
        """

    def check_gradients(self, tolerance: float, step_size: float) -> None:
        """
        Checks the implemented gradients of the problem description with the values inside self.param 
        
        Args: 
            tolerance (float): tolerance for the relative error between analytical and finite differences.
            step_size (float): step-size of the finite differences.
        """

    def set_rws_u(self, u_new: np.ndarray) -> None:
        """
        Sets the data in rws.u with u_new

        Args:
            u_new (np.ndarray): New data.

        Raises:
            ValueError: if the dimensions of u_new don't match rws.u
        """

    def set_rws_multiplier(self, multiplier_new: np.ndarray) -> None:
        """
        Sets the data in rws.mult with multiplier_new

        Args:
            multiplier_new (np.ndarray): New data.

        Raises:
            ValueError: if the dimensions of multiplier_new don't match rws.mult
        """

    def set_rws_penalty(self, penalty_new: np.ndarray) -> None:
        """
        Sets the data in rws.pen with penalty_new

        Args:
            penalty_new (np.ndarray): New data.

        Raises:
            ValueError: if the dimensions of penalty_new don't match rws.pen
        """

    def ffct(self, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, param: GrampcParam) -> np.ndarray: ...
    def lfct(self, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, param: GrampcParam) -> float: ...
    def Vfct(self, T: float, x: np.ndarray, p: np.ndarray, param: GrampcParam) -> float: ...

    def gTfct(self, T: float, x: np.ndarray, p: np.ndarray, param: GrampcParam) -> np.ndarray: ...
    def gfct(self, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, param: GrampcParam) -> np.ndarray: ...
    def hTfct(self, T: float, x: np.ndarray, p: np.ndarray, param: GrampcParam) -> np.ndarray: ...
    def hfct(self, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, param: GrampcParam) -> np.ndarray: ...

    def print_opts(self) -> None:
        """
        Prints the options of the underlying grampc_opt struct.
        """

    def print_params(self) -> None:
        """
        Prints the parameters of the underlying grampc_param struct.
        """

    def print_status(self) -> None:
        """
        Prints the current status of GRAMPC.
        """

class ProblemDescription:
    """
    Provides an interface for problem descriptions. Every problem description must inherit from this class.
    """
    Nx: int
    Nu: int
    Np: int
    Ng: int
    NgT: int
    Nh: int
    NhT: int

    def __init__(self, Nx, Nu, Np, Ng, Nh, NgT, NhT) -> None: 
        "Constructor for the problem description"
    def ffct(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, param: GrampcParam) -> None: 
        "System function f(x, u, p, t) with out[Nx]. Must be implemented"
    def dfdx_vec(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, vec: np.ndarray, param: GrampcParam) -> None: 
        "Jacobi vector product df/dx^T * vec with out[Nx] and vec[Nx]. Must be implemented"
    def dfdu_vec(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, vec: np.ndarray, param: GrampcParam) -> None: 
        "Jacobi vector product df/du^T * vec with out[Nu] and vec[Nx]"
    def dfdp_vec(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, vec: np.ndarray, param: GrampcParam) -> None: 
        "Jacobi vector product df/dp^T * vec with out[Np] and vec[Nx]"

    def lfct(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, param: GrampcParam) -> None: 
        "Integral cost with out[1]"
    def dldx(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, param: GrampcParam) -> None: 
        "Integral cost derivative dl/dx with out[Nx]"
    def dldu(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, param: GrampcParam) -> None: 
        "Integral cost derivative dl/du with out[Nu]"
    def dldp(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, param: GrampcParam) -> None: 
        "Integral cost derivative dl/dp with out[Np]"

    def Vfct(self, out: np.ndarray, T: float, x: np.ndarray, p: np.ndarray, param: GrampcParam) -> None: 
        "Terminal cost with out[1]"
    def dVdx(self, out: np.ndarray, T: float, x: np.ndarray, p: np.ndarray, param: GrampcParam) -> None: 
        "Terminal cost derivative dV/dx with out[Nx]"
    def dVdp(self, out: np.ndarray, T: float, x: np.ndarray, p: np.ndarray, param: GrampcParam) -> None: 
        "Terminal cost derivative dV/dp with out[Np]"
    def dVdT(self, out: np.ndarray, T: float, x: np.ndarray, p: np.ndarray, param: GrampcParam) -> None: 
        "Terminal cost derivative dV/dT with out[1]"

    def gfct(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, param: GrampcParam) -> None: 
        "Equality constraints with out[Ng]"
    def dgdx_vec(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, vec: np.ndarray, param: GrampcParam) -> None: 
        "Jacobi vector product dg/dx^T * vec with out[Nx] and vec[Ng]"
    def dgdu_vec(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, vec: np.ndarray, param: GrampcParam) -> None: 
        "Jacobi vector product dg/du^T * vec with out[Nu] and vec[Ng]"
    def dgdp_vec(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, vec: np.ndarray, param: GrampcParam) -> None: 
        "Jacobi vector product dg/dp^T * vec with out[Np] and vec[Ng]"

    def hfct(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, param: GrampcParam) -> None: 
        "Inequality constraints with out[Nh]"
    def dhdx_vec(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, vec: np.ndarray, param: GrampcParam) -> None: 
        "Jacobi vector product dh/dx^T * vec with out[Nx] and vec[Nh]"
    def dhdu_vec(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, vec: np.ndarray, param: GrampcParam) -> None: 
        "Jacobi vector product dh/du^T * vec with out[Nu] and vec[Nh]"
    def dhdp_vec(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, vec: np.ndarray, param: GrampcParam) -> None: 
        "Jacobi vector product dh/dp^T * vec with out[Np] and vec[Nh]"

    def gTfct(self, out: np.ndarray, T: float, x: np.ndarray, u: np.ndarray) -> None: 
        "Terminal equality constraints with out[NgT]"
    def dgTdx_vec(self, out: np.ndarray, T: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, param: GrampcParam) -> None: 
        "Jacobi vector product dgT/dx^T * vec with out[Nx] and vec[NgT]"
    def dgTdp_vec(self, out: np.ndarray, T: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, param: GrampcParam) -> None: 
        "Jacobi vector product dgT/dp^T * vec with out[Np] and vec[NgT]"
    def dgTdT_vec(self, out: np.ndarray, T: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, param: GrampcParam) -> None: 
        "Jacobi vector product dgT/dT^T * vec with out[1] and vec[NgT]"

    def hTfct(self, out: np.ndarray, T: float, x: np.ndarray, u: np.ndarray) -> None: 
        "Terminal inequality constraints with out[NhT]"
    def dhTdx_vec(self, out: np.ndarray, T: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, param: GrampcParam) -> None: 
        "Jacobi vector product dhT/dx^T * vec with out[Nx] and vec[NhT]"
    def dhTdp_vec(self, out: np.ndarray, T: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, param: GrampcParam) -> None: 
        "Jacobi vector product dhT/dp^T * vec with out[Np] and vec[NhT]"
    def dhTdT_vec(self, out: np.ndarray, T: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, param: GrampcParam) -> None: 
        "Jacobi vector product dhT/dT^T * vec with out[1] and vec[NhT]"

    def dfdx(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, param: GrampcParam) -> None: 
        "System jacobian for Rodas with out[Nx * (MLJAC + MUJAC + 1)]"
    def dfdxtrans(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, param: GrampcParam) -> None: 
        "Transposed system jacobian for Rodas with out[Nx * (MLJAC + MUJAC + 1)]"
    def dfdt(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, param: GrampcParam) -> None: 
        "df/dt for Rodas with out[Nx]"
    def dHdxdt(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, vec: np.ndarray, param: GrampcParam) -> None: 
        "Time derivative of dH/dx for Rodas with out[Nx]"
    def Mfct(self, out: np.ndarray) -> None: 
        "Mass matrix for Rodas with out[Nx * (MLMAS + MUMAS + 1)]"
    def Mtrans(self, out: np.ndarray) -> None: 
        "Transposed mass matrix for Rodas with out[Nx * (MLMAS + MUMAS + 1)]"
