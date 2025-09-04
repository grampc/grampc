function [grampc,Tsim,grampc_sdata] = initData()
% This function initializes a grampc struct in MATLAB and sets parameters 
% and options. In case of three output arguments the struct grampc_sdata for
% the use in Simulink is created as well. Define all options and parameters
% for the use of GRAMPC in MATLAB here.
%
% This file is part of GRAMPC - (https://github.com/grampc/grampc)
%
% GRAMPC -- A software framework for embedded nonlinear model predictive
% control using a gradient-based augmented Lagrangian approach
%
% Copyright 2014-2025 by Knut Graichen, Andreas Voelz, Thore Wietzke,
% Tobias Englert (<v2.3), Felix Mesmer (<v2.3), Soenke Rhein (<v2.3),
% Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
% All rights reserved.
%
% GRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
%

%% Parameter definition
% Initial values and setpoints of the states
user.param.x0    = [1.0000,    0.9970,    0.9876,    0.9717,    0.9491,    0.9195,    0.8829,    0.8390,    0.7878,    0.7294,    0.6637];
user.param.xdes  = [2.0000,    1.9929,    1.9711,    1.9342,    1.8816,    1.8130,    1.7283,    1.6273,    1.5101,    1.3768,    1.2276];

% Initial values, setpoints and limits of the inputs
user.param.u0    = -0.6935;
user.param.udes  = -1.5695;
user.param.umax  = 2.0;
user.param.umin  = -2.0;

% Time variables
user.param.Thor  = 0.4;         % Prediction horizon
user.param.dt    = 0.005;       % Sampling time

%% Option definition
% Basic algorithmic options
user.opt.Nhor  = 60;

% System integration
user.opt.Integrator = 'rodas';
user.opt.IntegratorRelTol = 1e-3;
user.opt.IntegratorAbsTol = 1e-4;
IFCN = 0;		% 0 --> right hand side independent of time t  
IDFX = 0;		% 0 --> DF/Dt is numerically computed 
IJAC = 1;		% 1(0) -> analytical (numerical) jacobian (partial derivatives of right hand side w.r.t. state) 
MLJAC = 1;		% no. of lower diagonals of jacobian ~	value must be between 1 and NX 
MUJAC = 1;		% no. of upper diagonals of jacobian ~	value must be between 1 and NX 
IMAS = 1;		% 1 --> mass matrix is supplied 
MLMAS = 1;		% no. of lower diagonals of mass matrix 
MUMAS = 1;		% no. of upper diagonals of mass matrix 
user.opt.FlagsRodas = [IFCN, IDFX, IJAC, IMAS, MLJAC, MUJAC, MLMAS, MUMAS];

% Line search
user.opt.LineSearchMax = 2.0;

% optional settings for a better performance
% user.opt.LineSearchInit = 2e-1;
% user.opt.LineSearchType ='adaptive';

%% Userparameter definition 
% e.g. system parameters or weights for the cost function
pCost = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,...
	  10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,...
	  1];
userparam = pCost;

%% Grampc initialization
grampc = CmexFiles.grampc_init_Cmex(userparam);

%% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc,user);

%% Simulation time
if nargout>1
    Tsim = 3;
end

%% Simulink data
% Converting the grampc struct to a code generation compatible format for
% the use in simulink
if nargout>2
    grampc_sdata = grampc_generate_sdata(grampc);
end
