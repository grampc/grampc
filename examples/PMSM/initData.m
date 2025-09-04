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
user.param.x0    = [0.0,0.0,0.0,0.0];
user.param.xdes  = [0.0,9.5,0,0.0]; 

% Initial values, setpoints and limits of the inputs
user.param.u0    = [0.0,0.0];
user.param.udes  = [0.0,0.0];
user.param.umax  = [1,1]*560/sqrt(3);
user.param.umin  = [-1,-1]*560/sqrt(3);

% Time variables
user.param.Thor  = 0.005;           % Prediction horizon
user.param.dt    = 0.000125;        % Sampling time
user.param.t0    = 0.0;             % time at the current sampling step

%% Option definition
% Basic algorithmic options
user.opt.Nhor        = 11;          % Number of steps for the system integration
user.opt.MaxGradIter = 3;           % Maximum number of gradient iterations
user.opt.MaxMultIter = 3;           % Maximum number of augmented Lagrangian iterations

% Constraint tolerances
user.opt.ConstraintsAbsTol = 1e-3*[1 1];
% user.opt.ConstraintsAbsTol = [0.1 104.5];

% optional settings for a better performance
% user.opt.PenaltyMin = 10000; % Comment line 83 (grampc = CmexFiles.grampc_estim_penmin_Cmex(grampc,1);) to use this option value
% user.opt.AugLagUpdateGradientRelTol = 1e0;
% user.opt.LineSearchMax = 1e1;

% System integration
user.opt.TerminalCost = 'off';

%% Userparameter definition 
% e.g. system parameters or weights for the cost function
pSys = [3.5,0.0175,0.0175,0.17,3,9e-4,4e-4,0,560^2/3,10^2];
pCost = [8.0,200.0,0.0,0.0,...
    0.001,0.001];
userparam = [pSys pCost];

% Scaling of the constraints (if not done by hand in the probfct_PMSM.c)
% user.opt.ScaleProblem = 'on';
% user.opt.cScale = userparam(9:10);

%% Grampc initialization
grampc = CmexFiles.grampc_init_Cmex(userparam);

%% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc,user);

%% Estimate and set PenaltyMin
[grampc, ~] = CmexFiles.grampc_estim_penmin_Cmex(grampc,1);

%% Simulation time
if nargout>1
    Tsim = 0.10;
end

%% Simulink data
% Converting the grampc struct to a code generation compatible format for
% the use in simulink
if nargout>2
    grampc_sdata = grampc_generate_sdata(grampc);
end
