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
user.param.x0    = [ 0.1, 0.01];
user.param.xdes  = [-0.2, 0.0];

% Initial values, setpoints and limits of the inputs
user.param.u0    = 0;
user.param.udes  = 0;
user.param.umax  =  0.0524;
user.param.umin  = -0.0524;

% Time variables
user.param.Thor  = 0.3;         % Prediction horizon

user.param.dt    = 0.01;        % Sampling time
user.param.t0    = 0.0;         % time at the current sampling step

%% Option definition
% Basic algorithmic options
user.opt.Nhor        = 10;          % Number of steps for the system integration
user.opt.MaxMultIter = 3;           % Maximum number of augmented Lagrangian iterations

% Do augmented Lagrangian update every MPC step
user.opt.AugLagUpdateGradientRelTol = 1e0;

% Constraints thresholds
user.opt.ConstraintsAbsTol = 1e-3*[1 1 1 1];

% optional settings for a better performance
% user.opt.LineSearchInit = 1e-3;
% user.opt.LineSearchExpAutoFallback = 'off';
% user.opt.PenaltyMin = 1000; % Comment line 75 (grampc = CmexFiles.grampc_estim_penmin_Cmex(grampc,1);) to use this option value

%% User parameter definition 
% e.g. system parameters or weights for the cost function
userparam = [100, 10, 180, 100, 10, -0.2, 0.2, -0.1, 0.1];

%% Grampc initialization
grampc = CmexFiles.grampc_init_Cmex(userparam);

%% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc,user);

%% Estimate and set PenaltyMin (optional)
[grampc, ~] = CmexFiles.grampc_estim_penmin_Cmex(grampc,1);

%% Simulation time
if nargout>1
    Tsim = 8;
end

%% Simulink data
% Converting the grampc struct to a code generation compatible format for
% the use in simulink
if nargout>2
    grampc_sdata = grampc_generate_sdata(grampc);
end
