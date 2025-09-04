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
user.param.x0    = [0.0, 1.5, 0.0, 0.0, 17.0];
user.param.xdes  = [0.0, 1.5, 0.0, 0.0, 0.0];

% Initial values, setpoints and limits of the inputs
user.param.u0    = [0.0,0.0];
user.param.udes  = [0.0,0.0];
user.param.umax  = [ 0.5, 5.0];
user.param.umin  = [-0.5,-5.0];

% Time variables
user.param.Thor  = 1;          % Prediction horizon

user.param.dt    = 0.01;       % Sampling time
user.param.t0    = 0.0;        % time at the current sampling step

%% Option definition
% Basic algorithmic options
user.opt.Nhor    = 20;          % Number of steps for the system integration
user.opt.MaxGradIter = 5;

% Constraints tolerances
user.opt.ConstraintsAbsTol = 1e-2;

% optional settings for a better performance
% user.opt.LineSearchType ='adaptive';
% user.opt.LineSearchMax = 1e2;
% user.opt.PenaltyMin = 2e2; % Comment line 77 (grampc = CmexFiles.grampc_estim_penmin_Cmex(grampc,1);) to use this option value
% user.opt.PenaltyIncreaseFactor = 1.1;

%% Userparameter definition 
% e.g. system parameters or weights for the cost function
pCost = [0.0, 1.0, 1.0, 1.0, 100.0,...
    0.0, 1.0, 1.0, 1.0, 100.0,...
    0.01, 0.01];

userparam = [3, 50, pCost];

%% Grampc initialization
grampc = CmexFiles.grampc_init_Cmex(userparam);

%% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc,user);

%% Estimate and set PenaltyMin
[grampc, ~] = CmexFiles.grampc_estim_penmin_Cmex(grampc,1);

%% Simulation time
if nargout>1
    Tsim = 4;
end

%% Simulink data
% Converting the grampc struct to a code generation compatible format for
% the use in simulink
if nargout>2
    grampc_sdata = grampc_generate_sdata(grampc);
end
