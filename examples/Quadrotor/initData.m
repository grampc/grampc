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
user.param.x0    = [ -3.0,0.0,-3.0,0.0,-3.0,0.0,0.0,0.0,0.0];
user.param.xdes  = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];

% Initial values, setpoints and limits of the inputs
user.param.u0    = [9.81,0.0,0.0,0.0];
user.param.udes  = [9.81,0.0,0.0,0.0];
user.param.umax  = [11.0,1.0,1.0,0.5];
user.param.umin  = [0.0,-1.0,-1.0,-0.5];

% Time variables
user.param.Thor  = 1.5;         % Prediction horizon
user.param.dt    = 0.002;       % Sampling time

%% Option definition
% Basic algorithmic options
user.opt.Nhor        = 15;          % Number of steps for the system integration

% optional settings for a better performance
% user.opt.LineSearchInit = 1e-2;
% user.opt.LineSearchExpAutoFallback = 'off';

%% Userparameter definition 
% e.g. system parameters or weights for the cost function
pCost = ones(1,22);
userparam = pCost;

%% Grampc initialization
grampc = CmexFiles.grampc_init_Cmex(userparam);

%% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc,user);

%% Simulation time
if nargout>1
    Tsim = 6.0;
end

%% Simulink data
% Converting the grampc struct to a code generation compatible format for
% the use in simulink
if nargout>2
    grampc_sdata = grampc_generate_sdata(grampc);
end
