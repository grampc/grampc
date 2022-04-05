function [grampc,Tsim,grampc_sdata] = initData()
% This function initializes a grampc struct in MATLAB and sets parameters 
% and options. In case of three output arguments the struct grampc_sdata for
% the use in Simulink is created as well. Define all options and parameters
% for the use of GRAMPC in MATLAB here.
%
% This file is part of GRAMPC - (https://sourceforge.net/projects/grampc/)
%
% GRAMPC -- A software framework for embedded nonlinear model predictive
% control using a gradient-based augmented Lagrangian approach
%
% Copyright 2014-2019 by Tobias Englert, Knut Graichen, Felix Mesmer,
% Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
% All rights reserved.
%
% GRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
%

%% Parameter definition
% Initial values and setpoints of the states
user.param.x0    = [0.5, 0.5, 0, 0, 0, 0];
user.param.xdes  = [0, 0, 0, 0, 0, 0];

% Initial values, setpoints and limits of the inputs
user.param.u0    = [0, 0,];
user.param.udes  = [0, 0];
user.param.umax  = 3 * ones(1, 2);
user.param.umin  = -1 * ones(1, 2);

% Time variables
user.param.Thor  = 3.0;         % Prediction horizon

user.param.dt    = 0.01;        % Sampling time
user.param.t0    = 0.0;         % time at the current sampling step

%% Option definition
% Basic algorithmic options
user.opt.MaxGradIter = 3;       % Maximum number of gradient iterations

% Constraints tolerances
user.opt.ConstraintsAbsTol = 1e-3*[1 1];

% optional settings for a better performance
% user.opt.LineSearchMax = 5e-3;
% user.opt.LineSearchInit = 2e-4;
% user.opt.LineSearchExpAutoFallback = 'off';
% user.opt.PenaltyIncreaseFactor = 1.75;
% user.opt.PenaltyMin = 5e4; % Comment line 75 (grampc = CmexFiles.grampc_estim_penmin_Cmex(grampc,1);) to use this option value

%% Userparameter definition 
% e.g. system parameters or weights for the cost function
userparam = [0.44, 0.6,...                         % psys
             100, 100, 10, 10, 400, 200,...        % P
          	 0.001, 0.001,...                      % R
             100, 100, 10, 10, 400, 200];          % Q];

%% Grampc initialization
grampc = CmexFiles.grampc_init_Cmex(userparam);

%% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc,user);

%% Estimate and set PenaltyMin
[grampc, ~] = CmexFiles.grampc_estim_penmin_Cmex(grampc,1);

%% Simulation time
if nargout>1
    Tsim = 11;
end

%% Simulink data
% Converting the grampc struct to a code generation compatible format for
% the use in simulink
if nargout>2
    grampc_sdata = grampc_generate_sdata(grampc);
end
