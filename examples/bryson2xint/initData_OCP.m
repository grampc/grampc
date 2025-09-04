function [grampc,Tsim,grampc_sdata] = initData_OCP()
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
user.param.x0    = [0;1];

% Initial values, setpoints and limits of the inputs
user.param.u0    = 0;

% Time variables
user.param.Thor  = 1.0;        % Prediction horizon
user.param.Tmax  = 10;         
user.param.Tmin  = 0.1;          

user.param.dt    = 0.001;        % Sampling time
user.param.t0    = 0;          % time at the current sampling step

%% Option definition
% Basic algorithmic options
user.opt.MaxGradIter = 250;    % Maximum number of gradient iterations
user.opt.MaxMultIter = 1;      % Maximum number of augmented Lagrangian iterations
user.opt.Nhor = 50;
user.opt.ShiftControl = 'off';

% Line search
user.opt.LineSearchInit = 5e-7;

% Penalties
user.opt.PenaltyIncreaseFactor = 1.5;
user.opt.PenaltyIncreaseThreshold = 0.75;
user.opt.PenaltyDecreaseFactor = 1.0;
user.opt.PenaltyMin = 2e4;
user.opt.PenaltyMax = 1e7;

% Convergence test
user.opt.ConvergenceCheck = 'on';
user.opt.ConstraintsAbsTol = 1e-6*ones(3,1);
user.opt.ConvergenceGradientRelTol = 1e-8;

%% Userparameter definition 
% e.g. system parameters or weights for the cost function
userparam = [];

%% Grampc initialization
grampc = CmexFiles.grampc_init_Cmex(userparam);

%% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc,user);

%% Estimate and set PenaltyMin
[grampc, ~] = CmexFiles.grampc_estim_penmin_Cmex(grampc,1);

%% Simulation time
if nargout>1
    Tsim = 2;
end

%% Simulink data
% Converting the grampc struct to a code generation compatible format for
% the use in simulink
if nargout>2
    grampc_sdata = grampc_generate_sdata(grampc);
end
