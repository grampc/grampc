function [grampc,Tsim,grampc_sdata] = initData_OCP()
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
user.param.x0    = [ pi/2, -pi/2, 0, -pi/2,  pi/2, 0];
user.param.xdes  = [-pi/2,  pi/2, 0,  pi/2, -pi/2, 0];

% Initial values, setpoints and limits of the inputs
user.param.u0    = [0, 0, 0, 0, 0, 0];
user.param.udes  = [0, 0, 0, 0, 0, 0];
user.param.umax  = 1 * ones(1, 6);
user.param.umin  = -1 * ones(1, 6);

% Time variables
user.param.Thor  = 10;         % Prediction horizon
user.param.dt    = 0.1;        % Sampling time

%% Option definition
% Basic algorithmic options
user.opt.Nhor        = 20;          % Number of steps for the system integration
user.opt.MaxGradIter = 1000;        % Maximum number of gradient iterations
user.opt.MaxMultIter = 1;           % Maximum number of augmented Lagrangian iterations
user.opt.ShiftControl = 'off';

% Line search
user.opt.LineSearchMax = 2;
user.opt.LineSearchExpAutoFallback = 'off';
user.opt.OptimTime = 'off';

% Type of considered constraints
user.opt.AugLagUpdateGradientRelTol = 1;

% Convergence test
user.opt.ConstraintsAbsTol = 1e-4 * [ones(1,3), ones(1,6)];
user.opt.ConvergenceCheck = 'off';
user.opt.ConvergenceGradientRelTol = 1e-6;

% Penalties
user.opt.PenaltyMax = 1e4;
user.opt.PenaltyMin = 50;
user.opt.PenaltyIncreaseFactor = 1.1;
user.opt.PenaltyDecreaseFactor = 1.0;

%% Userparameter definition 
% e.g. system parameters or weights for the cost function
weights = [0, 0, 0, 0, 1, 1];   % weights for cost function
lengths = [0.5, 0.3, 0.2, ...   % length left arm
           0.5, 0.3, 0.2];      % length right arm
offsets = [0, 0, 0, ...         % offset left arm
           1, 0, pi];           % offset right arm
userparam = [weights, lengths, offsets, user.param.xdes];

%% Grampc initialization
grampc = CmexFiles.grampc_init_Cmex(userparam);

%% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc,user);

%% Estimate and set PenaltyMin
% [grampc, ~] = CmexFiles.grampc_estim_penmin_Cmex(grampc,1);

%% Simulation time
if nargout>1
    Tsim = 400;
end

%% Simulink data
% Converting the grampc struct to a code generation compatible format for
% the use in simulink
if nargout>2
    grampc_sdata = grampc_generate_sdata(grampc);
end
