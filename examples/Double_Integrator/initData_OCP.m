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
user.param.x0    = [-1, -1];
user.param.xdes  = [0, 0];

% Initial values, setpoints and limits of the inputs
user.param.u0    = 0;
user.param.udes  = 0;
user.param.umax  = 1.0;
user.param.umin  = -1.0;

% Time variables
user.param.Thor  = 5.25;        % Prediction horizon
% user.param.Thor  = 4.0;        % Prediction horizon
user.param.Tmax  = 10;         
user.param.Tmin  = 1.0;          

user.param.dt    = 0.01;       % Sampling time
user.param.t0    = 0;          % time at the current sampling step

%% Option definition
% Basic algorithmic options
user.opt.MaxGradIter = 200;    % Maximum number of gradient iterations
user.opt.MaxMultIter = 1;      % Maximum number of augmented Lagrangian iterations
user.opt.Nhor = 50;
user.opt.ShiftControl = 'off';

% user.opt.InequalityConstraints = 'off';

% System integration
user.opt.Integrator = 'euler';

% Line search
user.opt.LineSearchMax = 1e2;

% Input and or parameter optimization 
user.opt.OptimTime = 'on';
user.opt.OptimTimeLineSearchFactor = 1.75;

% Penalties
user.opt.PenaltyMin = 1e1;
user.opt.PenaltyIncreaseFactor = 1.25;
user.opt.PenaltyDecreaseFactor = 1.0;

% Convergence test
user.opt.ConvergenceCheck = 'on';
user.opt.ConstraintsAbsTol = 1e-6 * ones(1,3);
user.opt.ConvergenceGradientRelTol = 1e-9;

%% Userparameter definition 
% e.g. system parameters or weights for the cost function
userparam = [0.1, ...    % integral cost
             1];        % horizon

%% Grampc initialization
grampc = CmexFiles.grampc_init_Cmex(userparam);

%% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc,user);

%% Simulation time
if nargout>1
    Tsim = 7;
end

%% Simulink data
% Converting the grampc struct to a code generation compatible format for
% the use in simulink
if nargout>2
    grampc_sdata = grampc_generate_sdata(grampc);
end
