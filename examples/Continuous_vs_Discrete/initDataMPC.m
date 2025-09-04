function [grampc,Tsim] = initDataMPC(use_discrete)
% This function initializes a grampc struct in MATLAB and sets parameters 
% and options. Define all options and parameters for the use of GRAMPC in
% MATLAB here.
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
user.param.x0    = [-1, -1];
user.param.xdes  = [0, 0];

% Initial values, setpoints and limits of the inputs
user.param.u0    = 0;
user.param.udes  = 0;
user.param.umax  = 1.0;
user.param.umin  = -1.0;

% Horizon parameterization
user.param.dt    = 0.2;     % Sampling time
user.opt.Nhor    = 6;       % Number of sampling steps
% compute time horizon consistent to dt and Nhor
user.param.Thor  = (user.opt.Nhor - 1) * user.param.dt;         

% weights
P = [1.0, 1.0]; % terminal state weight
Q = [1.0, 1.0]; % integral state weight
R = 0.01;       % integral control weight

%% Option definition
user.opt.MaxGradIter = 20;

if use_discrete
    % formulation as discrete-time system
    user.opt.Integrator = 'discrete';
    user.opt.IntegratorCost = 'discrete';
else
    % formulation as continuous-time system
    user.opt.Integrator = 'erk1';
    user.opt.IntegratorCost = 'trapezoidal';
    % scale integral cost to get same cost as in discrete version
    Q = Q / user.param.dt;
    R = R / user.param.dt;
end

%% Userparameter definition 
% e.g. system parameters or weights for the cost function
userparam = [P, Q, R];

%% Grampc initialization
if use_discrete
    grampc = CmexFiles.grampc_init_Cmex_MPC_discrete(userparam);
else
    grampc = CmexFiles.grampc_init_Cmex_MPC_continuous(userparam);
end

%% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc,user);

%% Simulation time
if nargout>1
    Tsim = 10;
end