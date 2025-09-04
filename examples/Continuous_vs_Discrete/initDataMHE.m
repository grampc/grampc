function [grampc, yfct] = initDataMHE(use_discrete, dt)
% This function initializes a grampc struct in MATLAB and sets parameters 
% and options. Define all options and parameter for the use of GRAMPC in
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

%% Measured output
yfct = @(x) x(1,:);

%% Parameter definition for MHE
% Initial values and limits of the parameters
user.param.p0    = [0, 0];

% Horizon parameterization
user.param.dt    = dt;  % Sampling time
user.opt.Nhor    = 6;   % Number of sampling steps
% compute time horizon consistent to dt and Nhor
user.param.Thor  = (user.opt.Nhor - 1) * user.param.dt;    

% weights
P = 1;  % terminal state weight
Q = 1;  % integral state weight

%% Option definition for MHE
user.opt.MaxGradIter = 10;

% Switch from control to parameter optimization
user.opt.ShiftControl = 'off';
user.opt.OptimControl = 'off';
user.opt.OptimParam = 'on';

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
end

%% Userparameter definition for MHE
y = nan * zeros(1, user.opt.Nhor);
userparam = [P, Q, y];

%% Grampc initialization
if use_discrete
    grampc = CmexFiles.grampc_init_Cmex_MHE_discrete(userparam);
else
    grampc = CmexFiles.grampc_init_Cmex_MHE_continuous(userparam);
end

%% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc, user);
