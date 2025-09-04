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
user.param.x0    = [0.5, 0.5, 0, 0, 0, 0];
user.param.xdes  = [0, 0, 0, 0, 0, 0];

% Initial values, setpoints and limits of the inputs
user.param.u0    = [0, 0];
user.param.udes  = [0, 0];
user.param.umax  = 3 * ones(1, 2);
user.param.umin  = -1 * ones(1, 2);

% Time variables
user.param.dt    = 0.1; % Sampling time
user.opt.Nhor    = 31;  % Number of prediction steps
% compute time horizon consistent to dt and Nhor
user.param.Thor  = (user.opt.Nhor - 1) * user.param.dt;    

%% Option definition
% Basic algorithmic options
user.opt.MaxGradIter = 30;       % Maximum number of gradient iterations

user.opt.Integrator = 'discrete'; % discrete-time dynamics in probfct
user.opt.IntegratorCost = 'discrete';

% Constraints tolerances
user.opt.ConstraintsAbsTol = 1e-3*[1 1];

%% Userparameter definition 
% e.g. system parameters or weights for the cost function
userparam = [0.44, 0.6,...                         % psys
             100, 100, 10, 10, 400, 200,...        % P
          	 0.001, 0.001,...                      % R
             100, 100, 10, 10, 400, 200];          % Q
% scale integral costs by sampling time for better comparison with
% continuous-time Helicopter problem
userparam(9:end) = userparam(9:end) * user.param.dt;

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
