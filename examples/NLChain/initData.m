function [grampc,Tsim,grampc_sdata] = initData(N)
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


if nargin ==0 || isempty(N)
    N = 4;
    disp('No chain length was given. Set N to 4');
end

xdim = 3*(2*N-1);

%% Parameter definition
% Initial values and setpoints of the states
x0     = zeros(1,xdim);
for kk = 0:(N-1)
    x0(kk * 3 + 1) = 7.5 * (kk +1) / N;
end
user.param.x0    = x0;
user.param.xdes  = zeros(1,xdim);

% Initial values, setpoints and limits of the inputs
user.param.u0    = [0.0,0.0,0.0];
user.param.udes  = [0.0,0.0,0.0];
user.param.umax  = [1.0,1.0,1.0];
user.param.umin  = [-1.0,-1.0,-1.0];

% Time variables
user.param.Thor  = 8.0;         % Prediction horizon

user.param.dt    = 0.1;         % Sampling time
user.param.t0    = 0.0;         % time at the current sampling step

%% Option definition
% Cost integration
user.opt.TerminalCost = 'off';

% Basic algorithmic options
user.opt.Nhor = 50;

% System integration
user.opt.Integrator = 'ruku45';
user.opt.IntegratorRelTol = 1e-2;
user.opt.IntegratorAbsTol = 1e-3;

% optional settings for a better performance
% user.opt.IntegratorRelTol = 1e-4;
% user.opt.IntegratorAbsTol = 1e-5;

%% Userparameter definition
% e.g. system parameters or weights for the cost function
pCost = [25, 2.5, 0.1, 10];
userparam = pCost;

%% Grampc initialization
grampc = CmexFiles.grampc_init_Cmex(userparam);

%% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc,user);

%% Simulation time
if nargout>1
    Tsim = 50;
end

%% Simulink data
% Converting the grampc struct to a code generation compatible format for
% the use in simulink
if nargout>2
    grampc_sdata = grampc_generate_sdata(grampc);
end
