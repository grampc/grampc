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
user.param.x0    = [0.7,0.7,-1.0472,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
user.param.xdes  = [0.2,0.25,1.0472,0.0,0.0,0.0,0.0,0.0,0.0,0.0];

% Initial values, setpoints and limits of the inputs
user.param.u0    = [0.0,0.0,0.0];
user.param.udes  = [0.0,0.0,0.0];
user.param.umax  = [2.0,2.0,2.0];
user.param.umin  = [-2.0,-2.0,-2.0];

% Time variables
user.param.Thor  = 1.5;           % Prediction horizon

user.param.dt    = 0.002;         % Sampling time

%% Option definition
% Default parameters

% optional settings for a better performance
% user.opt.LineSearchType = 'adaptive';
% user.opt.LineSearchMax = 1e1;
% user.opt.LineSearchInit = 1e-1;

%% Userparameter definition 
% e.g. system parameters or weights for the cost function
pCost  = [1.0,1.0,1.0,1.0,1.0,0.1,0.1,0.1,0.1,1.0,...   % Q
	  0.01,0.01,0.01,...                                % R
	  1.0,1.0,1.0,1.0,1.0,0.1,0.1,0.1,0.1,1.0];         % P
  
userparam = pCost;

%% Grampc initialization
grampc = CmexFiles.grampc_init_Cmex(userparam);

%% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc,user);

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
