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
user.param.x0    = [2.02e3,1.07e3,100.0,97.1];
user.param.xdes  = [1.37e3,0.95e3,110.0,108.6];

% Initial values, setpoints and limits of the inputs
user.param.u0    = [5.0,-2540.0];
user.param.udes  = [5.0,-1190.0];
user.param.umax  = [35.0,0.0];
user.param.umin  = [3.0,-9000.0];

% Time variables
user.param.Thor  = 1200.0/3600;     % Prediction horizon
user.param.dt    = 1.0/3600;        % Sampling time

%% Option definition
% Basic algorithmic options
user.opt.Nhor        = 40;          % Number of steps for the system integration
user.opt.MaxGradIter = 3;           % Maximum number of gradient iterations

% Line search
 user.opt.LineSearchMax = 1e-6;

% Scaling values for the states, inputs and parameters
user.opt.ScaleProblem = 'on';
user.opt.xScale  = [500.0,500.0,50.0,50.0];
user.opt.xOffset = [500.0,500.0,50.0,50.0];
user.opt.uScale  = [16.0,4500.0];
user.opt.uOffset = [19.0,-4500.0];

% optional settings for a better performance
% user.opt.LineSearchMax = 1e-2;
% user.opt.LineSearchInit = 1e-6;
% user.opt.LineSearchType ='adaptive';
% user.opt.LineSearchExpAutoFallback = 'off';

%% Userparameter definition 
% e.g. system parameters or weights for the cost function
pCost  = [0.2,1.0,0.5,0.2,...
  0.2,1.0,0.5,0.2,...
  0.5,5.0e-3];

pSys   = [1.287e12,...
  9.043e6,...
  9758.3,...
  8560.0,...
  30.828,...
  86.688,...
  0.1,...
  3.522e-4,...
  104.9,...
  5.1e3,...
  4.2,...
  -11.0,...
  -41.85,...
  1];

userparam = [pSys, pCost];

%% Grampc initialization
grampc = CmexFiles.grampc_init_Cmex(userparam);

%% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc,user);

%% Simulation time
if nargout>1
    Tsim = 1.6;
end

%% Simulink data
% Converting the grampc struct to a code generation compatible format for
% the use in simulink
if nargout>2
    grampc_sdata = grampc_generate_sdata(grampc);
end
