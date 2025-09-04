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

%% Userparameter definition 
% e.g. system parameters or weights for the cost function
pCost = [1.0,2.0,2.0,1.0,1.0,4.0,...
    0.05,0.05];
pSys = [0.2 1.25 0.3];
userparam = [pCost pSys];

%% Grampc initialization
grampc = CmexFiles.grampc_init_Cmex(userparam);


%% Read configuration from external config file
configFile = 'config_CRANE_2D.cfg';
grampc = CmexFiles.grampc_get_config_from_file_Cmex(grampc, configFile);

%% Update grampc struct while ensuring correct data types
% optional settings for a better performance
% user.opt.LineSearchMax = 1e1;
% user.opt.LineSearchInit = 1e-1;
% user.opt.LineSearchExpAutoFallback = 'off';
% user.opt.PenaltyMin = 4e1; % Comment line 44 (grampc = CmexFiles.grampc_estim_penmin_Cmex(grampc,1);) to use this option value
% grampc = grampc_update_struct_grampc(grampc,user);


%% Estimate and set PenaltyMin
[grampc, ~] = CmexFiles.grampc_estim_penmin_Cmex(grampc,1);

%% Simulation time
if nargout>1
    Tsim = 12.5;
end

%% Simulink data
% Converting the grampc struct to a code generation compatible format for
% the use in simulink
if nargout>2
    grampc_sdata = grampc_generate_sdata(grampc);
end
