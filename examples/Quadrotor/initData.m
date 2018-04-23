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
% Copyright (C) 2014-2018 by Tobias Englert, Knut Graichen, Felix Mesmer,
% Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0). 
% Developed at the Institute of Measurement, Control, and Microtechnology,
% Ulm University. All rights reserved.
%
% GRAMPC is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as 
% published by the Free Software Foundation, either version 3 of 
% the License, or (at your option) any later version.
%
% GRAMPC is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public 
% License along with GRAMPC. If not, see <http://www.gnu.org/licenses/>.
%

%% Parameter definition
% Initial values and setpoints of the states
user.param.x0    = [ -3.0,0.0,-3.0,0.0,-3.0,0.0,0.0,0.0,0.0];
user.param.xdes  = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];

% Initial values, setpoints and limits of the inputs
user.param.u0    = [9.81,0.0,0.0,0.0];
user.param.udes  = [9.81,0.0,0.0,0.0];
user.param.umax  = [11.0,1.0,1.0,0.5];
user.param.umin  = [0.0,-1.0,-1.0,-0.5];

% Time variables
user.param.Thor  = 1.5;         % Prediction horizon
user.param.dt    = 0.002;       % Sampling time

%% Option definition
% Basic algorithmic options
user.opt.Nhor        = 15;          % Number of steps for the system integration

% optional settings for a better performance
% user.opt.LineSearchInit = 1e-2;
% user.opt.LineSearchExpAutoFallback = 'off';

%% Userparameter definition 
% e.g. system parameters or weights for the cost function
pCost = ones(1,22);
userparam = pCost;

%% Grampc initialization
grampc = CmexFiles.grampc_init_Cmex(userparam);

%% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc,user);

%% Simulation time
if nargout>1
    Tsim = 6.0;
end

%% Simulink data
% Converting the grampc struct to a code generation compatible format for
% the use in simulink
if nargout>2
    grampc_sdata = grampc_generate_sdata(grampc);
end
