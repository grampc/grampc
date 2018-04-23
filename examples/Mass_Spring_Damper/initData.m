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
m = 1;
c = 1;
d = 0.2;

NN = 5;
NX = 2*NN;
NU = 2;

% Initial values and setpoints of the states
user.param.x0    = zeros(NX,1);
user.param.xdes  = zeros(NX,1);
user.param.x0(1) = 1.0;
user.param.x0(NN) = 1.0;

% Initial values, setpoints and limits of the inputs
user.param.u0    = [0.0,0.0]';
user.param.udes  = [0.0,0.0]';
user.param.umax  = [1,1];
user.param.umin  = [-1,-1];

% Time variables
user.param.Thor  = 10.0;         % Prediction horizon
user.param.dt    = 0.005;        % Sampling time

%% Option definition
% Basic algorithmic options
user.opt.MaxGradIter     = 5;    % Maximum number of gradient iterations

% optional settings for a better performance
% user.opt.LineSearchType ='adaptive';
% user.opt.LineSearchMax = 1e1;
% user.opt.LineSearchInit = 1e-1;
% user.opt.LineSearchExpAutoFallback = 'off';

%% Userparameter definitin 
% e.g. system parameters or weights for the cost function
NpCost      = 2*NX+NU;
pCost       = ones(1,NpCost);
pCost(NX+1) = 0.01;
pCost(NX+2) = 0.01;  
pSys  = [m,c,d];

userparam = [pSys, pCost];

%% Grampc initialization
grampc = CmexFiles.grampc_init_Cmex(userparam);

%% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc,user);

%% Simulation time
if nargout>1
    Tsim = 12.0;
end

%% Simulink data
% Converting the grampc struct to a code generation compatible format for
% the use in simulink
if nargout>2
    grampc_sdata = grampc_generate_sdata(grampc);
end
