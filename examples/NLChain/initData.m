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
% Basic algorithmic options
user.opt.MaxGradIter = 5;           % Maximum number of gradient iterations

% Cost integration
user.opt.TerminalCost = 'off';

% System integration
user.opt.Integrator = 'ruku45';
user.opt.IntegratorRelTol = 1e-3;
user.opt.IntegratorAbsTol = 1e-4;

% optional settings for a better performance
% user.opt.LineSearchType ='adaptive';
% user.opt.LineSearchMax = 1e2;
% user.opt.LineSearchInit = 1e-1;
% user.opt.LineSearchExpAutoFallback = 'off';

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
