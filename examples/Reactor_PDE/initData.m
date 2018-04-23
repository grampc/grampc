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
user.param.x0    = [1.0000,    0.9970,    0.9876,    0.9717,    0.9491,    0.9195,    0.8829,    0.8390,    0.7878,    0.7294,    0.6637];
user.param.xdes  = [2.0000,    1.9929,    1.9711,    1.9342,    1.8816,    1.8130,    1.7283,    1.6273,    1.5101,    1.3768,    1.2276];

% Initial values, setpoints and limits of the inputs
user.param.u0    = -0.6935;
user.param.udes  = -1.5695;
user.param.umax  = 2.0;
user.param.umin  = -2.0;

% Time variables
user.param.Thor  = 0.4;         % Prediction horizon
user.param.dt    = 0.005;       % Sampling time

%% Option definition
% Basic algorithmic options
user.opt.Nhor  = 60;

% System integration
user.opt.Integrator = 'rodas';
user.opt.IntegratorRelTol = 1e-3;
user.opt.IntegratorAbsTol = 1e-4;
IFCN = 0;		% 0 --> right hand side independent of time t  
IDFX = 0;		% 0 --> DF/Dt is numerically computed 
IJAC = 1;		% 1(0) -> analytical (numerical) jacobian (partial derivatives of right hand side w.r.t. state) 
MLJAC = 1;		% no. of lower diagonals of jacobian ~	value must be between 1 and NX 
MUJAC = 1;		% no. of upper diagonals of jacobian ~	value must be between 1 and NX 
IMAS = 1;		% 1 --> mass matrix is supplied 
MLMAS = 1;		% no. of lower diagonals of mass matrix 
MUMAS = 1;		% no. of upper diagonals of mass matrix 
user.opt.FlagsRodas = [IFCN, IDFX, IJAC, IMAS, MLJAC, MUJAC, MLMAS, MUMAS];

% Line search
user.opt.LineSearchMax = 2.0;

% optional settings for a better performance
% user.opt.LineSearchInit = 2e-1;
% user.opt.LineSearchType ='adaptive';

%% Userparameter definition 
% e.g. system parameters or weights for the cost function
pCost = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,...
	  10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,...
	  1];
userparam = pCost;

%% Grampc initialization
grampc = CmexFiles.grampc_init_Cmex(userparam);

%% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc,user);

%% Simulation time
if nargout>1
    Tsim = 3;
end

%% Simulink data
% Converting the grampc struct to a code generation compatible format for
% the use in simulink
if nargout>2
    grampc_sdata = grampc_generate_sdata(grampc);
end
