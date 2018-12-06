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
Nx = 3;

% Initial values and setpoints of the states
user.param.x0    = [1.0, 0.0, 1.0];
user.param.xdes  = [0.0, 1.0, 1.0];

% Initial values, setpoints and limits of the inputs
user.param.u0    = [-2.0, 2.0];
user.param.udes  = [0.0, 0.0];
user.param.umax  = [ 2.0,  2.0];
user.param.umin  = [-2.0, -2.0];

% Time variables
user.param.Thor  = 0.4;          % Prediction horizon

user.param.dt    = 0.01;          % Sampling time
user.param.t0    = 0.0;         % time at the current sampling step

%% Option definition
% Basic algorithmic options
user.opt.Nhor        = 30;          % Number of steps for the system integration
user.opt.MaxGradIter = 10;           % Maximum number of gradient iterations
user.opt.MaxMultIter = 3;           % Maximum number of augmented Lagrangian iterations

% Cost integration
user.opt.IntegralCost = 'on';
user.opt.TerminalCost = 'off';

% Line search
user.opt.LineSearchMax = 1e-1;

% System integration
user.opt.Integrator = 'rodas';
user.opt.IntegratorRelTol = 1e-4;
user.opt.IntegratorAbsTol = 1e-5;
IFCN = 0;		% 0 --> right hand side independent of time t  
IDFX = 0;		% 0 --> DF/Dt is numerically computed 
IJAC = 0;		% 1(0) -> analytical (numerical) jacobian (partial derivatives of right hand side w.r.t. state) 
IMAS = 1;		% 1 --> mass matrix is supplied 
MLJAC = Nx;		% no. of lower diagonals of jacobian ~	value must be between 1 and NX 
MUJAC = Nx;		% no. of upper diagonals of jacobian ~	value must be between 1 and NX 
MLMAS = 0;		% no. of lower diagonals of mass matrix 
MUMAS = 0;		% no. of upper diagonals of mass matrix 
user.opt.FlagsRodas = [IFCN, IDFX, IJAC, IMAS, MLJAC, MUJAC, MLMAS, MUMAS];

% Constraints
user.opt.ConstraintsAbsTol = [1e-4];
user.opt.PenaltyIncreaseFactor = 1.1;


%% Userparameter definition 
% e.g. system parameters or weights for the cost function
userparam = [500, 0, 1];   % weights

%% Grampc initialization
grampc = CmexFiles.grampc_init_Cmex(userparam);

%% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc,user);

%% Estimate and set PenaltyMin (only needed for constraint optimization)
grampc = CmexFiles.grampc_estim_penmin_Cmex(grampc,1);

%% Simulation time
if nargout>1
    Tsim = 2.9;
end

%% Simulink data
% Converting the grampc struct to a code generation compatible format for
% the use in simulink
if nargout>2
    grampc_sdata = grampc_generate_sdata(grampc);
end
