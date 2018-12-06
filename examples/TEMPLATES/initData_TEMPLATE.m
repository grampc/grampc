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
Nx = 1;
Nc = 1;

% Initial values and setpoints of the states
user.param.x0    = [ ];
user.param.xdes  = [ ];

% Initial values, setpoints and limits of the inputs
user.param.u0    = [ ];
user.param.udes  = [ ];
user.param.umax  = [ ];
user.param.umin  = [ ];

% Initial values and limits of the parameters
user.param.p0    = [ ];
user.param.pmax  = [ ];
user.param.pmin  = [ ];

% Time variables
user.param.Thor  = -1;          % Prediction horizon
user.param.Tmax  = 1e8;
user.param.Tmin  = 1e-8;

user.param.dt    = -1;          % Sampling time
user.param.t0    = 0.0;         % time at the current sampling step

%% Option definition
% Basic algorithmic options
user.opt.Nhor        = 30;          % Number of steps for the system integration
user.opt.MaxGradIter = 2;           % Maximum number of gradient iterations
user.opt.MaxMultIter = 1;           % Maximum number of augmented Lagrangian iterations
user.opt.ShiftControl = 'on';

% Cost integration
user.opt.IntegralCost = 'on';
user.opt.TerminalCost = 'on';
user.opt.IntegratorCost = 'trapezodial';

% System integration
user.opt.Integrator = 'heun';
user.opt.IntegratorRelTol = 1e-6;
user.opt.IntegratorAbsTol = 1e-8;
user.opt.IntegratorMinStepSize = 1e-16;
user.opt.IntegratorMaxSteps = 1e8;
user.opt.FlagsRodas = [0,0,0,Nx,Nx,0,Nx,Nx];

% Line search
user.opt.LineSearchType = 'explicit2';
user.opt.LineSearchExpAutoFallback = 'on';
user.opt.LineSearchMax = 0.75;
user.opt.LineSearchMin = 1e-10;
user.opt.LineSearchInit = 1e-4;
user.opt.LineSearchAdaptAbsTol = 1e-6;
user.opt.LineSearchAdaptFactor = 3.0/ 2.0;
user.opt.LineSearchIntervalTol = 1e-1;
user.opt.LineSearchIntervalFactor = 0.85;

% Input and or parameter optimization 
user.opt.OptimControl = 'on';
user.opt.OptimParam = 'off';
user.opt.OptimParamLineSearchFactor = 1.0;
user.opt.OptimTime = 'off';
user.opt.OptimTimeLineSearchFactor = 1.0;

% Scaling values for the states, inputs and parameters
user.opt.ScaleProblem = 'off';
user.opt.xScale  = [];
user.opt.xOffset = [];
user.opt.uScale  = [];
user.opt.uOffset = [];
user.opt.pScale  = [];
user.opt.pOffset = [];
user.opt.TScale  = 1;
user.opt.TOffset = 0;
user.opt.JScale  = 1;
user.opt.cScale  = [];

% Type of considered constraints
user.opt.EqualityConstraints = 'on';
user.opt.InequalityConstraints = 'on';
user.opt.TerminalEqualityConstraints = 'on';
user.opt.TerminalInequalityConstraints = 'on';
user.opt.ConstraintsHandling = 'auglag';    
user.opt.ConstraintsAbsTol = 1e-4 * ones(1,Nc);

% Multipliers & penalties
user.opt.MultiplierMax = 1e6;
user.opt.MultiplierDampingFactor = 0;
user.opt.PenaltyMax = 1e6;
user.opt.PenaltyMin = 1;
user.opt.PenaltyIncreaseFactor = 1.05;
user.opt.PenaltyDecreaseFactor = 0.95;
user.opt.PenaltyIncreaseThreshold = 1.0;
user.opt.AugLagUpdateGradientRelTol = 1e-2;

% Convergences test
user.opt.ConvergenceCheck = 'off';
user.opt.ConvergenceGradientRelTol = 1e-6;

%% Userparameter definition 
% e.g. system parameters or weights for the cost function
pSys = [ ];    % system parameters
pCost = [ ];   % weights
    
userparam = [pSys, pCost];

%% Grampc initialization
grampc = CmexFiles.grampc_init_Cmex(userparam);

%% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc,user);

%% Estimate and set PenaltyMin (only needed for constraint optimization)
grampc = CmexFiles.grampc_estim_penmin_Cmex(grampc,1);

%% Simulation time
if nargout>1
    Tsim = 20;
end

%% Simulink data
% Converting the grampc struct to a code generation compatible format for
% the use in simulink
if nargout>2
    grampc_sdata = grampc_generate_sdata(grampc);
end
