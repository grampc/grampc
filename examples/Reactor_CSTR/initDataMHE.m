function [grampcMPC,grampcMHE,Tsim] = initDataMHE()
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

%% Parameter definition MPC
% Initial values and setpoints of the states
userMPC.param.x0    = [2.02e3,1.07e3,100.0,97.1];
userMPC.param.xdes  = [1.37e3,0.95e3,110.0,108.6];

% Initial values, setpoints and limits of the inputs
userMPC.param.u0    = [5.0,-2540.0];
userMPC.param.udes  = [5.0,-1190.0];
userMPC.param.umax  = [35.0,0.0];
userMPC.param.umin  = [3.0,-9000.0];

% Time variables
userMPC.param.Thor  = 1200.0/3600;     % Prediction horizon

userMPC.param.dt    = 1.0/3600;        % Sampling time

%% Option definition MPC
% Basic algorithmic options
userMPC.opt.Nhor         = 40;         % Number of steps for the system integration
userMPC.opt.MaxGradIter  = 3;          % Maximum number of gradient iterations

% Line search
userMPC.opt.LineSearchMax = 1e-6;
userMPC.opt.LineSearchInit = 1e-8;

% Scaling values for the states, inputs and parameters
userMPC.opt.ScaleProblem = 'on';
userMPC.opt.xScale  = [500.0,500.0,50.0,50.0];
userMPC.opt.xOffset = [500.0,500.0,50.0,50.0];
userMPC.opt.uScale  = [16.0,4500.0];
userMPC.opt.uOffset = [19.0,-4500.0];

%% Parameter definition MHE
% Initial values and setpoints of the states
userMHE.param.x0    = [0,0,0,0];
userMHE.param.xdes  = [0,0,0,0];

% Initial values, setpoints and limits of the inputs
userMHE.param.u0    = [0,0];
userMHE.param.udes  = [0,0];
userMHE.param.umax  = [35.0,0.0];
userMHE.param.umin  = [3.0,-9000.0];

% Initial values and limits of the parameters
userMHE.param.p0      = [2.02e3,1.07e3,100.0,97.1]+[100, 100, 5, -7]; % add some disturbances to the initial guess

% Time variables
userMHE.param.Thor  = 10.0/3600;      % Prediction horizon

userMHE.param.dt    = 1.0/3600;       % Sampling time

%% Option definition MHE
% Basic algorithmic options
userMHE.opt.Nhor         = userMHE.param.Thor / userMHE.param.dt + 1;    % Number of steps for the system integration
userMHE.opt.MaxGradIter  = 1;                                        % Maximum number of gradient iterations
userMHE.opt.ShiftControl = 'off';

% Cost integration
userMHE.opt.IntegralCost = 'on';
userMHE.opt.TerminalCost = 'off';
userMHE.opt.IntegratorCost = 'trapezodial';

% System integration
userMHE.opt.Integrator = 'heun';

% Line search
userMHE.opt.LineSearchType = 'explicit2';
userMHE.opt.LineSearchMax = 1e-2;
userMHE.opt.LineSearchMin = 1e-10;
userMHE.opt.LineSearchInit = 1e-8;

% Input and or parameter optimization 
userMHE.opt.OptimControl = 'off';
userMHE.opt.OptimParam = 'on';

% Scaling values for the states, inputs and parameters
userMHE.opt.ScaleProblem = 'on';
userMHE.opt.xScale  = [500.0,500.0,50.0,50.0];
userMHE.opt.xOffset = [500.0,500.0,50.0,50.0];
userMHE.opt.uScale  = [16.0,4500.0];
userMHE.opt.uOffset = [19.0,-4500.0];
userMHE.opt.pScale  = [500.0,500.0,50.0,50.0];
userMHE.opt.pOffset = [500.0,500.0,50.0,50.0];

%% userMPCparameter definition 
% e.g. system parameters or weights for the cost function
pCost  = [0.2,1.0,0.5,0.2,...
  0.2,1.0,0.5,0.2,...
  0.5,5.0e-3]*1e-1;
pCost2  = [1,1]*1e-1;
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

measDummy = repmat(userMPC.param.x0(3:4), 1, userMHE.opt.Nhor);

userMPCparam = [pSys, pCost];
userMHEparam = [pSys, pCost2, measDummy];

%% Grampc initialization
grampcMPC = CmexFiles.grampc_init_Cmex(userMPCparam);
grampcMHE = CmexFiles.grampc_init_Cmex_MHE(userMHEparam);

%% Update grampc struct while ensuring correct data types
grampcMPC = grampc_update_struct_grampc(grampcMPC,userMPC);
grampcMHE = grampc_update_struct_grampc(grampcMHE,userMHE);

%% Simulation time
Tsim = 2.5;

