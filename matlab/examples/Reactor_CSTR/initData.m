function [grampc,Tsim] = initData()
%
% This file is part of GRAMPC.
%
% GRAMPC - a gradient-based MPC software for real-time applications
%
% Copyright (C) 2014 by Bartosz Kaepernick, Knut Graichen, Tilman Utz
% Developed at the Institute of Measurement, Control, and
% Microtechnology, University of Ulm. All rights reserved.
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
% File: initDate.m
% Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
% Date: February 2014
% Version: v1.0
%
% MATLAB/Simulink/GUI initialization file of GRAMPC.
%

Tsim = 0.6;

grampc = grampc_init_Cmex();

%% Options (optional)
% grampc_setopt_Cmex(grampc,'MaxIter',2);
% grampc_setopt_Cmex(grampc,'ShiftControl','on');
grampc_setopt_Cmex(grampc,'ScaleProblem','on');
% grampc_setopt_Cmex(grampc,'CostIntegrator','simpson');
% grampc_setopt_Cmex(grampc,'Integrator','euler');
% grampc_setopt_Cmex(grampc,'IntegratorRelTol',1e-3);
% grampc_setopt_Cmex(grampc,'IntegratorAbsTol',1e-3);
% grampc_setopt_Cmex(grampc,'LineSearchType','explicit2');
grampc_setopt_Cmex(grampc,'LineSearchMax',2.0);
grampc_setopt_Cmex(grampc,'LineSearchMin',1e-9);
grampc_setopt_Cmex(grampc,'LineSearchInit',1e-6);
% grampc_setopt_Cmex(grampc,'LineSearchIntervalFactor',0.5);
grampc_setopt_Cmex(grampc,'LineSearchAdaptFactor',2.0);
% grampc_setopt_Cmex(grampc,'LineSearchIntervalTol',1e-3);
% grampc_setopt_Cmex(grampc,'JacobianX','sysjacx');
% grampc_setopt_Cmex(grampc,'JacobianU','sysjacu');
% grampc_setopt_Cmex(grampc,'IntegralCost','on');
% grampc_setopt_Cmex(grampc,'FinalCost','off');

%% Parameter
% mandatory paramter
x0     = [2.02e3,1.07e3,100.0,97.1];
u0     = [5.0,-2540.0];
xdes   = [1.37e3,0.95e3,110.0,108.6];
udes   = [5.0,-1190.0];
dt     = 1.0/3600;
Thor   = 1200.0/3600;
% optional parameter
xScale  = [500.0,500.0,50.0,50.0];
xOffset = [500.0,500.0,50.0,50.0];
uScale  = [16.0,4500.0];
uOffset = [19.0,-4500.0];

umax   = [35.0,0.0];
umin   = [3.0,-9000.0];
Nhor   = 40;
NpCost = 10;
NpSys  = 14;
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

grampc_setparam_Cmex(grampc,'umax',umax);
grampc_setparam_Cmex(grampc,'umin',umin);
grampc_setparam_Cmex(grampc,'xScale',xScale);
grampc_setparam_Cmex(grampc,'uScale',uScale);
grampc_setparam_Cmex(grampc,'xOffset',xOffset);
grampc_setparam_Cmex(grampc,'uOffset',uOffset);
grampc_setparam_Cmex(grampc,'xk',x0);
grampc_setparam_Cmex(grampc,'u0',u0);
grampc_setparam_Cmex(grampc,'xdes',xdes);
grampc_setparam_Cmex(grampc,'udes',udes);
grampc_setparam_Cmex(grampc,'Thor',Thor);
grampc_setparam_Cmex(grampc,'Nhor',Nhor);
grampc_setparam_Cmex(grampc,'dt',dt);
grampc_setparam_Cmex(grampc,'NpCost',NpCost);
grampc_setparam_Cmex(grampc,'pCost',pCost);
grampc_setparam_Cmex(grampc,'NpSys',NpSys);
grampc_setparam_Cmex(grampc,'pSys',pSys);

% --- END FUNCTION
end
