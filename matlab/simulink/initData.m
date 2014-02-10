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

Tsim = 4.0;

grampc = grampc_init_Cmex();

%% Options (optional)
% grampc_setopt_Cmex(grampc,'MaxIter',2);
% grampc_setopt_Cmex(grampc,'ShiftControl','on');
% grampc_setopt_Cmex(grampc,'ScaleProblem','on');
% grampc_setopt_Cmex(grampc,'CostIntegrator','simpson');
% grampc_setopt_Cmex(grampc,'Integrator','euler');
% grampc_setopt_Cmex(grampc,'IntegratorRelTol',1e-3);
% grampc_setopt_Cmex(grampc,'IntegratorAbsTol',1e-3);
% grampc_setopt_Cmex(grampc,'LineSearchType','explicit2');
% grampc_setopt_Cmex(grampc,'LineSearchMax',2.0);
% grampc_setopt_Cmex(grampc,'LineSearchMin',1e-3);
% grampc_setopt_Cmex(grampc,'LineSearchInit',0.1);
% grampc_setopt_Cmex(grampc,'LineSearchIntervalFactor',0.5);
% grampc_setopt_Cmex(grampc,'LineSearchAdaptFactor',2.0);
% grampc_setopt_Cmex(grampc,'LineSearchIntervalTol',1e-3);
% grampc_setopt_Cmex(grampc,'JacobianX','sysjacx');
% grampc_setopt_Cmex(grampc,'JacobianU','sysjacu');
% grampc_setopt_Cmex(grampc,'IntegralCost','on');
% grampc_setopt_Cmex(grampc,'FinalCost','off');

%% Parameter
% mandatory parameter
x0      = ;
u0      = ;
xdes    = ;
udes    = ;
Thor    = ;
dt      = ;
% optional parameter
% tk      = ;
% umax    = ;
% umin    = ;
% Nhor    = ;
% NpCost  = ;
% pCost   = ;
% NpSys   = ;
% pSys    = ;
% xScale  = ;
% xOffset = ;
% uScale  = ;
% uOffset = ;

grampc_setparam_Cmex(grampc,'xk',x0);
grampc_setparam_Cmex(grampc,'u0',u0);
grampc_setparam_Cmex(grampc,'xdes',xdes);
grampc_setparam_Cmex(grampc,'udes',udes);
% grampc_setparam_Cmex(grampc,'dt',dt);
% grampc_setparam_Cmex(grampc,'tk',tk);
% grampc_setparam_Cmex(grampc,'umax',umax);
% grampc_setparam_Cmex(grampc,'umin',umin);
% grampc_setparam_Cmex(grampc,'Thor',Thor);
% grampc_setparam_Cmex(grampc,'Nhor',Nhor);
% grampc_setparam_Cmex(grampc,'NpCost',NpCost);
% grampc_setparam_Cmex(grampc,'pCost',pCost);
% grampc_setparam_Cmex(grampc,'NpSys',NpSys);
% grampc_setparam_Cmex(grampc,'pSys',pSys);
% grampc_setparam_Cmex(grampc,'xScale',xScale);
% grampc_setparam_Cmex(grampc,'uScale',uScale);
% grampc_setparam_Cmex(grampc,'xOffset',xScale);
% grampc_setparam_Cmex(grampc,'uScale',uScale);
							
% --- END FUNCTION
end
