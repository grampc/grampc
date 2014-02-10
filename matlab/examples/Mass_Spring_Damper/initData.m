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

Tsim = 10.0;

NN = 5;
NX = 2*NN;
NU = 2;

m = 1;
c = 1;
d = 0.2;

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
% mandatory paramter 
u0   = [0.0,0.0]';
udes = [0.0,0.0]';
x0   = zeros(NX,1);
xdes = zeros(NX,1);
x0(1)  = 1.0;
x0(NN) = 1.0;
Thor   = 10.0;
dt     = 0.005;

% optional parameter
umax        = [1,1];
umin        = [-1,-1];
Nhor        = 30;
NpCost      = 2*NX+NU;
pCost       = ones(1,NpCost);
pCost(NX+1) = 0.01;
pCost(NX+2) = 0.01;  
NpSys = 3;
pSys  = [m,c,d];

grampc_setparam_Cmex(grampc,'xk',x0);
grampc_setparam_Cmex(grampc,'u0',u0);
grampc_setparam_Cmex(grampc,'xdes',xdes);
grampc_setparam_Cmex(grampc,'udes',udes);
grampc_setparam_Cmex(grampc,'umax',umax);
grampc_setparam_Cmex(grampc,'umin',umin);
grampc_setparam_Cmex(grampc,'Thor',Thor);
grampc_setparam_Cmex(grampc,'Nhor',Nhor);
grampc_setparam_Cmex(grampc,'dt',dt);
grampc_setparam_Cmex(grampc,'NpCost',NpCost);
grampc_setparam_Cmex(grampc,'pCost',pCost);
grampc_setparam_Cmex(grampc,'NpSys',NpSys);
grampc_setparam_Cmex(grampc,'pSys',pSys);

% Ac = Ac_func(NN,m,c,d);
% Bc = Bc_func(NN,m);
% 
% sys.cont = ss(Ac,Bc,eye(NX,NX),zeros(NX,NU));
% sys.disc = c2d(sys.cont,dt);

% --- END FUNCTION
end

function out = Ac_func(NN,m,c,d)

Ac = zeros(NN);
Ad = zeros(NN);

for i = 1:NN
  Ac(i,i) = -2*c/m;
  Ad(i,i) = -2*d/m;
end
for i = 1:NN-1
  Ac(i,i+1) = c/m;
  Ac(i+1,i) = c/m;
  Ad(i,i+1) = d/m;
  Ad(i+1,i) = d/m;
end

out = [ zeros(NN) , eye(NN)
        Ac        , Ad      ];

% --- EOF ---
end


function out = Bc_func(NN,m)

out = zeros(2*NN,2);

out(NN+1,1) = 1/m;
out(2*NN,2) = -1/m;

% --- EOF ---
end
