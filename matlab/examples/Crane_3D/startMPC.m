function [vec,grampc] = startMPC()
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
  % File: startMPC.m
  % Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
  % Date: February 2014
  % Version: v1.0
  %
  % MATLAB/Simulink file for running GRAMPC example.
  %


  %% Compilation
  PROBFCT = 'probfct_CRANE_3D.c';
  COMPILE = 0;
  DEBUG   = 0;
  VERBOSE = 0;
  CLEAN   = 0;

  if COMPILE || CLEAN
    if exist('make','file') ~= 2
      curPath = cd;
      cd('..');
      addpath(cd);
      cd(curPath);
    end
  end

  if CLEAN
    make clean;
    vec    = [];
    grampc = [];
    return;
  end

  if COMPILE
    if VERBOSE && DEBUG
      make(PROBFCT,'verbose','debug');
    elseif VERBOSE
      make(PROBFCT,'verbose');
    elseif DEBUG
      make(PROBFCT,'debug');
    else
      make(PROBFCT);
    end
  end


  %% Initialization
  grampc = grampc_init_Cmex();


  %% user options (all options are optional)
  % grampc_setopt_Cmex(grampc,'MaxIter',5);
  % grampc_setopt_Cmex(grampc,'ShiftControl','on');
  % grampc_setopt_Cmex(grampc,'ScaleProblem','on');
  % grampc_setopt_Cmex(grampc,'CostIntegrator','simpson');
  % grampc_setopt_Cmex(grampc,'Integrator','euler');
  % grampc_setopt_Cmex(grampc,'IntegratorRelTol',1e-3);
  % grampc_setopt_Cmex(grampc,'IntegratorAbsTol',1e-3);
  grampc_setopt_Cmex(grampc,'LineSearchType','adaptive');
  % grampc_setopt_Cmex(grampc,'LineSearchMax',2.0);
  % grampc_setopt_Cmex(grampc,'LineSearchMin',1e-3);
  % grampc_setopt_Cmex(grampc,'LineSearchInit',0.1);
  % grampc_setopt_Cmex(grampc,'LineSearchIntervalFactor',0.5);
  % grampc_setopt_Cmex(grampc,'LineSearchAdaptFactor',2.0);
  % grampc_setopt_Cmex(grampc,'LineSearchIntervalTol',0.2);
  grampc_setopt_Cmex(grampc,'JacobianX','sysjacxadj');
  grampc_setopt_Cmex(grampc,'JacobianU','sysjacuadj');
  % grampc_setopt_Cmex(grampc,'IntegralCost','on');
  % grampc_setopt_Cmex(grampc,'FinalCost','off');


  %% User parameter
  % mandatory paramter 
  user.param.x0     = [0.7,0.7,-1.0472,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
  user.param.u0     = [0.0,0.0,0.0];
  user.param.xdes   = [0.2,0.25,1.0472,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
  user.param.udes   = [0.0,0.0,0.0];
  user.param.Thor   = 1.5;
  user.param.dt     = 0.002;
  % optional parameter
  user.param.umax   = [2.0,2.0,2.0];
  user.param.umin   = [-2.0,-2.0,-2.0];
  user.param.Nhor   = 40;
  user.param.NpCost = 23;
  user.param.pCost  = [1.0,1.0,1.0,1.0,1.0,0.1,0.1,0.1,0.1,1.0,...
		       0.01,0.01,0.01,...
		       1.0,1.0,1.0,1.0,1.0,0.1,0.1,0.1,0.1,1.0];

  grampc_setparam_Cmex(grampc,'xk',user.param.x0);
  grampc_setparam_Cmex(grampc,'u0',user.param.u0);
  grampc_setparam_Cmex(grampc,'xdes',user.param.xdes);
  grampc_setparam_Cmex(grampc,'udes',user.param.udes);
  grampc_setparam_Cmex(grampc,'umax',user.param.umax);
  grampc_setparam_Cmex(grampc,'umin',user.param.umin);
  grampc_setparam_Cmex(grampc,'Thor',user.param.Thor);
  grampc_setparam_Cmex(grampc,'Nhor',user.param.Nhor);
  grampc_setparam_Cmex(grampc,'dt',user.param.dt);
  grampc_setparam_Cmex(grampc,'NpCost',user.param.NpCost);
  grampc_setparam_Cmex(grampc,'pCost',user.param.pCost);


  %% MPC loop

  Tsim = 4;

  i=1;
  vec.t = 0:grampc.param.dt:Tsim;  
  vec.u = zeros(grampc.param.Nu,length(vec.t));
  vec.x = zeros(grampc.param.Nx,length(vec.t));
  vec.x(:,1) = grampc.param.xk;
  vec.J = zeros(1,length(vec.t));
  vec.CPUtime = zeros(1,length(vec.t)-1);

  while (vec.t(i)<Tsim)

    tic;
    [xnext,unext,J] = grampc_run_Cmex(grampc);
    grampc_setparam_Cmex(grampc,'xk',xnext);
    vec.CPUtime(i) = toc;

    vec.x(:,i+1) = xnext;
    vec.u(:,i+1) = unext;
    vec.J(i+1)   = J;

    i = i + 1;

  end

% ******* END OF FUNCTION *******
end
