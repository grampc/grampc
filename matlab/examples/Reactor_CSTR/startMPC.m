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
  PROBFCT = 'probfct_REACTOR_CSTR.c';
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
  % grampc_setopt_Cmex(grampc,'MaxIter',2);
  % grampc_setopt_Cmex(grampc,'ShiftControl','on');
  grampc_setopt_Cmex(grampc,'ScaleProblem','on');
  % grampc_setopt_Cmex(grampc,'CostIntegrator','simpson');
  % grampc_setopt_Cmex(grampc,'Integrator','euler');
  % grampc_setopt_Cmex(grampc,'IntegratorRelTol',1e-3);
  % grampc_setopt_Cmex(grampc,'IntegratorAbsTol',1e-3);
  grampc_setopt_Cmex(grampc,'LineSearchType','explicit1');
  % grampc_setopt_Cmex(grampc,'LineSearchMax',2.0);
  % grampc_setopt_Cmex(grampc,'LineSearchMin',1e-9);
  grampc_setopt_Cmex(grampc,'LineSearchInit',1e-6);
  % grampc_setopt_Cmex(grampc,'LineSearchIntervalFactor',0.5);
  % grampc_setopt_Cmex(grampc,'LineSearchAdaptFactor',2.0);
  % grampc_setopt_Cmex(grampc,'LineSearchIntervalTol',0.2);
  % grampc_setopt_Cmex(grampc,'JacobianX','sysjacx');
  % grampc_setopt_Cmex(grampc,'JacobianU','sysjacu');
  % grampc_setopt_Cmex(grampc,'IntegralCost','on');
  % grampc_setopt_Cmex(grampc,'FinalCost','on');



  %% User parameter
  % mandatory paramter
  user.param.x0     = [2.02e3,1.07e3,100.0,97.1];
  user.param.u0     = [5.0,-2540.0];
  user.param.xdes   = [1.37e3,0.95e3,110.0,108.6];
  user.param.udes   = [5.0,-1190.0];
  user.param.Thor   = 1200.0/3600;
  user.param.dt     = 1.0/3600;
  % optional parameter
  user.param.xScale  = [500.0,500.0,50.0,50.0];
  user.param.xOffset = [500.0,500.0,50.0,50.0];
  user.param.uScale  = [16.0,4500.0];
  user.param.uOffset = [19.0,-4500.0];
  user.param.umax   = [35.0,0.0];
  user.param.umin   = [3.0,-9000.0];
  user.param.Nhor   = 40;
  user.param.NpCost = 10;
  user.param.NpSys  = 14;
  user.param.pCost  = [0.2,1.0,0.5,0.2,...
		       0.2,1.0,0.5,0.2,...
		       0.5,5.0e-3];
  user.param.pSys   = [1.287e12,...
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

  grampc_setparam_Cmex(grampc,'umax',user.param.umax);
  grampc_setparam_Cmex(grampc,'umin',user.param.umin);
  grampc_setparam_Cmex(grampc,'xScale',user.param.xScale);
  grampc_setparam_Cmex(grampc,'uScale',user.param.uScale);
  grampc_setparam_Cmex(grampc,'xOffset',user.param.xOffset);
  grampc_setparam_Cmex(grampc,'uOffset',user.param.uOffset);
  grampc_setparam_Cmex(grampc,'Thor',user.param.Thor);
  grampc_setparam_Cmex(grampc,'dt',user.param.dt);
  grampc_setparam_Cmex(grampc,'NpCost',user.param.NpCost);
  grampc_setparam_Cmex(grampc,'pCost',user.param.pCost);
  grampc_setparam_Cmex(grampc,'NpSys',user.param.NpSys);
  grampc_setparam_Cmex(grampc,'pSys',user.param.pSys);
  grampc_setparam_Cmex(grampc,'xk',user.param.x0);
  grampc_setparam_Cmex(grampc,'u0',user.param.u0);
  grampc_setparam_Cmex(grampc,'xdes',user.param.xdes);
  grampc_setparam_Cmex(grampc,'udes',user.param.udes);
  grampc_setparam_Cmex(grampc,'Nhor',user.param.Nhor);

  %% MPC loop

  Tsim = 0.6;

  i=1;
  vec.t = 0:grampc.param.dt:Tsim;
  vec.u = zeros(grampc.param.Nu,length(vec.t));
  vec.x = zeros(grampc.param.Nx,length(vec.t));
  vec.x(:,1) = grampc.param.xk;
  vec.J = zeros(1,length(vec.t));
  vec.CPUtime = zeros(1,length(vec.t)-1);

  while (vec.t(i)<Tsim)
	
    tic;
    [xnext,unext,J]           = grampc_run_Cmex(grampc);
    grampc_setparam_Cmex(grampc,'xk',xnext);
    vec.CPUtime(i) = toc;
    
    vec.x(:,i+1) = xnext;
    vec.u(:,i+1) = unext;
    vec.J(i+1)   = J;
    
    i = i + 1;
    
  end

% ******* END OF FUNCTION *******
end
