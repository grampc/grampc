function [vec,grampcMPC,figNr] = startMHE(figNr,compile,varargin)
% This function runs the Reactor CSTR MHE example. It compiles the c-Code,
% initializes the grampc struct, runs the simulation and plots the results.
%
% Input arguments are:
% 1) figNr - number of the first plot
% 2) compile - flag whether to compile the whole toolbox or/and the problem function: 
%              1: only the problem function is compiled
%              2: the whole toolbox and the problem function is compiled
%              else or empty input: nothing is compiled
% 3 - end) - flags for the compilation (e.g. 'debug' or 'verbose') see 
%            make.m in the matlab folder for more details
%
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

%% Check input arguments
if nargin < 1 || isempty(figNr)
    figNr = 1;
end
if nargin < 2 || isempty(compile)
    compile = 0;
end


%% Parameters
% path to grampc root
grampc_root_path = '../../';
addpath([grampc_root_path 'matlab/mfiles']);
% name of problem function
probfctMPC = 'probfct_REACTOR_CSTR.c';
probfctMHE = 'probfct_REACTOR_CSTR_MHE.c';

% plot predicted trajectories
PLOT_PRED = 1;
% plot solution trajectories
PLOT_TRAJ = 1;
% plot optimization statistics
PLOT_STAT = 1;
% update plots after N steps
PLOT_STEPS = 100;
% pause after each plot
PLOT_PAUSE = 0;

% problem-specific plot
PLOT_MHE = 1;

% Options for the reference simulation
odeopt =  [];%odeset('RelTol',1e-6,'AbsTol',1e-8);


%% Compilation
% compile toolbox
if compile > 1 || ~exist([grampc_root_path 'matlab/bin'], 'dir')
    grampc_make_toolbox(grampc_root_path, varargin{:});
end
% compile problem
if compile > 0 || ~exist('+CmexFiles', 'dir')
    % compile MHE
    grampc_make_probfct(grampc_root_path, probfctMHE, varargin{:});
    copyfile(['+CmexFiles/grampc_init_Cmex.' mexext], ['+CmexFiles/grampc_init_Cmex_MHE.' mexext]);
    copyfile(['+CmexFiles/grampc_run_Cmex.' mexext], ['+CmexFiles/grampc_run_Cmex_MHE.' mexext]);
    % compile MPC
    grampc_make_probfct(grampc_root_path, probfctMPC, varargin{:});
end


%% Initialization
% init GRAMPC and print options and parameters
[grampcMPC, grampcMHE, Tsim] = initDataMHE;
CmexFiles.grampc_printopt_Cmex(grampcMPC);
CmexFiles.grampc_printparam_Cmex(grampcMPC);

% init solution structure
vec = grampc_init_struct_sol(grampcMPC, Tsim);

% init additional fields for MHE
Nsim = length(vec.t);
vec.pMHE        = nan * zeros(grampcMHE.param.Np, Nsim);
vec.xEst        = nan * zeros(grampcMPC.param.Nx, Nsim);
vec.xMeas       = nan * zeros(2, Nsim);
vec.JMHE        = nan * zeros(2, Nsim);
vec.CPUtimeMHE  = nan * zeros(1, Nsim);

% init plots and store figure handles
if PLOT_PRED
    phpP = grampc_init_plot_pred(grampcMPC,figNr);
    figNr = figNr+1;
end
if PLOT_TRAJ
    phpT = grampc_init_plot_sim(vec,figNr);
    figNr = figNr+1;
end
if PLOT_STAT
    phpS = grampc_init_plot_stat(vec,grampcMPC,figNr);
    figNr = figNr+1;
end
if PLOT_MHE
    figNrMHE = figNr+1;
end

%% MPC loop
% init array of last MHE-Nhor measurements of the two temperatures
xMeas_array = repmat(grampcMPC.param.x0(3:4), 1, grampcMHE.opt.Nhor);
i = 1;
while 1
    % get current state estimate from MHE
%     vec.xEst(:,i) = vec.x(:,i);
    vec.xEst(:,i) = (grampcMHE.rws.p + grampcMHE.rws.x(:,end) - grampcMHE.rws.x(:,1)) .* grampcMHE.opt.xScale + grampcMHE.opt.xOffset;
    
    % set current time and current state estimate
    grampcMPC = CmexFiles.grampc_setparam_Cmex(grampcMPC,'t0',vec.t(i));
    grampcMPC = CmexFiles.grampc_setparam_Cmex(grampcMPC,'x0',vec.xEst(:,i));
        
    % switch setpoint at 1.5 hours
    if (((vec.t(i) - grampcMPC.param.dt) < 1.5) && ((vec.t(i) + grampcMPC.param.dt) > 1.5))
        userMPC.param.xdes  = [2.02e3,1.07e3,100.0,97.1];
        userMPC.param.udes  = [5.0,-2540.0];
        grampcMPC = CmexFiles.grampc_setparam_Cmex(grampcMPC,'xdes',userMPC.param.xdes);
        grampcMPC = CmexFiles.grampc_setparam_Cmex(grampcMPC,'udes',userMPC.param.udes);
    end
        
    % run MPC and save results
    [grampcMPC,vec.CPUtime(i)] = CmexFiles.grampc_run_Cmex(grampcMPC);
    vec = grampc_update_struct_sol(grampcMPC, vec, i);
        
    % print solver status
    printed = CmexFiles.grampc_printstatus_Cmex(grampcMPC.sol.status,'Error');
    if printed
        fprintf('at simulation time %f.\n --------\n', vec.t(i));    
    end
    
    % check for end of simulation
    if i+1 > length(vec.t)
        break;
    end
    
    % simulate system
    [~,xtemp] = ode45(@CmexFiles.grampc_ffct_Cmex,vec.t(i)+[0 (grampcMPC.param.dt)],vec.x(:,i),odeopt,(grampcMPC.sol.unext),(grampcMPC.sol.pnext),(grampcMPC.userparam));
%     xtemp(end,1:2) = xtemp(end,1:2) + randn(1,2)*10; % system noise
    vec.x(:,i+1) = xtemp(end,:);
            
    % set values of last MHE-Nhor measurements of the two temperatures
    xMeas_temp = xtemp(end,3:4) + randn(1,2)*4; % measurement noise
    xMeas_array = [xMeas_array(3:end), xMeas_temp];
    grampcMHE.userparam(end-2*grampcMHE.opt.Nhor+1:end) = xMeas_array;
    % set values of last MHE-Nhor controls or use initial values 
    if i > 10
        grampcMHE.rws.u = (vec.u(:,i-10:i) - repmat(grampcMHE.opt.uOffset,1,11))./repmat(grampcMHE.opt.uScale,1,11);
    else
        grampcMHE.rws.u(:,10-i+1:10) = (vec.u(:,1:i) - repmat(grampcMHE.opt.uOffset,1,i))./repmat(grampcMHE.opt.uScale,1,i);
    end    
    % run MHE and save results
    [grampcMHE,vec.CPUtimeMHE(i)] = CmexFiles.grampc_run_Cmex_MHE(grampcMHE);
    vec.xMeas(:,i)  = xMeas_temp';
    vec.pMHE(:,i)   = grampcMHE.sol.pnext;
    vec.JMHE(:,i)   = grampcMHE.sol.J;
    % predict next state to initialize p in the next step
%     grampcMHE.rws.p = grampcMHE.rws.p + grampcMHE.rws.x(:,2) + grampcMHE.opt.xOffset ./ grampcMHE.opt.xScale;
    grampcMHE.rws.p = grampcMHE.rws.p + grampcMHE.rws.x(:,2) - grampcMHE.rws.x(:,1);
    
    % update iteration counter
    i = i + 1;
    
    % plot data
    if mod(i,PLOT_STEPS) == 0 || i == length(vec.t)
        if PLOT_PRED
            grampc_update_plot_pred(grampcMPC,phpP);
        end
        if PLOT_TRAJ
            grampc_update_plot_sim(vec,phpT);
        end
        if PLOT_STAT
            grampc_update_plot_stat(vec,grampcMPC,phpS);
        end
        % problem-specific plot for MHE
        if PLOT_MHE
            figure(figNrMHE)
            clf
            subplot(3,2,1)
            plot(0,0)
            hold on
            plot(vec.t(1:PLOT_STEPS:i), vec.xEst(1,1:PLOT_STEPS:i)', 'r')
            plot(vec.t(1:PLOT_STEPS:i), vec.x(1,1:PLOT_STEPS:i)', 'b')
            xlim([0, Tsim]);
            ylabel('Monomer concentration c_A [kmol m^{-3}]');
            subplot(3,2,2)
            plot(0,0)
            hold on
            plot(vec.t(1:PLOT_STEPS:i), vec.xEst(2,1:PLOT_STEPS:i)', 'r')
            plot(vec.t(1:PLOT_STEPS:i), vec.x(2,1:PLOT_STEPS:i)', 'b')
            xlim([0, Tsim]);
            ylabel('Product concentration c_B [kmol m^{-3}]');
            subplot(3,2,3)
            plot(vec.t(1:PLOT_STEPS:i), vec.xMeas(1,1:PLOT_STEPS:i)','x')
            hold on
            plot(vec.t(1:PLOT_STEPS:i), vec.xEst(3,1:PLOT_STEPS:i)', 'r')
            plot(vec.t(1:PLOT_STEPS:i), vec.x(3,1:PLOT_STEPS:i)', 'b')
            xlim([0, Tsim]);
            ylabel('Reactor temperature T [°C]');
            subplot(3,2,4)
            plot(vec.t(1:PLOT_STEPS:i), vec.xMeas(2,1:PLOT_STEPS:i)','x')
            hold on
            plot(vec.t(1:PLOT_STEPS:i), vec.xEst(4,1:PLOT_STEPS:i)', 'r')
            plot(vec.t(1:PLOT_STEPS:i), vec.x(4,1:PLOT_STEPS:i)', 'b')
            legend('Measured state','Estimated state','Simulated state','Location','NorthEast');
            xlim([0, Tsim]);
            ylabel('Cooling temperature T_C [°C]');
            subplot(3,2,5)
            plot(vec.t(1:PLOT_STEPS:i), vec.u(1,1:PLOT_STEPS:i)', 'b')
            ylabel('Normalized flow rate u_1 [h^{-1}]');
            xlabel('Time [h]');
            xlim([0, Tsim]);
            subplot(3,2,6)
            plot(vec.t(1:PLOT_STEPS:i), vec.u(2,1:PLOT_STEPS:i)', 'b')
            ylabel('Cooling power u_2 [kJh^{-1}]');
            xlabel('Time [h]');
            xlim([0, Tsim]);
%             subplot(3,1,3)
%             plot(vec.t(1:i), vec.JMHE(:,1:i))
        end
        drawnow
        if PLOT_PAUSE
            pause;
        end
    end
end

%% Additional code
disp(['Average computation time for the MPC: ', num2str(mean(vec.CPUtime(1:end-1))), ' and for the MHE: ', num2str(mean(vec.CPUtimeMHE(1:end-1)))]);

end
