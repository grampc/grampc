function [vec,grampc,figNr] = startMPC(figNr,compile,varargin)
% This function runs the VEHICLE example. It compiles the c-Code,
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
% This file is part of GRAMPC - (https://github.com/grampc/grampc)
%
% GRAMPC -- A software framework for embedded nonlinear model predictive
% control using a gradient-based augmented Lagrangian approach
%
% Copyright 2014-2025 by Knut Graichen, Andreas Voelz, Thore Wietzke,
% Tobias Englert (<v2.3), Felix Mesmer (<v2.3), Soenke Rhein (<v2.3),
% Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
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
probfct = 'probfct_VEHICLE.c';

% plot predicted trajectories
PLOT_PRED = 1;
% plot solution trajectories
PLOT_TRAJ = 1;
% plot optimization statistics
PLOT_STAT = 1;
% update plots after N steps
PLOT_STEPS = 10;
% pause after each plot
PLOT_PAUSE = 0;

% problem-specific plot
PLOT_MOVE = 1;

% Options for the reference simulation
odeopt =  [];%odeset('RelTol',1e-6,'AbsTol',1e-8);


%% Compilation
% compile toolbox
if compile > 1 || ~exist([grampc_root_path 'matlab/bin'], 'dir')
    grampc_make_toolbox(grampc_root_path, varargin{:});
end
% compile problem
if compile > 0 || ~exist('+CmexFiles', 'dir')
    grampc_make_probfct(grampc_root_path, probfct, varargin{:});
end


%% Initialization
% init GRAMPC and print options and parameters
[grampc,Tsim] = initData;
CmexFiles.grampc_printopt_Cmex(grampc);
CmexFiles.grampc_printparam_Cmex(grampc);

% init solution structure
vec = grampc_init_struct_sol(grampc, Tsim);

% init plots and store figure handles
if figNr > 0 && PLOT_PRED
    phpP = grampc_init_plot_pred(grampc,figNr);
    figNr = figNr+1;
end
if figNr > 0 && PLOT_TRAJ
    phpT = grampc_init_plot_sim(vec,figNr);
    figNr = figNr+1;
end
if figNr > 0 && PLOT_STAT
    phpS = grampc_init_plot_stat(vec,grampc,figNr);
    figNr = figNr+1;
end


%% MPC loop
i = 1;
while 1
    % set current time and current state
    grampc = CmexFiles.grampc_setparam_Cmex(grampc,'t0',vec.t(i));
    grampc = CmexFiles.grampc_setparam_Cmex(grampc,'x0',vec.x(:,i));
    
    % run MPC and save results
    [grampc,vec.CPUtime(i)] = CmexFiles.grampc_run_Cmex(grampc);
    vec = grampc_update_struct_sol(grampc, vec, i);
    
    % print solver status
    printed = CmexFiles.grampc_printstatus_Cmex(grampc.sol.status,'Error');
    if printed
        fprintf('at simulation time %f.\n --------\n', vec.t(i));    
    end
    
    % check for end of simulation
    if i+1 > length(vec.t)
        break;
    end
    
    % simulate system
    [~,xtemp] = ode45(@CmexFiles.grampc_ffct_Cmex,vec.t(i)+[0 double(grampc.param.dt)],vec.x(:,i),odeopt,vec.u(:,i),vec.p(:,i),grampc.param,grampc.userparam);
    vec.x(:,i+1) = xtemp(end,:);
        
    % evaluate time-dependent constraints 
    % to obtain h(x,u,p) instead of max(0,h(x,u,p))
    vec.constr(:,i) = CmexFiles.grampc_ghfct_Cmex(vec.t(i), vec.x(:,i), vec.u(:,i), vec.p(:,i), grampc.param, grampc.userparam);
    
    % update iteration counter
    i = i + 1;
    
    % plot data
    if figNr > 0 && (mod(i,PLOT_STEPS) == 0 || i == length(vec.t))
        if PLOT_PRED
            grampc_update_plot_pred(grampc,phpP);
        end
        if PLOT_TRAJ
            grampc_update_plot_sim(vec,phpT);
        end
        if PLOT_STAT
            grampc_update_plot_stat(vec,grampc,phpS);
        end
        drawnow
        if PLOT_PAUSE
            pause;
        end
    end
end

%% Additional code
if figNr > 0 && PLOT_MOVE == 1
    figure(figNr)
    clf
    figNr = figNr+1;
    steps = 8;
    plot(vec.x(1,1:steps:end), vec.x(2,1:steps:end),'.-','MarkerSize',10);
    hold on;
    plot([-50;50], [6 0;6 0],'k--');
    pos = [18 -1.5 4 4];
    rectangle('Position',pos,'Curvature',[1 1])
    hold off, axis equal,ylim([-10,20]),xlim([-2,vec.x(1,end)+2])
    xlabel('x_1'),ylabel('x_2')
end
end