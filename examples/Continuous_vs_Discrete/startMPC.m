function [vec,grampc,figNr] = startMPC(figNr,use_discrete,compile,varargin)
% This function runs the double integrator example. It compiles the c-Code,
% initializes the grampc struct, runs the simulation and plots the results.
%
% Input arguments are:
% 1) figNr - number of the first plot
% 2) use_discrete - flag whether to use continuous or discrete formulation
% 3) compile - flag whether to compile the whole toolbox or/and the problem function: 
%              1: only the problem function is compiled
%              2: the whole toolbox and the problem function is compiled
%              else or empty input: nothing is compiled
% 4 - end) - flags for the compilation (e.g. 'debug' or 'verbose') see 
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
if nargin < 2 || isempty(use_discrete)
    use_discrete = false;
end
if nargin < 3 || isempty(compile)
    compile = 0;
end


%% Parameters
% path to grampc root
grampc_root_path = '../../';
addpath([grampc_root_path 'matlab/mfiles']);
% name of problem function
probfct_MPC_continuous = 'probfct_MPC_continuous.c';
probfct_MPC_discrete = 'probfct_MPC_discrete.c';

% plot predicted trajectories
PLOT_PRED = 0;
% plot solution trajectories
PLOT_TRAJ = 1;
% plot optimization statistics
PLOT_STAT = 0;
% update plots after N steps
PLOT_STEPS = 1;
% pause after each plot
PLOT_PAUSE = 0;

% Options for the reference simulation
odeopt =  [];%odeset('RelTol',1e-6,'AbsTol',1e-8);


%% Compilation
% compile toolbox
if compile > 1 || ~exist([grampc_root_path 'matlab/bin'], 'dir')
    grampc_make_toolbox(grampc_root_path, varargin{:});
end
% compile problem
if compile > 0 || ~exist('+CmexFiles', 'dir')
    clear mex;
    grampc_make_probfct(grampc_root_path, probfct_MPC_continuous, varargin{:});
    % rename to avoid overwriting
    movefile(['+CmexFiles/grampc_init_Cmex.' mexext], ['+CmexFiles/grampc_init_Cmex_MPC_continuous.' mexext]);
    movefile(['+CmexFiles/grampc_run_Cmex.' mexext], ['+CmexFiles/grampc_run_Cmex_MPC_continuous.' mexext]);
    movefile(['+CmexFiles/grampc_ffct_Cmex.' mexext], ['+CmexFiles/grampc_ffct_Cmex_MPC_continuous.' mexext]);

    grampc_make_probfct(grampc_root_path, probfct_MPC_discrete, varargin{:});
    % rename to avoid overwriting
    movefile(['+CmexFiles/grampc_init_Cmex.' mexext], ['+CmexFiles/grampc_init_Cmex_MPC_discrete.' mexext]);
    movefile(['+CmexFiles/grampc_run_Cmex.' mexext], ['+CmexFiles/grampc_run_Cmex_MPC_discrete.' mexext]);
    movefile(['+CmexFiles/grampc_ffct_Cmex.' mexext], ['+CmexFiles/grampc_ffct_Cmex_MPC_discrete.' mexext]);
end


%% Initialization
% init GRAMPC and print options and parameters
[grampc,Tsim] = initDataMPC(use_discrete);
CmexFiles.grampc_printopt_Cmex(grampc);
CmexFiles.grampc_printparam_Cmex(grampc);

% init solution structure
vec = grampc_init_struct_sol(grampc, Tsim);

% init plots and store figure handles
if PLOT_PRED
    phpP = grampc_init_plot_pred(grampc,figNr);
    figNr = figNr+1;
end
if PLOT_TRAJ
    phpT = grampc_init_plot_sim(vec,figNr);
    figNr = figNr+1;
end
if PLOT_STAT
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
    if use_discrete
        [grampc,vec.CPUtime(i)] = CmexFiles.grampc_run_Cmex_MPC_discrete(grampc);
    else
        [grampc,vec.CPUtime(i)] = CmexFiles.grampc_run_Cmex_MPC_continuous(grampc);
    end    
    vec = grampc_update_struct_sol(grampc, vec, i);

    % print solver status
    printed = CmexFiles.grampc_printstatus_Cmex(grampc.sol.status,'Error');
    if printed
        fprintf('at simulation time %f.\n --------\n', vec.t(i));    
    end
    
    % check for end of simulation
    if i+1 > length(vec.t) || vec.t(i+1) >Tsim
        break;
    end
    
    % simulate system
    if use_discrete
        vec.x(:,i+1) = CmexFiles.grampc_ffct_Cmex_MPC_discrete(vec.t(i), vec.x(:,i), vec.u(:,i), vec.p(:,i), grampc.param, grampc.userparam);
    else
        [~,xtemp] = ode45(@CmexFiles.grampc_ffct_Cmex_MPC_continuous,vec.t(i)+[0 double(grampc.param.dt)],vec.x(:,i),odeopt,vec.u(:,i),vec.p(:,i),grampc.param,grampc.userparam);
        vec.x(:,i+1) = xtemp(end,:);
    end
            
    % update iteration counter
    i = i + 1;
    
    % plot data
    if mod(i,PLOT_STEPS) == 0 || vec.t(i)+grampc.param.dt >= Tsim
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

end

