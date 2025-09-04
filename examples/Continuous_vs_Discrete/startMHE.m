function [vec, grampcMPC, grampcMHE, figNr] = startMHE(figNr, use_discrete, compile, varargin)
% This function runs a simple MHE example for a double integrator to
% compare continuous-time and discrete-time implementations. It compiles
% the C-code, initializes the grampc struct, runs the simulation and plots
% the results.
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
probfct_MHE_continuous = 'probfct_MHE_continuous.c';
probfct_MPC_discrete   = 'probfct_MPC_discrete.c';
probfct_MHE_discrete   = 'probfct_MHE_discrete.c';

% plot predicted trajectories
PLOT_PRED = 0;
% plot solution trajectories
PLOT_TRAJ = 0;
% plot optimization statistics
PLOT_STAT = 0;
% plot estimation over MHE horizon
PLOT_EST = 0;
% plot estimation over simulation time
PLOT_MHE = 1;
% update plots after N steps
PLOT_STEPS = 1;
% pause after each plot
PLOT_PAUSE = 0;

% Options for the reference simulation
odeopt =  []; %odeset('RelTol',1e-6,'AbsTol',1e-8);


%% Compilation
% compile toolbox
if compile > 1 || ~exist([grampc_root_path 'matlab/bin'], 'dir')
    grampc_make_toolbox(grampc_root_path, varargin{:});
end
% compile problem
if compile > 0 || ~exist('+CmexFiles', 'dir')
    clear mex
    grampc_make_probfct(grampc_root_path, probfct_MPC_continuous, varargin{:});
    % rename to avoid overwriting
    movefile(['+CmexFiles/grampc_init_Cmex.' mexext], ['+CmexFiles/grampc_init_Cmex_MPC_continuous.' mexext]);
    movefile(['+CmexFiles/grampc_run_Cmex.' mexext], ['+CmexFiles/grampc_run_Cmex_MPC_continuous.' mexext]);
    movefile(['+CmexFiles/grampc_ffct_Cmex.' mexext], ['+CmexFiles/grampc_ffct_Cmex_MPC_continuous.' mexext]);

    grampc_make_probfct(grampc_root_path, probfct_MHE_continuous, varargin{:});
    % rename to avoid overwriting
    movefile(['+CmexFiles/grampc_init_Cmex.' mexext], ['+CmexFiles/grampc_init_Cmex_MHE_continuous.' mexext]);
    movefile(['+CmexFiles/grampc_run_Cmex.' mexext], ['+CmexFiles/grampc_run_Cmex_MHE_continuous.' mexext]);

    grampc_make_probfct(grampc_root_path, probfct_MPC_discrete, varargin{:});
    % rename to avoid overwriting
    movefile(['+CmexFiles/grampc_init_Cmex.' mexext], ['+CmexFiles/grampc_init_Cmex_MPC_discrete.' mexext]);
    movefile(['+CmexFiles/grampc_run_Cmex.' mexext], ['+CmexFiles/grampc_run_Cmex_MPC_discrete.' mexext]);
    movefile(['+CmexFiles/grampc_ffct_Cmex.' mexext], ['+CmexFiles/grampc_ffct_Cmex_MPC_discrete.' mexext]);

    grampc_make_probfct(grampc_root_path, probfct_MHE_discrete, varargin{:});
    % rename to avoid overwriting
    movefile(['+CmexFiles/grampc_init_Cmex.' mexext], ['+CmexFiles/grampc_init_Cmex_MHE_discrete.' mexext]);
    movefile(['+CmexFiles/grampc_run_Cmex.' mexext], ['+CmexFiles/grampc_run_Cmex_MHE_discrete.' mexext]);
end


%% Initialization
% init MPC and print options and parameters
[grampcMPC, Tsim] = initDataMPC(use_discrete);
% CmexFiles.grampc_printopt_Cmex(grampcMPC);
% CmexFiles.grampc_printparam_Cmex(grampcMPC);

% init solution structure
vec = grampc_init_struct_sol(grampcMPC, Tsim);

% init MHE and print options and parameters
[grampcMHE, yfct] = initDataMHE(use_discrete, grampcMPC.param.dt);
% CmexFiles.grampc_printopt_Cmex(grampcMHE);
% CmexFiles.grampc_printparam_Cmex(grampcMHE);

% init estimation structure
est.t       = grampcMHE.rws.t;
est.u       = zeros(grampcMHE.param.Nu, grampcMHE.opt.Nhor);
est.x       = zeros(grampcMHE.param.Nx, grampcMHE.opt.Nhor);
est.xest    = est.x;
est.y       = yfct(est.x);
est.yest    = est.y;
est.ymeas   = est.y;

% add solution fields for MHE
Nsim = length(vec.t);
vec.xest       = nan * zeros(grampcMHE.param.Nx, Nsim);
vec.y          = yfct(vec.xest);
vec.yest       = vec.y;
vec.ymeas      = vec.y;
vec.CPUtimeMHE = nan * zeros(1, Nsim);

% init plots and store figure handles
if figNr > 0 && PLOT_PRED
    phpP = grampc_init_plot_pred(grampcMPC, figNr);
    figNr = figNr+1;
end
if figNr > 0 && PLOT_TRAJ
    phpT = grampc_init_plot_sim(vec, figNr);
    figNr = figNr+1;
end
if figNr > 0 && PLOT_STAT
    phpS = grampc_init_plot_stat(vec, grampcMPC, figNr);
    figNr = figNr+1;
end
if figNr > 0 && PLOT_EST
    phpE = init_plot_mhe(est, figNr);
    figNr = figNr+1;
end
if figNr > 0 && PLOT_MHE
    phpM = init_plot_mhe(vec, figNr);
    figNrMHE = figNr+1;
end

%% Combined MHE and MPC loop
i = 1;
while 1
    % true state and output
    x = vec.x(:,i);
    y = yfct(x);
    % measured output
    ymeas = y + randn(1,1)*0.01;

    est.x = [est.x(:,2:end), x];
    est.y = [est.y(:,2:end), y];
    est.ymeas = [est.ymeas(:,2:end), ymeas];
    
    % run moving horizon estimation
    grampcMHE = CmexFiles.grampc_setparam_Cmex(grampcMHE, 't0', vec.t(i));
    grampcMHE.rws.u = est.u;
    % the first two parameters are weights for the cost function
    % therefore, the measured outputs start at index 3
    grampcMHE.userparam(3:end) = est.ymeas;

    if use_discrete
        [grampcMHE, vec.CPUtimeMHE(i)] = CmexFiles.grampc_run_Cmex_MHE_discrete(grampcMHE);
    else
        [grampcMHE, vec.CPUtimeMHE(i)] = CmexFiles.grampc_run_Cmex_MHE_continuous(grampcMHE);
    end

    % estimated state and output
    xest = grampcMHE.rws.p + grampcMHE.rws.x(:,end);
    yest = yfct(xest);

    est.xest = [est.xest(:,2:end), xest];
    est.yest = [est.yest(:,2:end), yest];

    % warmstart MHE by predicting the next initial state
    grampcMHE.rws.pprev = grampcMHE.rws.pprev + grampcMHE.rws.x(:,2);
    grampcMHE.rws.p = grampcMHE.rws.p + grampcMHE.rws.x(:,2);
            
    % run MPC
    grampcMPC = CmexFiles.grampc_setparam_Cmex(grampcMPC, 't0', vec.t(i));
    % grampcMPC = CmexFiles.grampc_setparam_Cmex(grampcMPC, 'x0', x);    % true state
    grampcMPC = CmexFiles.grampc_setparam_Cmex(grampcMPC, 'x0', xest); % estimated state

    if use_discrete
        [grampcMPC, vec.CPUtime(i)] = CmexFiles.grampc_run_Cmex_MPC_discrete(grampcMPC);
    else
        [grampcMPC, vec.CPUtime(i)] = CmexFiles.grampc_run_Cmex_MPC_continuous(grampcMPC);
    end

    % update array of previous controls
    % the first two elements are needed since discrete and erk1 do not use
    % the last element of the control trajectory
    est.u = [est.u(:,2:end-1), grampcMPC.rws.u(:,1), grampcMPC.rws.u(:,2)];

    % save results
    vec = grampc_update_struct_sol(grampcMPC, vec, i);
    vec.xest(:,i)  = xest;
    vec.y(:,i)     = y;
    vec.yest(:,i)  = yest;
    vec.ymeas(:,i) = ymeas;
    
    % print solver status
    printed = CmexFiles.grampc_printstatus_Cmex(grampcMPC.sol.status, 'Error');
    if printed
        fprintf('at simulation time %f.\n --------\n', vec.t(i));    
    end
    
    % check for end of simulation
    if i+1 > length(vec.t)
        break;
    end
    
    % simulate system
    if use_discrete
        vec.x(:,i+1) = CmexFiles.grampc_ffct_Cmex_MPC_discrete(...
            vec.t(i), vec.x(:,i), vec.u(:,i), vec.p(:,i), grampcMPC.param, grampcMPC.userparam);
    else
        [~, xtemp] = ode45(@CmexFiles.grampc_ffct_Cmex_MPC_continuous, ...
            vec.t(i)+[0 double(grampcMPC.param.dt)], vec.x(:,i), odeopt, vec.u(:,i), vec.p(:,i), grampcMPC.param, grampcMPC.userparam);
        vec.x(:,i+1) = xtemp(end,:);
    end
    
    % update iteration counter
    i = i + 1;
    
    % plot data
    if figNr > 0 && (mod(i,PLOT_STEPS) == 0 || i == length(vec.t))
        if PLOT_PRED
            grampc_update_plot_pred(grampcMPC, phpP);
        end
        if PLOT_TRAJ
            grampc_update_plot_sim(vec, phpT);
        end
        if PLOT_STAT
            grampc_update_plot_stat(vec, grampcMPC, phpS);
        end
        if PLOT_EST
            update_plot_mhe(est, phpE);
        end
        if PLOT_MHE
            update_plot_mhe(vec, phpM);
        end
        drawnow
        if PLOT_PAUSE
            pause;
        end
    end
end

end