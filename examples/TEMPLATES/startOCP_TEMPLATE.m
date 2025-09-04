function [vec,grampc,figNr] = startOCP(figNr,compile,varargin)
% This function runs the TEMPLATE problem. It compiles the c-Code,
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
probfct = 'probfct_TEMPLATE.c';

% plot predicted trajectories
PLOT_PRED = 1;
% plot optimization statistics
PLOT_STAT = 1;
% update plots after N steps
PLOT_STEPS = 10;
% pause after each plot
PLOT_PAUSE = 0;


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
% init GRAMPC, estimate option PenaltyMin and print options and parameters
% Init with options and parameter defines in initData.m
[grampc,~] = initData_TEMPLATE;
% Init with options and parameter from configuration file
% [grampc,~] = initData_TEMPLATE_CONFIG;

CmexFiles.grampc_printopt_Cmex(grampc);
CmexFiles.grampc_printparam_Cmex(grampc);

% init solution structure
MaxMultIter = 100;
vec = grampc_init_struct_sol(grampc, [], MaxMultIter);

% init plots and store figure handles
if figNr > 0 && PLOT_PRED
    phpP = grampc_init_plot_pred(grampc,figNr);
    figNr = figNr+1;
end
if figNr > 0 && PLOT_STAT
    phpS = grampc_init_plot_stat(vec,grampc,figNr);
    figNr = figNr+1;
end


%% OCP solution (Loop over MultIter in matlab to improve debug possibilities)
i = 1;
while i <= MaxMultIter
    % optimize and save results
    [grampc,vec.CPUtime(i)] = CmexFiles.grampc_run_Cmex(grampc);
    vec = grampc_update_struct_sol(grampc, vec, i);
    
    % print solver status
    printed = CmexFiles.grampc_printstatus_Cmex(grampc.sol.status,'Error');
    if printed
        fprintf('at iteration %d.\n --------\n', vec.MultIter(i));
    end
        
    % run convergence check
    converged = 0;
    if grampc.opt.ConvergenceCheck == 1
        converged = 1;
        if grampc.opt.OptimControl == 1 && norm(grampc.rws.u-grampc.rws.uprev)/norm(grampc.rws.u) > grampc.opt.ConvergenceGradientRelTol
            converged = 0;
        elseif grampc.opt.OptimParam == 1 && norm(grampc.rws.p-grampc.rws.pprev)/norm(grampc.rws.p) > grampc.opt.ConvergenceGradientRelTol 
            converged = 0;
        elseif grampc.opt.OptimTime == 1 && norm(grampc.rws.T-grampc.rws.Tprev)/norm(grampc.rws.T) > grampc.opt.ConvergenceGradientRelTol
            converged = 0;
        elseif any(max(abs(grampc.rws.cfct([vec.idx_g,vec.idx_gT],:)),[],2) > grampc.opt.ConstraintsAbsTol([vec.idx_g,vec.idx_gT])')
            converged = 0;
        elseif any(max(grampc.rws.cfct([vec.idx_h,vec.idx_hT],:),[],2) >grampc.opt.ConstraintsAbsTol([vec.idx_h,vec.idx_hT])')
            converged = 0;
        end
    end
    
    % plot data
    if figNr > 0 && (mod(i,PLOT_STEPS) == 0 || i == MaxMultIter || converged)
        if PLOT_PRED
            grampc_update_plot_pred(grampc,phpP);
        end
        if PLOT_STAT
            grampc_update_plot_stat(vec,grampc,phpS);
        end
        drawnow
        if PLOT_PAUSE
            pause;
        end
    end
    
    % break if converged
    if converged
        fprintf('Converged after %i multiplier iterations\n', i);
        break
    end
    
    % update iteration counter
    i = i + 1;
end

fprintf('Computation time: %.3f ms\n', sum(vec.CPUtime(1:i-1)));
fprintf('Cost: %.3f Constraints: %.3f\n', vec.J(2,i-1), vec.Nconstr(i-1));

end
