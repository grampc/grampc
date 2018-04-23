function grampc_update_plot_stat(vec,grampc,ph,decimation)
% This function updates the statistics plot of GRAMPC. The figure
% specified via the plot handle ph must have been generated by the
% function grampc_init_plot_stat. 
%
% Input arguments are:
% 1) vec - solution struct
% 2) grampc - GRAMPC struct
% 3) ph - plot handles generated by the grampc_initStatisticsPlot function
% 4) decimation - factor for downsampling the plot
%
%
% This file is part of GRAMPC - (https://sourceforge.net/projects/grampc/)
%
% GRAMPC -- A software framework for embedded nonlinear model predictive
% control using a gradient-based augmented Lagrangian approach
%
% Copyright (C) 2014-2018 by Tobias Englert, Knut Graichen, Felix Mesmer, 
% Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0). 
% Developed at the Institute of Measurement, Control, and Microtechnology, 
% Ulm University. All rights reserved.
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

% Check input arguments
if nargin < 4 || isempty(decimation)
    decimation = 1;
end

% Determine number of terminal constraints
NTerminalConstraints = (grampc.param.NgT + grampc.param.NhT);

% Update handles
% plot cost
for i = 1:size(vec.J,1)
    set(ph.s1(i),'YData',vec.J(i,1:decimation:end))
end

% Computation time
set(ph.s2,'YData',vec.CPUtime(1:decimation:end))

% plot iterations
if grampc.opt.ConvergenceCheck
    for i = 1:size(vec.iter,1)
        set(ph.s3(i),'YData',vec.iter(i,1:decimation:end))
    end
end

% plot adaptive LineSearch interval
if grampc.opt.LineSearchType == 0
    for i = 1:4
    set(ph.s4(i),'YData',vec.lsAdapt(i,1:decimation:end))
    end
else
    set(ph.s4,'YData',vec.lsExpl(1:decimation:end))
end

% plot Thor
if grampc.opt.OptimTime == 1
    set(ph.s5,'YData',vec.T(1:decimation:end))
end

% plot param
if grampc.opt.OptimParam == 1
    for i = 1:size(vec.p,1)
        set(ph.s6(i),'YData',vec.p(i,1:decimation:end))
    end
end

if grampc.param.Nc ~= 0
    % plot cNorm
    set(ph.s7,'YData',vec.Nconstr(1:decimation:end))
    
    % plot penNorm
    set(ph.s8,'YData',vec.Npen(1:decimation:end))
    
    for i = 1:NTerminalConstraints
        % Terminal constraints
        set(ph.s9(i),'YData',vec.constrT(i,1:decimation:end))
        
        % Terminal Lagrangian multiplier
        set(ph.s10(i),'YData',vec.multT(i,1:decimation:end))
        
        % Terminal Penalties
        set(ph.s11(i),'YData',vec.penT(i,1:decimation:end))
    end
end

% Adjust x axis in OCP plots
if isfield(vec,'MultIter')
    subplot_comptime = get(ph.s2,'Parent');
    iter = ceil(find(isnan(vec.CPUtime),1)/50)*50;
    if(isempty(iter))
        iter = vec.MultIter(end);
    end
    xlim(subplot_comptime,[vec.MultIter(1), iter])
end

end

