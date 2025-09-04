function ph_out = grampc_init_plot_sim(vec,figNr,idx,decimation)
% This function initializes the simulation plot of GRAMPC. A plot
% handle is returned to be passed to grampc_update_plot_sol
% for successive update of the simulation plot over the MPC
% runtime. Call this function once before the simulation loop.
%
% Input arguments are:
% 1) vec - solution struct
% 2) figNr - number of the figure
% 3) idx - index struct that specifies the quantities to be plotted (optional)
% 4) decimation - factor for downsampling the plot (optional)
%
% Output argument is:
% 1) ph_out - plot handles which serve as input for grampc_update_plot_sim funtion
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

if nargin < 2 || isempty(figNr)
    figNr = 2;
end
if nargin < 3 || isempty(idx)
    idx.States = 1:size(vec.x,1);
    idx.Controls = 1:size(vec.u,1);
    if isfield(vec,'constr')
        idx.Constraints = 1:size(vec.constr,1);
    else
        idx.Constraints = [];
    end
end
if nargin < 4 || isempty(decimation)
    decimation = 1;
end

% Create figure
figure(figNr)
set(figNr,'Name','Simulation plot of GRAMPC')
clf(figNr)

% If there are constraints two plot rows are needed
if ~isempty(idx.Constraints)
    sprownr = 2;
else
    add_axis = [];
    sprownr = 1;
end

% plot states
subplot_States = subplot(sprownr,3,1,'Parent',figNr);
ph_out.s1 = plot(vec.t(1:decimation:end),vec.x(idx.States,(1:decimation:end)),'Parent',subplot_States);
title(subplot_States,'States');
xlabel(subplot_States,'Time');

% plot adjoint states
subplot_adjStates = subplot(sprownr,3,2,'Parent',figNr);
ph_out.s2 = plot(vec.t(1:decimation:end),vec.adj(idx.States,(1:decimation:end)),'Parent',subplot_adjStates);
title(subplot_adjStates,'Adjoint states');
xlabel(subplot_adjStates,'Time');

% plot controls
subplot_Controls = subplot(sprownr,3,3,'Parent',figNr);
ph_out.s3 = plot(vec.t(1:decimation:end),vec.u(idx.Controls,(1:decimation:end)),'Parent',subplot_Controls);
title(subplot_Controls,'Controls');
xlabel(subplot_Controls,'Time');

if ~isempty(idx.Constraints)
    % Constraints
    subplot_cfct = subplot(2,3,4,'Parent',figNr);
    plot(vec.t([1 end]),[0 0],'k--','Parent',subplot_cfct),hold on
    ph_out.s4 = plot(vec.t(1:decimation:end),vec.constr(idx.Constraints,(1:decimation:end)),'Parent',subplot_cfct);
    hold off
    title(subplot_cfct,'Constraints');
    xlabel(subplot_cfct,'Time');
    
    % mult
    subplot_mult = subplot(2,3,5,'Parent',figNr);
    ph_out.s5 = plot(vec.t(1:decimation:end),vec.mult(idx.Constraints,(1:decimation:end)),'Parent',subplot_mult);
    title(subplot_mult,'Lagrange multipliers');
    xlabel(subplot_mult,'Time');
    
    % pen
    subplot_pen = subplot(2,3,6,'Parent',figNr);
    ph_out.s6 = plot(vec.t(1:decimation:end),vec.pen(idx.Constraints,(1:decimation:end)),'Parent',subplot_pen);
    title(subplot_pen,'Penalty parameters');
    xlabel(subplot_pen,'Time');
    
    add_axis = [subplot_cfct,subplot_mult,subplot_pen];
end


linkaxes([subplot_States,subplot_adjStates,subplot_Controls,add_axis],'x')
xlim([vec.t(1) vec.t(end)])
end
