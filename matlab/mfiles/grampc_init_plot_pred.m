function ph_out = grampc_init_plot_pred(grampc,figNr,idx,decimation)
% This function initializes the prediction plot of GRAMPC. A plot
% handle is returned to be passed to grampc_update_plot_pred
% for successive update of the prediction plot over the MPC
% runtime. Call this function once before the simulation loop.
%
% Input arguments are:
% 1) grampc - GRAMPC struct
% 2) figNr - number of the figure
% 3) idx - index struct that specifies the quantities to be plotted (optional)
% 4) decimation - factor for downsampling the plot (optional)
%
% Output argument is:
% 1) ph_out - plot handles which serve as input for grampc_update_plot_pred funtion
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

% Check input arguments
if nargin < 2 || isempty(figNr)
    figNr = 1;
end
if nargin < 3 || isempty(idx)
    idx.States = 1:size(grampc.rws.x,1);
    idx.Controls = 1:size(grampc.rws.u,1);
    idx.Constraints = 1:size(grampc.rws.cfct,1);
end
if nargin < 4 || isempty(decimation)
    decimation = 1;
end

% Specify terminal constraints
idx.TConstraints = idx.Constraints(idx.Constraints>(grampc.param.Ng+grampc.param.Nh));
idx.gConstraints = idx.Constraints(idx.Constraints<=(grampc.param.Ng+grampc.param.Nh));

% If there are constraints two plot rows are needed
if ~isempty(idx.Constraints)
    sprownr = 2;
else
    sprownr = 1;
end

% If new handle needed clear figure
figure(figNr)
set(figNr,'Name','Prediction plot of GRAMPC')
clf(figNr)

% Determine the number of columns
if grampc.opt.LineSearchType == 0
    spcolnr = 4;
else
    spcolnr = 3;
end

add_axis =[];

% plot states
subplot_States   = subplot(sprownr,spcolnr,1,'Parent',figNr);
ph_out.s1 = plot(grampc.rws.t(1:decimation:end),grampc.rws.x(idx.States,1:decimation:end),'Parent',subplot_States);
title(subplot_States,'Predicted states');
xlabel(subplot_States,'Time');
%     legendcell = cell(length(idx.States),1);
%     for i = 1:length(idx.States)
%         legendcell{i} = sprintf('x_{%d}',idx.States(i));
%     end
%     legend(legendcell)

% plot adjoint states
subplot_adjStates   = subplot(sprownr,spcolnr,2,'Parent',figNr);
ph_out.s2 = plot(grampc.rws.t(1:decimation:end),grampc.rws.adj(idx.States,1:decimation:end),'Parent',subplot_adjStates);
title(subplot_adjStates,'Predicted adjoint states');
xlabel(subplot_adjStates,'Time');

% plot controls
subplot_Controls = subplot(sprownr,spcolnr,3,'Parent',figNr);
ph_out.s3 = plot(grampc.rws.t(1:decimation:end),grampc.rws.u(idx.Controls,1:decimation:end),'Parent',subplot_Controls);
title(subplot_Controls,'Predicted controls');
xlabel(subplot_Controls,'Time');
%     legendcell = cell(length(idx.Controls),1);
%     for i = 1:length(idx.Controls)
%         legendcell{i} = sprintf('u_{%d}',idx.Controls(i));
%     end
%     legend(legendcell)


% line search
if grampc.opt.LineSearchType == 0
    subplot_ls = subplot(sprownr,spcolnr,4,'Parent',figNr);
    title(subplot_ls,'Line search iterations');
    xlabel(subplot_ls,'stepsize');
    ylabel(subplot_ls,'costs')
    
    hold(subplot_ls,'on')
    for i = 1:grampc.opt.MaxGradIter
        ph_out.s4a(i) = plot(1:19,nan*(1:19),'Color',[0 0.4470 0.7410],'Parent',subplot_ls);
        ph_out.s4b(i) = plot(1:3,nan*(1:3),'Marker','.','MarkerSize',17,'MarkerEdgeColor',[0 0.4470 0.7410],'LineStyle','none','Parent',subplot_ls);
        ph_out.s4c(i) = plot(1,nan,'Marker','.','MarkerSize',15,'MarkerEdgeColor',[0.8500 0.3250 0.0980],'Parent',subplot_ls);
    end
    hold(subplot_ls,'off')
end

% Constraints
if ~isempty(idx.Constraints)
    subplot_cfct = subplot(2,3,4,'Parent',figNr);
    plot(grampc.rws.t([1 end]),[0 0],'k--','Parent',subplot_cfct), hold on
    if ~isempty(idx.gConstraints)
        ph_out.s5a = plot(grampc.rws.t(1:decimation:end),[grampc.rws.cfct(1:grampc.param.Ng,1:decimation:end);max(grampc.rws.cfct(grampc.param.Ng + (1:grampc.param.Nh),1:decimation:end),0)],'Parent',subplot_cfct);
    end
    if ~isempty(idx.TConstraints)
        hold on,
        ph_out.s5b = plot(grampc.rws.t(end),[grampc.rws.cfct(grampc.param.Ng + grampc.param.Nh + (1:grampc.param.NgT),end);max(grampc.rws.cfct(grampc.param.Ng + grampc.param.Nh + grampc.param.NgT + (1:grampc.param.NhT),end),0)],'o','Parent',subplot_cfct);
        hold off,
    end
    hold off
    title(subplot_cfct,'Predicted constraints');
    xlabel(subplot_cfct,'Time');
%     legendcell = cell(length(idx.Constraints),1);
%     for i = 1:length(idx.Constraints)
%         if idx.Constraints(i) <= grampc.param.Ng
%             legendcell{i} = sprintf('g_{%d}',idx.Constraints(i));
%         elseif idx.Constraints(i) <= grampc.param.Ng+grampc.param.Nh
%             legendcell{i} = sprintf('h_{%d}',idx.Constraints(i)-grampc.param.Ng);
%         elseif idx.Constraints(i) <= grampc.param.Ng+grampc.param.Nh+grampc.param.NgT
%             legendcell{i} = sprintf('g_{T_{%d}}',idx.Constraints(i)-grampc.param.Ng-grampc.param.Nh);
%         else
%             legendcell{i} = sprintf('h_{T_{%d}}',idx.Constraints(i)-grampc.param.Ng-grampc.param.Nh-grampc.param.NgT);
%         end
%     end
%     legend(legendcell)
    
    % mult
    subplot_mult      = subplot(2,3,5,'Parent',figNr);
    if ~isempty(idx.gConstraints)
        ph_out.s6a = plot(grampc.rws.t(1:decimation:end),grampc.rws.mult(idx.gConstraints,1:decimation:end),'Parent',subplot_mult);
    end
    if ~isempty(idx.TConstraints)
        hold on,
        ph_out.s6b = plot(grampc.rws.t(end),grampc.rws.mult(idx.TConstraints,end),'o','Parent',subplot_mult);
        hold off,
    end
    title(subplot_mult,'Predicted Lagrange multipliers');
    xlabel(subplot_mult,'Time');
    
    % pen
    subplot_pen      = subplot(2,3,6,'Parent',figNr);
    if ~isempty(idx.gConstraints)
        ph_out.s7a = plot(grampc.rws.t(1:decimation:end),grampc.rws.pen(idx.gConstraints,1:decimation:end),'Parent',subplot_pen);
    end
    if ~isempty(idx.TConstraints)
        hold on,
        ph_out.s7b = plot(grampc.rws.t(end),grampc.rws.pen(idx.TConstraints,end),'o','Parent',subplot_pen);
        hold off,
    end
    title(subplot_pen,'Predicted penalty parameters');
    xlabel(subplot_pen,'Time');
    
    add_axis = [subplot_cfct,subplot_mult,subplot_pen];
end

linkaxes([subplot_States,subplot_adjStates,subplot_Controls,add_axis],'x')
xlim([grampc.rws.t(1),grampc.rws.t(end)])
end