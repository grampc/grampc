function ph_out = grampc_init_plot_stat(vec,grampc,figNr,decimation)
% This function initializes the statistics plot of GRAMPC. A plot
% handle is returned that is to be passed to grampc_update_plot_stat
% for successive update of the simulation plot over the MPC
% runtime. Call this function once before the simulation loop.
%
% Input arguments are:
% 1) vec - solution struct
% 2) grampc - GRAMPC struct
% 3) figNr - number of the figure
% 4) decimation - factor for downsampling the plot
%
% Output argument is:
% 1) ph_out - plot handles which serve as input for grampc_update_plot_stat funtion
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
if nargin < 3 || isempty(figNr)
    figNr = 3;
end
if nargin < 4 || isempty(decimation)
    decimation = 1;
end

% Determine number of terminal constraints
NTerminalConstraints = (grampc.param.NgT + grampc.param.NhT);

% Create figure
figure(figNr)
set(figNr,'Name','Statistics plot of GRAMPC')
clf(figNr)

% Depending whether an OCP is solved, different xdata and labels must be used
if isfield(vec,'MultIter')
    xlab = 'Iterations';
    xlabct = 'Iterations';
    xdata = vec.MultIter;
else
    xlab = 'Time';
    xlabct = 'Iteration \cdot dt';
    xdata = vec.t;
end

% Depnding on the options and constraints different column numbers are needed
spnr = 3 + 2*double(grampc.param.Nc ~= 0)+ double(grampc.opt.ConvergenceCheck) + double(grampc.opt.OptimTime == 1) + double(grampc.opt.OptimParam == 1);
switch spnr
    case {2,3}
        spcolnr(1)= spnr;
        spcolnr(2) = 0;
    case {4,6,8}
        spcolnr(1) = spnr/2;
        spcolnr(2) = spnr/2;
    case 5
        spcolnr(1) = 3;
        spcolnr(2) = 2;
    case 7
        spcolnr(1) = 4;
        spcolnr(2) = 3;
end

% Depending on the options and constraints three, two or one rows are needed
if (grampc.param.NgT + grampc.param.NhT) ~=0
    sprownr = 3;
elseif spcolnr(2) ~= 0
    sprownr = 2;
else
    sprownr = 1;
end

add_axis = [];

% Plot cost
subplot_Cost = subplot(sprownr,spcolnr(1),1,'Parent',figNr);
ph_out.s1(1) = plot(xdata(1:decimation:end),vec.J(1,1:decimation:end),'Parent',subplot_Cost);
if(size(vec.J)>1)
    hold on,
    ph_out.s1(2) = plot(xdata(1:decimation:end),vec.J(1,1:decimation:end),'--','Parent',subplot_Cost);
    hold off,
    legend('Original costs','Augmented costs')
end
title(subplot_Cost,'Costs');
xlabel(subplot_Cost,xlab);

% Computation time
subplot_comp_time = subplot(sprownr,spcolnr(1),2,'Parent',figNr);
ph_out.s2 = plot(xdata(1:decimation:end),vec.CPUtime(1:decimation:end),'Parent',subplot_comp_time);
title(subplot_comp_time,'Computation time (ms)');
xlabel(subplot_comp_time,xlabct);


% plot iterations
if grampc.opt.ConvergenceCheck
    subplot_iter = subplot(sprownr,spcolnr(1),3,'Parent',figNr);
    ph_out.s3 = stem(xdata(1:decimation:end)',vec.iter(:,1:decimation:end)','Parent',subplot_iter);
    title(subplot_iter,'Gradient iterations');
    xlabel(subplot_iter,xlab);
    add_axis = [add_axis subplot_iter];
end

% plot adaptive LineSearch interval
    plot_nr = 3+(grampc.opt.ConvergenceCheck==1);
    if(plot_nr)>spcolnr(1) && grampc.param.Nc ~= 0
        plot_nr = plot_nr+2;
    end
    subplot_als = subplot(sprownr,spcolnr(ceil(plot_nr/spcolnr(1))),plot_nr+diff(spcolnr)*(plot_nr>spcolnr(1)),'Parent',figNr);
    if grampc.opt.LineSearchType == 0
        ph_out.s4 = semilogy(xdata(1:decimation:end),vec.lsAdapt(:,1:decimation:end),'Parent',subplot_als);
        legend('Lower interval bound','Interval center','Upper interval bound','Used stepsize')
    else
        ph_out.s4 = semilogy(xdata(1:decimation:end),vec.lsExpl(1:decimation:end),'Parent',subplot_als);
    end
    title(subplot_als,'Line search step size');
    xlabel(subplot_als,xlab);
    add_axis = [add_axis subplot_als];

% plot Thor
if grampc.opt.OptimTime == 1
    plot_nr = 4+(grampc.opt.ConvergenceCheck==1);
    if(plot_nr)>spcolnr(1) && grampc.param.Nc ~= 0
        plot_nr = plot_nr+2;
    end
    subplot_Thor = subplot(sprownr,spcolnr(ceil(plot_nr/spcolnr(1))),plot_nr+diff(spcolnr)*(plot_nr>spcolnr(1)),'Parent',figNr);
    ph_out.s5 = plot(xdata(1:decimation:end),vec.T(1:decimation:end),'Parent',subplot_Thor);
    title(subplot_Thor,'Prediction horizon');
    xlabel(subplot_Thor,xlab);
    add_axis = [add_axis subplot_Thor];
end

% plot param
if grampc.opt.OptimParam == 1
    plot_nr =  4+(grampc.opt.ConvergenceCheck==1)+(grampc.opt.OptimTime == 1);
    if(plot_nr)>spcolnr(1) && grampc.param.Nc ~= 0
        plot_nr = plot_nr+2;
    end
    subplot_param = subplot(sprownr,spcolnr(ceil(plot_nr/spcolnr(1))),plot_nr+diff(spcolnr)*(plot_nr>spcolnr(1)),'Parent',figNr);
    ph_out.s6 = plot(xdata(1:decimation:end),vec.p(:,1:decimation:end),'Parent',subplot_param);
    title(subplot_param,'Parameter');
    xlabel(subplot_param,xlab);
    add_axis = [add_axis subplot_param];
end

if grampc.param.Nc ~= 0
    % plot cNorm
    subplot_cNorm = subplot(sprownr,spcolnr(2),spcolnr(2)+1,'Parent',figNr);
    ph_out.s7 = plot(xdata(1:decimation:end),vec.Nconstr(1:decimation:end),'Parent',subplot_cNorm);
    title(subplot_cNorm,'Norm of constraints over horizon');
    xlabel(subplot_cNorm,xlab);
    add_axis = [add_axis subplot_cNorm];
    
    % plot penNorm
    subplot_penNorm = subplot(sprownr,spcolnr(2),spcolnr(2)+2,'Parent',figNr);
    ph_out.s8 = plot(xdata(1:decimation:end),vec.Npen(1:decimation:end),'Parent',subplot_penNorm);
    title(subplot_penNorm,'Norm of penalty parameters over horizon');
    xlabel(subplot_penNorm,xlab);
    add_axis = [add_axis subplot_penNorm];
    
    if  NTerminalConstraints~=0
        % Terminal constraints
        subplot_Tcfct = subplot(3,3,7,'Parent',figNr);
        ph_out.s9 = plot(xdata(1:decimation:end),vec.constrT(:,1:decimation:end),'Parent',subplot_Tcfct);
        title(subplot_Tcfct,'Terminal constraints');
        xlabel(subplot_Tcfct,xlab);
        
        % Terminal Lagrangian multiplier
        subplot_Tmult = subplot(3,3,8,'Parent',figNr);
        ph_out.s10 = plot(xdata(1:decimation:end),vec.multT(:,1:decimation:end),'Parent',subplot_Tmult);
        title(subplot_Tmult,'Terminal Lagrangian multipliers');
        xlabel(subplot_Tmult,xlab);
        
        % Terminal Penalties
        subplot_Tpen = subplot(3,3,9,'Parent',figNr);
        ph_out.s11 = plot(xdata(1:decimation:end),vec.penT(:,1:decimation:end),'Parent',subplot_Tpen);
        title(subplot_Tpen,'Terminal penalty parameters');
        xlabel(subplot_Tpen,xlab);
        
        add_axis = [add_axis subplot_Tcfct,subplot_Tmult,subplot_Tpen];
    end
end

linkaxes([subplot_Cost,subplot_comp_time,add_axis],'x')
if isfield(vec,'MultIter')
    xlim([xdata(1) min(xdata(end),50)])
else
    xlim([xdata(1) xdata(end)])
end
