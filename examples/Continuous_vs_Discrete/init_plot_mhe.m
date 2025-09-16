function ph_out = init_plot_mhe(vec, figNr)
% This function initializes the MHE plot. A plot handle is returned to be
% passed to update_plot_mhe for successive updates.
% Call this function once before the simulation loop.
%
% Input arguments are:
% 1) vec - struct with time, inputs, states and outputs
% 2) figNr - number of the figure
%
% Output argument is:
% 1) ph_out - plot handles which serve as input for update_plot_mhe function
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

% Check input arguments
if nargin < 2 || isempty(figNr)
    figNr = 1;
end

% If new handle needed clear figure
figure(figNr)
set(figNr, 'Name', 'Moving horizon estimation')
clf(figNr)

% Inputs
subplot_inputs = subplot(1, 3, 1, 'Parent', figNr);
ph_out.s1 = plot(vec.t, vec.u, 'Parent', subplot_inputs);
title(subplot_inputs, 'Inputs');
xlabel(subplot_inputs, 'Time');

% States
subplot_states = subplot(1, 3, 2, 'Parent', figNr);
ph_out.s2 = plot(vec.t, vec.x, 'Parent', subplot_states, 'DisplayName', 'True');
hold on
set(gca, 'ColorOrderIndex', 1);
ph_out.s3 = plot(vec.t, vec.xest, '--', 'Parent', subplot_states, 'DisplayName', 'Estimated');
hold off
title(subplot_states, 'States');
xlabel(subplot_states, 'Time');
legend();

% Outputs
subplot_outputs = subplot(1, 3, 3, 'Parent', figNr);
ph_out.s4 = plot(vec.t, vec.y, 'Parent', subplot_outputs, 'DisplayName', 'True');
hold on
set(gca, 'ColorOrderIndex', 1);
ph_out.s5 = plot(vec.t, vec.yest, '--', 'Parent', subplot_outputs, 'DisplayName', 'Estimated');
set(gca, 'ColorOrderIndex', 1);
ph_out.s6 = plot(vec.t, vec.ymeas, 'x', 'Parent', subplot_outputs, 'DisplayName', 'Measured');
hold off
title(subplot_outputs, 'Outputs');
xlabel(subplot_outputs, 'Time');
legend();

linkaxes([subplot_inputs, subplot_states, subplot_outputs], 'x');
xlim([vec.t(1) vec.t(end)]);

end
