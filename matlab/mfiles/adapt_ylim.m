% This function adapts ylim
%
% Input arguments are:
% 1) data - data of the plot
% 2) line_handle - line handle to the desired plot
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

function adapt_ylim(data,line_handle)
ymax = max(max(data));
ymin = min(min(data));
ycenter = (ymax+ymin)/2;
yrange = (ymax-ymin)*1.1;
if yrange == 0
    if ycenter == 0
        yrange = 2;
    else
        yrange = abs(ycenter)/2;
    end
end
ylim(get(line_handle(1),'Parent'),[-yrange yrange]/2+ycenter);
end