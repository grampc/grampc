function [] = grampc_make_toolbox(grampc_root_path, varargin)
% Compiles the GRAMPC toolbox.
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

% check input arguments
if nargin < 1
    error('grampc_make_toolbox: not enough arguments');
end

if grampc_root_path(end) ~= '/'
    grampc_root_path = [grampc_root_path '/'];
end

% switch to grampc/matlab directory and compile toolbox
curdir = cd;
cd([grampc_root_path 'matlab/']);
try
    make(varargin{:});
catch err
    cd(curdir);
    rethrow(err)
end
cd(curdir);

end