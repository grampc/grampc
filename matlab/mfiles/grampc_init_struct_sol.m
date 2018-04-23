function vec = grampc_init_struct_sol(grampc, Tsim, MaxMultIter)
% Create solution structure either for duration Tsim or for Nsim steps.
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
    
    if nargin < 2
        error('grampc_init_struct_sol: not enough arguments');
    end
    
    % index sets for constraints
    vec.idx_g        = 1:grampc.param.Ng;
    vec.idx_h        = grampc.param.Ng + (1:grampc.param.Nh);
    vec.idx_gT       = grampc.param.Ng + grampc.param.Nh + (1:grampc.param.NgT);
    vec.idx_hT       = grampc.param.Ng + grampc.param.Nh + grampc.param.NgT + (1:grampc.param.NhT);

    % time
    if nargin > 1 && ~isempty(Tsim)
        vec.t        = double(0 : grampc.param.dt : Tsim);
    end
    
    % use length of vec.t or (if given) MaxMultIter for number of steps
    if nargin < 3
        Nsteps       = length(vec.t);
    else
        Nsteps       = MaxMultIter;
        vec.MultIter = 1 : MaxMultIter;
    end
    
    % states, set initial state
    vec.x           = nan * zeros(grampc.param.Nx, Nsteps);
    vec.x(:,1)      = grampc.param.x0;
    % adjoint states
    vec.adj         = nan * zeros(grampc.param.Nx, Nsteps);
    % inputs
    vec.u           = nan * zeros(grampc.param.Nu, Nsteps);
    
    if grampc.param.Ng + grampc.param.Nh > 0
        % time-dependent constraints
        vec.constr  = nan * zeros(grampc.param.Ng+grampc.param.Nh, Nsteps);
        % time-dependent multipliers
        vec.mult    = nan * zeros(grampc.param.Ng+grampc.param.Nh, Nsteps);
        % time-dependent penalties
        vec.pen     = nan * zeros(grampc.param.Ng+grampc.param.Nh, Nsteps);
    end
    
    % augmented and original cost
    if grampc.param.Nc > 0
        vec.J       = nan * zeros(2, Nsteps);
    else
        vec.J       = nan * zeros(1, Nsteps);
    end
    % computation time
    vec.CPUtime     = nan * zeros(1, Nsteps);
    
    if grampc.param.Nc > 0
        % norm of all constraints
        vec.Nconstr = nan * zeros(1, Nsteps);
        % norm of all penalties
        vec.Npen    = nan * zeros(1, Nsteps);
    end
    
    % number of gradient iterations per multiplier iteration
    vec.iter        = nan * zeros(grampc.opt.MaxMultIter, Nsteps);
    % parameters of adaptive or explicit line search
    if ~isempty(grampc.rws.lsAdapt)
        vec.lsAdapt = nan * zeros(4, Nsteps);
    else
        vec.lsExpl  = nan * zeros(1, Nsteps);
    end
    
    % parameters
    vec.p           = nan * zeros(grampc.param.Np, Nsteps);
    % final time
    vec.T           = nan * zeros(1, Nsteps);
    
    if grampc.param.NgT + grampc.param.NhT > 0
        % terminal constraints
        vec.constrT = nan * zeros(grampc.param.NgT+grampc.param.NhT, Nsteps);
        % terminal multipliers
        vec.multT   = nan * zeros(grampc.param.NgT+grampc.param.NhT, Nsteps);
        % terminal penalties
        vec.penT    = nan * zeros(grampc.param.NgT+grampc.param.NhT, Nsteps);
    end
    
    % status flags
    vec.status      = zeros(1, Nsteps);
end

