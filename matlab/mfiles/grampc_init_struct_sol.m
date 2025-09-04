function vec = grampc_init_struct_sol(grampc, Tsim, MaxMultIter)
% Create solution structure either for duration Tsim or for Nsim steps.
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

