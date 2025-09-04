function vec = grampc_update_struct_sol(grampc, vec, i)
% Update solution structure for step i.
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

    % adjoint states
    vec.adj(:,i)        = grampc.rws.adj(:,1);
    % inputs
    vec.u(:,i)          = grampc.sol.unext;
    
    if grampc.param.Ng + grampc.param.Nh > 0
        % time-dependent constraints
        vec.constr(:,i) = [grampc.rws.cfct(vec.idx_g,1); max(0,grampc.rws.cfct(vec.idx_h,1))];
        % time-dependent  multipliers
        vec.mult(:,i)   = grampc.rws.mult([vec.idx_g,vec.idx_h],1);
        % time-dependent penalties
        vec.pen(:,i)    = grampc.rws.pen([vec.idx_g,vec.idx_h],1);
    end
    
    % augmented and original cost
    vec.J(1,i)          = grampc.sol.J(1);
    if grampc.param.Nc > 0
        vec.J(2,i)      = grampc.sol.J(2);
    end
    % computation time
%     vec.CPUtime(i)    = 0;
    
    if grampc.param.Nc > 0
        % norm of all constraints
        vec.Nconstr(i)  = grampc.sol.cfct;
        % norm of all penalties
        vec.Npen(i)     = grampc.sol.pen;
    end
    
    % number of gradient iterations per multiplier iteration
    vec.iter(:,i)       = grampc.sol.iter;
    % parameters of adaptive or explicit line search
    if ~isempty(grampc.rws.lsAdapt)
        vec.lsAdapt(:,i) = grampc.rws.lsAdapt((end-7):(end-4));
    else
        vec.lsExpl(i)   = grampc.rws.lsExplicit(3);
    end
        
    % parameters
    vec.p(:,i)          = grampc.sol.pnext;
    % final time
    vec.T(:,i)          = grampc.sol.Tnext;
    
    if grampc.param.NgT + grampc.param.NhT > 0
        % terminal constraints
        vec.constrT(:,i) = [grampc.rws.cfct(vec.idx_gT,end); max(0,grampc.rws.cfct(vec.idx_hT,end))];
        % terminal multipliers
        vec.multT(:,i)   = grampc.rws.mult([vec.idx_gT,vec.idx_hT],end);
        % terminal penalties
        vec.penT(:,i)    = grampc.rws.pen([vec.idx_gT,vec.idx_hT],end);
    end
    
    % status flags
    vec.status(i)       = grampc.sol.status;
end

