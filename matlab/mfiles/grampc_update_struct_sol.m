function vec = grampc_update_struct_sol(grampc, vec, i)
% Update solution structure for step i.
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

