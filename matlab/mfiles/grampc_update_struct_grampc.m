function [ grampc ] = grampc_update_struct_grampc(grampc,user)
% This function sets all parameters and options specified by the user struct
% to the grampc structure ensuring the correct data types.
%
% The user struct can contain a param and an opt struct. The fields of 
% the opt and param struct can be any GRAMPC option or parameter.
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

% Get the floating point datatype (single or double)
typeRNum = class(grampc.opt.LineSearchMax);
% Update the datatype of the userparams, as they are floating point values
grampc.userparam = grampc_set_datatype(typeRNum,grampc.userparam);

% Set the options and parameters with correct data type
try
    param_fields = fields(user.param);
catch
    param_fields = [];
end
try
opt_fields = fields(user.opt);
catch
    opt_fields = [];
end
for i = 1:length(param_fields)
    if ~isempty(user.param.(param_fields{i}))
        user.param.(param_fields{i}) = double(user.param.(param_fields{i}));
        grampc = CmexFiles.grampc_setparam_Cmex(grampc,param_fields{i},user.param.(param_fields{i}));
    end
end
for i = 1:length(opt_fields)
    if ~isempty(user.opt.(opt_fields{i}))
        if isnumeric(user.opt.(opt_fields{i}))
            user.opt.(opt_fields{i}) = double(user.opt.(opt_fields{i}));
        end
        grampc = CmexFiles.grampc_setopt_Cmex(grampc,opt_fields{i},user.opt.(opt_fields{i}));
    end
end
end

