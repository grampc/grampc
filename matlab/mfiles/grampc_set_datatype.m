function [ convdata ] = grampc_set_datatype( dataType, data )
% This function casts data to the type specified by dataType. This can
% be either 'single' (='float'), 'double', 'int64, 'int32' (='int'), 'int16' or 'int8'.
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

switch(dataType)
    case {'single','float'}
        convdata = single(data);
    case 'double'
        convdata = double(data);
    case 'int64'
        convdata = int64(data);
    case {'int32','int'}
        convdata = int32(data);
    case 'int16'
        convdata = int16(data);
    case 'int8'
        convdata = int8(data);
    otherwise
        error(['setDataType: unknown data tye: ' dataType]);
end
end

