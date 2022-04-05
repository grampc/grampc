function [ convdata ] = grampc_set_datatype( dataType, data )
% This function casts data to the type specified by dataType. This can
% be either 'single' (='float'), 'double', 'int64, 'int32' (='int'), 'int16' or 'int8'.
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

