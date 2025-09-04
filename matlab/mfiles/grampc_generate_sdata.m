function [ grampc_sdata ] = grampc_generate_sdata( grampc )
% This function converts a grampc struct to a struct consisting
% of double arrays which is required for the Simulink interface.
% 
% For this function it is important that the option and parameter
% sequence of the grampc struct matches with the selection in
% grampc_sfun_Cmex.c
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


% Initialization of the output structs' fields
grampc_sdata.intopt = [];
grampc_sdata.intoptidx = 0;
grampc_sdata.numopt = [];
grampc_sdata.numoptidx =0;
grampc_sdata.numparam = [];
grampc_sdata.numparamidx =0;

% Options
optfields = fields(grampc.opt);
for i = 1:length(optfields)
    if strfind(class(grampc.opt.(optfields{i})),'int')
        grampc_sdata.intopt = [grampc_sdata.intopt grampc.opt.(optfields{i})(:)'];
        grampc_sdata.intoptidx = [grampc_sdata.intoptidx grampc_sdata.intoptidx(end)+int32(length(grampc.opt.(optfields{i})))];
    else
        grampc_sdata.numopt = [grampc_sdata.numopt grampc.opt.(optfields{i})(:)'];
        grampc_sdata.numoptidx = [grampc_sdata.numoptidx grampc_sdata.numoptidx(end)+int32(length(grampc.opt.(optfields{i})))];
    end
end
grampc_sdata.intopt = double([length(grampc_sdata.intopt) grampc_sdata.intopt]);
grampc_sdata.intoptidx = double([length(grampc_sdata.intoptidx) grampc_sdata.intoptidx]);
grampc_sdata.numopt = double([length(grampc_sdata.numopt) grampc_sdata.numopt]);
grampc_sdata.numoptidx = double([length(grampc_sdata.numoptidx) grampc_sdata.numoptidx]);

% Parameters
param_fields = fields(grampc.param);
for i = 1:length(param_fields)
    if isempty(strfind(class(grampc.param.(param_fields{i})),'int'))
        grampc_sdata.numparam = [grampc_sdata.numparam grampc.param.(param_fields{i})(:)'];
        grampc_sdata.numparamidx = [grampc_sdata.numparamidx grampc_sdata.numparamidx(end)+int32(length(grampc.param.(param_fields{i})))];
    end
end
grampc_sdata.numparam = double([length(grampc_sdata.numparam) grampc_sdata.numparam]);
grampc_sdata.numparamidx = double([length(grampc_sdata.numparamidx) grampc_sdata.numparamidx]);
end

