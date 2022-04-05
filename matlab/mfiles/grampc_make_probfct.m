function [] = grampc_make_probfct(grampc_root_path, probfct, varargin)
% Compiles the GRAMPC problem in MATLAB/Simulink.
%
% Input arguments are:
% 1) grampc_root_path - relative or absolute path to GRAMPC
% 2) probfct - name of problem description C-file
% 3 - end) optional arguments
%     'debug':   Creates the executables with additional information for
%                use in debugging. See also MEX flag -g.
%     'verbose': Create the executables with with verbose option, i.e. 
%                print commands during build. See also MEX flag -v.
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
if nargin < 2
    error('grampc_make_probfct: not enough arguments');
end

if grampc_root_path(end) ~= '/'
    grampc_root_path = [grampc_root_path '/'];
end

% path to resulting Cmex-files
RESPATH = '+CmexFiles';
% path to include folder
CHEADERPATH = [grampc_root_path 'include'];
% path to bin folder
BINPATH = [grampc_root_path 'matlab/bin'];

% create destination directory
if ~exist(RESPATH, 'dir')
    mkdir(RESPATH);
end

% set flags
DEBUG   = '';
VERBOSE = '';
for i = 1 : nargin-2
    if ~isempty(varargin{i})
        switch(varargin{i})
            case 'debug'
                DEBUG = ' -g';
            case 'verbose'
                VERBOSE = ' -v';
        end
    end
end

% set file extension for objects
if strcmp(computer,'PCWIN') || strcmp(computer,'PCWIN64')
    OBJEXT = 'obj';
else
    OBJEXT = 'o';
end

% names of object files
ObjFiles = {['grampc_run.',OBJEXT],...
    ['grampc_alloc.',OBJEXT],...
    ['grampc_fixedsize.',OBJEXT],...
    ['grampc_init.',OBJEXT],...
    ['grampc_mess.',OBJEXT],...
    ['grampc_setopt.',OBJEXT],...
    ['grampc_setparam.',OBJEXT],...
    ['grampc_util.',OBJEXT],...
    ['euler1.',OBJEXT],...
    ['eulermod2.',OBJEXT],...
    ['heun2.',OBJEXT],...
    ['ruku45.',OBJEXT],...
    ['rodas.',OBJEXT],...
    ['trapezodial.',OBJEXT],...
    ['simpson.',OBJEXT]};

% names of executables
EXE = {'grampc_init_Cmex',...
    'grampc_setopt_Cmex',...
    'grampc_setparam_Cmex',...
    'grampc_run_Cmex',...
    'grampc_run_Sfct',...
    'grampc_ghfct_Cmex',...
    'grampc_gThTfct_Cmex',...
    'grampc_lfct_Cmex',...
    'grampc_Vfct_Cmex',...
    'grampc_ffct_Cmex',...
    'grampc_ffct_Sfct',...
    'grampc_estim_penmin_Cmex',...
    'grampc_printopt_Cmex',...
    'grampc_printparam_Cmex',...
    'grampc_printstatus_Cmex'};

% clean output files 
delete([RESPATH '/probfct.',OBJEXT]);
for i = 1:length(EXE)
    if ~strcmp(EXE{i},'grampc_run_Sfct') || ~strcmp(EXE{i},'grampc_ffct_Sfct')
        delete([RESPATH '/' EXE{i},'.*']);
    else
        % sfunction must be in the same folder as the slx file or in the matlab path
        delete([EXE{i},'.*']);
    end
end

% list of all object files
OBJ = '';
for i = 1:length(ObjFiles)
    OBJ = [OBJ,' ',BINPATH,'/',ObjFiles{i}];
end
% compile probfct
try
    disp('Building probfct ...');
    eval(['mex -c ', probfct,' -I',CHEADERPATH,VERBOSE,DEBUG]);
    % remove the path at the beginning of the file name, if present
    probfctfilename = strsplit(probfct,'/');
    probfctfilename = strsplit(probfctfilename{end},'\');
    movefile([probfctfilename{end}(1:end-1),OBJEXT],[RESPATH '/probfct.',OBJEXT]);
catch e
    disp([probfct,': Invalid problem function.']);
    rethrow(e)
end
% compile executables
for i = 1:length(EXE)
    disp(['Building ',EXE{i},' ...']);
    try
        if ~strcmp(EXE{i},'grampc_run_Sfct') && ~strcmp(EXE{i},'grampc_ffct_Sfct')
            eval(['mex -output ',RESPATH '/' EXE{i},' ',BINPATH,'/',EXE{i},'.',OBJEXT,OBJ,' ' RESPATH '/probfct.',OBJEXT,' -I',CHEADERPATH,VERBOSE,DEBUG]);
        else
            % sfunction must be in the same folder als the slx file or in the matlab path
            eval(['mex -output ',EXE{i},' ',BINPATH,'/',EXE{i},'.',OBJEXT,OBJ,' ' RESPATH '/probfct.',OBJEXT,' -I',CHEADERPATH,VERBOSE,DEBUG]);
        end
    catch e
        disp([EXE{i},': Error during building process.'])
        rethrow(e)
    end
end
disp('Building process successively finished.');

% Updated matlab path to prevent the error "Undefined variable "CmexFiles ..."
rehash path;
end