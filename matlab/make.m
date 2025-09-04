function [] = make(varargin)
% MATLAB/Simulink makefile to build the necessary executables for
% running GRAMPC in MATLAB/Simulink.
%
% >> MAKE
% >> MAKE all
% Creates necessary object files from the sources to use the grampc tool.
%
% >> MAKE clean
% Removes all previously generated object files.
%
% >> MAKE all debug
% >> MAKE debug
% Creates the necessary object files with additional information for use in
% debugging. See also MEX command with flag -g.
%
% >> MAKE all verbose
% >> MAKE verbose
% Creates the necessary object files with verbose option, i.e. additional
% information are printed on the command window during the building process.
% See also MEX command with flag -v.
%
% >> MAKE all timer N
% >> MAKE timer N
% Creates the necessary object files with the define N_TIMER N. This causes
% at every call of the grampc_run_Cmex function a time measurement of N
% identical MPC iterations. See also grampc_run_Cmex.c.
%
% >> MAKE all verbose debug
% >> MAKE all debug verbose
% >> MAKE verbose debug
% >> MAKE debug verbose
% >> MAKE all verbose debug timer N
% >> MAKE all debug verbose timer N
% >> MAKE verbose debug timer N
% >> MAKE debug verbose timer N
% Creates the object files with additional information for use with debug as
% well as with verbose options. See also description above.
%
%
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
%

% paths
CSRCPATH = '../../src';
% path of corresponding headers
CHEADERPATH = '../../include';
% path of mex headers
MEXHEADERPATH = '../include';
% path of mex sources
MEXSRCPATH = '../src';
% path of bin folder
BINPATH = './bin';

% current path
destPath = cd;

% flags for compilation and cleaning
CLEAN     = 0;
DEBUG     = '';
VERBOSE   = '';
TIMER     = '';

for i = 1:nargin
    switch varargin{i}
        case 'clean'
            CLEAN = 1;
        case 'debug'
            DEBUG = ' -g';
        case 'verbose'
            VERBOSE = ' -v';
        case 'timer'
            TIMER = [' -DN_TIMER=' num2str(varargin{i+1},'%d') ' '];
        otherwise
            if(~strcmp(varargin{i-1},'timer'))
                error('Undefined input argument.');
            end
    end
end

CFiles    = {'grampc_run.c',...
    'grampc_alloc.c',...
    'grampc_erk.c',...
    'grampc_fixedsize.c',...
    'grampc_init.c',...
    'grampc_mess.c',...
    'grampc_configreader.c',...
    'grampc_setopt.c',...
    'grampc_setparam.c',...
    'grampc_util.c',...
    'discrete.c',...
    'ruku45.c',...
    'rodas.c',...
    'simpson.c',...
    'timing.c',...
    'finite_diff.c',...
    'trapezoidal.c'};
CmexFiles = {'grampc_init_Cmex.c',...
    'grampc_get_config_from_file_Cmex.c',...
    'grampc_setopt_Cmex.c',...
    'grampc_setparam_Cmex.c',...
    'grampc_run_Cmex.c',...
    'grampc_run_Sfct.c',...
    'grampc_ghfct_Cmex.c',...
    'grampc_gThTfct_Cmex.c',...
    'grampc_Vfct_Cmex.c',...
    'grampc_lfct_Cmex.c',...
    'grampc_ffct_Cmex.c',...
    'grampc_estim_penmin_Cmex.c',...
    'grampc_check_gradients_Cmex.c',...
    'grampc_printopt_Cmex.c',...
    'grampc_printparam_Cmex.c',...
    'grampc_printstatus_Cmex.c'};

if ~exist(BINPATH,'dir')
    mkdir(BINPATH);
end

if CLEAN
    cd(BINPATH);
    for i = 1:length(CFiles)
        try %#ok<*TRYNC>
            delete([CFiles{i}(1:end-1),'*']);
        end
    end
    for i = 1:length(CmexFiles)
        try
            delete([CmexFiles{i}(1:end-1),'*']);
        end
    end
    cd(destPath);
end


cd(BINPATH);
% copying of corresponding source files into bin folder
% and subsequent compilation
for i=1:length(CFiles)
    copyfile([CSRCPATH,'/',CFiles{i}]);
end
for i=1:length(CmexFiles)
    copyfile([MEXSRCPATH,'/',CmexFiles{i}]);
end
% compilation process
for i = 1:length(CFiles)
    disp(['Building ',CFiles{i}(1:end-2),' ...']);
    try
        eval(['mex -c ',CFiles{i},' -DMEXCOMPILE -I',CHEADERPATH,VERBOSE,DEBUG,TIMER]);
    catch e %#ok<*CTCH>
        delete('*');
        cd(destPath);
        disp('Error during building process.');
        rethrow(e)
    end
end
for i = 1:length(CmexFiles)
    disp(['Building ',CmexFiles{i}(1:end-2),' ...']);
    try
        eval(['mex -c ',CmexFiles{i},' -I',CHEADERPATH,' -I',MEXHEADERPATH,VERBOSE,DEBUG,TIMER]);
    catch e %#ok<*CTCH>
        delete('*');
        cd(destPath);
        disp('Error during building process.');
        rethrow(e)
    end
end
disp('Building process successfully finished.');
for i=1:length(CFiles)
    delete(CFiles{i});
end
for i=1:length(CmexFiles)
    delete(CmexFiles{i});
end
% go back to destination
cd(destPath);

end
