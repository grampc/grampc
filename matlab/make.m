function [] = make(varargin)
%
% This file is part of GRAMPC.
%
% GRAMPC - a gradient-based MPC software for real-time applications
%
% Copyright (C) 2014 by Bartosz Kaepernick, Knut Graichen, Tilman Utz
% Developed at the Institute of Measurement, Control, and
% Microtechnology, University of Ulm. All rights reserved.
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
%
%
% File: make.m
% Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
% Date: February 2014
% Version: v1.0
%
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
% >> MAKE all verbose debug
% >> MAKE all debug verbose
% >> MAKE verbose debug
% >> MAKE debug verbose
% Creates the object files with additional information for use with debug as 
% well as with verbose options. See also description above.
%

% paths
CSRCPATH = '../../src';
% path of corresponding headers
CHEADERPATH = '../../include';
% path of mex sources
MEXSRCPATH = '../src';
% path of bin folder
BINPATH = './bin';

% current path
destPath = cd;

% flags for compilation and cleaning
CLEAN     = 0;
ALL       = 0;
DEBUG     = '';
VERBOSE   = '';

if isempty(varargin)
  ALL = 1;
else
  for i = 1:nargin
    switch varargin{i}
      case 'all'
        ALL = 1;
      case 'clean'
        CLEAN = 1;
      case 'debug'
        DEBUG = ' -g';
        ALL   = 1;
      case 'verbose'
        VERBOSE = ' -v';
        ALL     = 1;
      otherwise
        error('Undefined input argument.');
    end
  end
end

CFiles    = {'grampc_run.c',...
             'grampc_init.c',...
             'grampc_mess.c',...
             'grampc_setopt.c',...
             'grampc_setparam.c',...
             'euler1.c',...
             'eulermod2.c',...
             'heun2.c',...
             'ruku45.c'};
CmexFiles = {'grampc_init_Cmex.c',...
             'grampc_setopt_Cmex.c',...
             'grampc_setparam_Cmex.c',...
             'grampc_run_Cmex.c',...
             'grampc_run_Sfct.c'};

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

if ALL 
  if ~exist(BINPATH,'var')
      mkdir(BINPATH);
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
      eval(['mex -c ',CFiles{i},' -DMEXCOMPILE -I',CHEADERPATH,VERBOSE,DEBUG]);
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
      eval(['mex -c ',CmexFiles{i},' -I',CHEADERPATH,VERBOSE,DEBUG]);
    catch e %#ok<*CTCH>
      delete('*');
      cd(destPath);
      disp('Error during building process.');
      rethrow(e)
    end
  end
  disp('Building process successively finished.');
  for i=1:length(CFiles)
    delete(CFiles{i});
  end
  for i=1:length(CmexFiles)
    delete(CmexFiles{i});
  end
  % go back to destination
  cd(destPath);
end


% ******* END OF FUNCTION *******
end




