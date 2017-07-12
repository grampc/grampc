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
% MATLAB/Simulink makefile to build the examples for running GRAMPC 
% in MATLAB/Simulink.
%
% >> MAKE PROBFCT
% Creates the executables according to the entered PROBFCT to run the
% grampc tool.
%
% >> MAKE clean
% Removes all previously generated executables.
%
% >> MAKE PROBFCT debug
% Creates the executables with additional information for use in 
% debugging. See also MEX command with flag -g.
%
% >> MAKE PROBFCT verbose
% Creates the executables with verbose option, i.e. additional 
% information are printed on the command window during the building process.
% See also MEX command with flag -v.
%
% >> MAKE PROBFCT verbose debug
% >> MAKE PROBFCT debug verbose
% Creates the executables with additional information for use in as well
% as with verbose option. See also description above.
%


% path of include folder
CHEADERPATH = '../../../include';
% path of bin folder
BINPATH = '../../bin';

% flags for compilation and cleaning
CLEAN   = 0;
DEBUG   = '';
VERBOSE = '';
PROBFCT = '';

if isempty(varargin)
  error('No probfct entered.');
else
  for i = 1:nargin
    switch varargin{i}
      case 'clean'
        CLEAN = 1;
      case 'debug'
        DEBUG = ' -g';
      case 'verbose'
        VERBOSE = ' -v';
      otherwise
        if strcmp('c',varargin{i}(end))
          PROBFCT = varargin{i};
        else
          error('Undefined input argument or invalid problem function.');
        end
    end
  end
end

if ~CLEAN && isempty(PROBFCT)
  error('First input argument must either be an probfct or the clean flag.');
end

if strcmp(computer,'PCWIN') || strcmp(computer,'PCWIN64')
  OBJEXT = 'obj';
else
  OBJEXT = 'o';
end

ObjFiles = {['grampc_run.',OBJEXT],...
            ['grampc_init.',OBJEXT],...
            ['grampc_mess.',OBJEXT],...
            ['grampc_setopt.',OBJEXT],...
            ['grampc_setparam.',OBJEXT],...
            ['euler1.',OBJEXT],...
            ['eulermod2.',OBJEXT],...
            ['heun2.',OBJEXT],...
            ['ruku45.',OBJEXT]};
          
EXE = {'grampc_init_Cmex',...
       'grampc_setopt_Cmex',...
       'grampc_setparam_Cmex',...
       'grampc_run_Cmex',...
       'grampc_run_Sfct'};

if CLEAN
  warning off; %#ok<*WNOFF>
  try %#ok<*TRYNC>
    delete(['probfct.',OBJEXT]);
  end
  for i = 1:length(EXE)
    try
      delete([EXE{i},'.*']);
    end
  end
  warning on;
end
     
if ~isempty(PROBFCT)
  OBJ = '';
  for i = 1:length(ObjFiles)
    OBJ = [OBJ,' ',BINPATH,'/',ObjFiles{i}]; %#ok<*AGROW>
  end
  try
    disp('Building probfct ...');
    eval(['mex -c ',PROBFCT,' -I',CHEADERPATH,VERBOSE,DEBUG]);
    if ~strcmp('probfct.c',PROBFCT)
      movefile([PROBFCT(1:end-1),OBJEXT],['probfct.',OBJEXT]);
    end
  catch e
    disp([PROBFCT,': Invalid problem function.']);
    rethrow(e)
  end
  for i = 1:length(EXE)
    disp(['Building ',EXE{i},' ...']);
    try
      eval(['mex -output ',EXE{i},' ',BINPATH,'/',EXE{i},'.',OBJEXT,OBJ,' probfct.',OBJEXT,' -I',CHEADERPATH,VERBOSE,DEBUG]);
    catch e %#ok<*CTCH>
      disp([EXE{i},': Error during building process.'])
      rethrow(e)
    end
  end
  disp('Building process successively finished.');
end

% ******* END OF FUNCTION *******
end




