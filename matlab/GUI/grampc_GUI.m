
function varargout = grampc_GUI(varargin)
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
% File: grampc_GUI.m
% Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
% Date: February 2014
% Version: v1.0
%
% MATLAB GUI interface of GRAMPC.
%

global CHEADERPATH;
global BINPATH;

% path of corresponding headers
CHEADERPATH = '../../include';
% path of bin folder
BINPATH = '../bin';

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @grampc_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @grampc_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before grampc_GUI is made visible.
function grampc_GUI_OpeningFcn(hObject, eventdata, handles, varargin)

global figRealTraj;
global figPredTraj;
global handleCheckboxRealTraj;
global handleCheckboxPredTraj;

global grampcGUI;

handles.output  = hObject;

% determine and save GUI as well as current path
handles.GUIpath  = mfilename('fullpath');
handles.GUIpath  = handles.GUIpath(1:end-11);
addpath(handles.GUIpath);

grampcGUI.START = 0;
grampcGUI.STOP  = 0;
grampcGUI.STEP  = 0;
grampcGUI.INIT  = 0;
grampcGUI.igrad = 0;

grampcGUI.vec.t       = [];
grampcGUI.vec.x       = [];
grampcGUI.vec.u       = [];
grampcGUI.vec.J       = [];
grampcGUI.vec.CPUtime = [];

handles.Tsim = 4.0;

handles.radiobutton.size = [22,21];
handles.textfield.size   = [95,21];

% MPC parameter default values of GUI
handles.GUI.param.Nx      = [];
handles.GUI.param.Nu      = [];
handles.GUI.param.Thor    = 0.5;
handles.GUI.param.dt      = 0.001;
handles.GUI.param.tk      = 0.0;
handles.GUI.param.Nhor    = 30;
handles.GUI.param.umax    = inf;
handles.GUI.param.umin    = -inf;
handles.GUI.param.xScale  = [];
handles.GUI.param.xOffset = [];
handles.GUI.param.uScale  = [];
handles.GUI.param.uOffset = [];
handles.GUI.param.xk      = [];
handles.GUI.param.u0      = [];
handles.GUI.param.xdes    = [];
handles.GUI.param.udes    = [];
handles.GUI.param.pSys    = [];
handles.GUI.param.pCost   = [];
handles.GUI.param.NpSys   = 0;
handles.GUI.param.NpCost  = 0;

% MPC options default values of GUI
handles.GUI.opt.MaxIter                  = 2;
handles.GUI.opt.ShiftControl             = 'on';
handles.GUI.opt.ScaleProblem             = 'off';
handles.GUI.opt.CostIntegrator           = 'trapezodial';
handles.GUI.opt.Integrator               = 'heun';
handles.GUI.opt.IntegratorRelTol         = 1e-6;
handles.GUI.opt.IntegratorAbsTol         = 1e-8;
handles.GUI.opt.LineSearchType           = 'adaptive';
handles.GUI.opt.LineSearchInit           = 5.5e-4;
handles.GUI.opt.LineSearchMin            = 1e-5;
handles.GUI.opt.LineSearchMax            = 0.75;
handles.GUI.opt.LineSearchIntervalFactor = 0.85;
handles.GUI.opt.LineSearchAdaptFactor    = 3.0/2.0;
handles.GUI.opt.LineSearchIntervalTol    = 1e-1;
handles.GUI.opt.JacobianX                = 'sysjacxadj';
handles.GUI.opt.JacobianU                = 'sysjacuadj';
handles.GUI.opt.IntegralCost             = 'on';
handles.GUI.opt.FinalCost                = 'on';

% MPC option and parameter values of grampcGUI
grampcGUI.param = handles.GUI.param;
grampcGUI.opt   = handles.GUI.opt;

% rws data of grampcGUI
grampcGUI.rws.t                  = 0;
grampcGUI.rws.x                  = 0;
grampcGUI.rws.xadj               = 0;
grampcGUI.rws.u                  = 0;
grampcGUI.rws.dHdu               = 0;
grampcGUI.rws.lsAdapt            = 0;
grampcGUI.rws.uls                = 0;
grampcGUI.rws.lsExplicit         = 0;
grampcGUI.rws.uprev              = 0;
grampcGUI.rws.dHduprev           = 0;
grampcGUI.rws.J                  = 0;
grampcGUI.rws.rwsScale           = 0;
grampcGUI.rws.rwsGradient        = 0;
grampcGUI.rws.rwsCostIntegration = 0;
grampcGUI.rws.rwsAdjIntegration  = 0;
grampcGUI.rws.rwsIntegration     = 0;

% update figure
popupmenu_MPCopt_Callback(handles.popupmenu_MPCopt, eventdata, handles);

% For close request during mpc run
set(hObject,'UserData',1);

% set position of GUI
set(0,'Units','pixels');
set(hObject,'Units','pixels');

ScreenSize = get(0,'ScreenSize');
P          = 30;
GUIsize   = get(hObject,'OuterPosition');
GUIwidth  = GUIsize(3);
GUIheight = GUIsize(4);
set(hObject,'OuterPosition',[P,ScreenSize(4)-GUIheight-P,GUIwidth,GUIheight]);

% For plotting results
figRealTraj                = createFigureRealTraj(handles);
figPredTraj                = createFigurePredTraj(handles);
handleCheckboxRealTraj     = handles.checkbox_RealTraj;
handleCheckboxPredTraj     = handles.checkbox_PredTraj;
handles.plot.idxStates     = '';
handles.plot.idxControls   = '';
handles.plot.idxLineSearch = '';
handles.plot.iSample       = str2num(get(handles.editText_PlotSamplerate,'String'));

% Update handles structure
guidata(hObject, handles);

% --- END FUNCTION


% --- Outputs from this function are returned to the command line.
function varargout = grampc_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- END FUNCTION


function plotResults(handles)

global figRealTraj;
global figPredTraj;
global grampcGUI;

if isempty(handles.plot.idxStates)
    x     = grampcGUI.vec.x;
    xpred = grampcGUI.rws.x;
else
    x     = grampcGUI.vec.x(str2num(handles.plot.idxStates),:); %#ok<*ST2NM>
    xpred = grampcGUI.rws.x(str2num(handles.plot.idxStates),:);
%     eval(['x = grampcGUI.vec.x(',handles.plot.idxStates,',:);']);
%     eval(['xpred = grampcGUI.rws.x(',handles.plot.idxStates,',:);']);
end

if isempty(handles.plot.idxControls)
    u     = grampcGUI.vec.u;
    upred = grampcGUI.rws.u;
else
    u     = grampcGUI.vec.u(str2num(handles.plot.idxControls),:);
    upred = grampcGUI.rws.u(str2num(handles.plot.idxControls),:);
%     eval(['u = grampcGUI.vec.u(',handles.plot.idxControls,',:);']);
%     eval(['upred = grampcGUI.rws.u(',handles.plot.idxControls,',:);']);
end

% plot real trajextories
if get(handles.checkbox_RealTraj,'Value') && ishandle(figRealTraj)
  % plot states
  subplot_States = subplot(2,2,1,'Parent',figRealTraj);
  plot(grampcGUI.vec.t,x,'Parent',subplot_States);
  xlim(subplot_States,[0,handles.Tsim]);
  title(subplot_States,'closed-loop states');
  xlabel(subplot_States,'time');
  % plot controls
  subplot_Controls = subplot(2,2,2,'Parent',figRealTraj);
  plot(grampcGUI.vec.t,u,'Parent',subplot_Controls);
  xlim(subplot_Controls,[0,handles.Tsim]);
  title(subplot_Controls,'closed-loop controls');
  xlabel(subplot_Controls,'time');
  % Plot cost
  subplot_Cost = subplot(2,2,3,'Parent',figRealTraj);
  plot(grampcGUI.vec.t,grampcGUI.vec.J,'Parent',subplot_Cost);
  xlim(subplot_Cost,[0,handles.Tsim]);
  title(subplot_Cost,'cost value');
  xlabel(subplot_Cost,'time');
  % Plot cpu time
  subplot_CPUtime = subplot(2,2,4,'Parent',figRealTraj);
  plot(grampcGUI.vec.t,grampcGUI.vec.CPUtime,'Parent',subplot_CPUtime);
  xlim(subplot_CPUtime,[0,handles.Tsim]);
  title(subplot_CPUtime,'computation time');
  xlabel(subplot_CPUtime,'time');
end

% plot predicted trajectories
if get(handles.checkbox_PredTraj,'Value') && ishandle(figPredTraj)
  if strcmp(handles.GUI.opt.LineSearchType,'adaptive')
    subplot_States   = subplot(1,3,1,'Parent',figPredTraj);
    subplot_Controls = subplot(1,3,2,'Parent',figPredTraj);
    subplot_ls       = subplot(1,3,3,'Parent',figPredTraj);
    % plot states
    plot(grampcGUI.rws.t,xpred,'Parent',subplot_States);
    xlim(subplot_States,[0,handles.GUI.param.Thor]);
    title(subplot_States,'predicted states');
    xlabel(subplot_States,'time');
    % plot controls
    plot(grampcGUI.rws.t,upred,'Parent',subplot_Controls);
    xlim(subplot_Controls,[0,handles.GUI.param.Thor]);
    title(subplot_Controls,'predicted controls');
    xlabel(subplot_Controls,'time');
    % plot line search process
    cla(subplot_ls);
    hold(subplot_ls,'all');
    if isempty(handles.plot.idxLineSearch)
        idx = 1:grampcGUI.opt.MaxIter;
    else
        idx = str2num(handles.plot.idxLineSearch);
    end    
    for i = 1:length(idx)
        [aa{i},jj{i}] = cost_approximation(grampcGUI.rws.lsAdapt(1+(idx(i)-1)*8:idx(i)*8));
        j{i}    = grampcGUI.rws.lsAdapt(5+(idx(i)-1)*8:idx(i)*8-1);
        a{i}    = grampcGUI.rws.lsAdapt(1+(idx(i)-1)*8:idx(i)*8-5);
        jopt{i} = grampcGUI.rws.lsAdapt(idx(i)*8);
        aopt{i} = grampcGUI.rws.lsAdapt(4+(idx(i)-1)*8);        
        plot(aa{i},jj{i},'b-','Parent',subplot_ls);
        plot(a{i},j{i},'Marker','.','MarkerSize',17,'MarkerEdgeColor','b','LineStyle','none','Parent',subplot_ls);
        plot(aopt{i},jopt{i},'Marker','.','MarkerSize',15,'MarkerEdgeColor','r','Parent',subplot_ls);
    end
%    xmin = a{1}(1)-abs(a{1}(1)-a{1}(2));
%    xmax = a{1}(end)+abs(a{1}(1)-a{1}(2));
%    ymin = jopt{1}-abs(j{1}(1)-j{1}(2));
%    ymax = max(j{1})+abs(j{1}(1)-j{1}(2));
%    for i = 2:length(idx)
%        xmin = min(xmin,a{i}(1)-abs(a{i}(1)-a{i}(2)));
%        xmax = max(xmax,a{i}(end)+abs(a{i}(1)-a{i}(2)));
%        ymin = min(ymin,jopt{i}-abs(j{i}(1)-j{i}(2)));
%        ymax = max(ymax,max(j{i})+abs(j{i}(1)-j{i}(2)));
%    end
    xmin = a{1}(1);
    xmax = a{1}(end);
%     ymin = jopt{1};
    ymin = min(j{1});
    ymax = max(j{1});
    for i = 2:length(idx)
        xmin = min(xmin,a{i}(1));
        xmax = max(xmax,a{i}(end));
%         ymin = min(ymin,jopt{i});
        ymin = min(ymin,min(jopt{i},min(j{i})));
        ymax = max(ymax,max(j{i}));
    end
    if ~(isnan(xmin) || isnan(xmax) || isnan(ymin) || isnan(ymax)) && ymax ~= ymin && xmax ~= xmin
      xlim(subplot_ls,[xmin-abs(xmax-xmin)*0.1 xmax+abs(xmax-xmin)*0.1]);
      ylim(subplot_ls,[ymin-abs(ymax-ymin)*0.1 ymax+abs(ymax-ymin)*0.1]);
    end
    title(subplot_ls,'line search iteration');
    xlabel(subplot_ls,'stepsize');
  else
    subplot_States   = subplot(1,2,1,'Parent',figPredTraj);
    subplot_Controls = subplot(1,2,2,'Parent',figPredTraj);
    % plot states
    plot(grampcGUI.rws.t,xpred,'Parent',subplot_States);
    xlim(subplot_States,[0,handles.GUI.param.Thor]);
    title(subplot_States,'predicted states');
    xlabel(subplot_States,'time');
    % plot controls
    plot(grampcGUI.rws.t,upred,'Parent',subplot_Controls);
    xlim(subplot_Controls,[0,handles.GUI.param.Thor]);
    title(subplot_Controls,'predicted controls');
    xlabel(subplot_Controls,'time');
  end
end

% --- END FUNCTION


% --- Executes during object creation, after setting all properties.
function popupmenu_MPCdata_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>

set(hObject,'BackgroundColor','white');

% --- END FUNCTION


function editText_MPCData_4_Callback(hObject, eventdata, handles)

handles.GUI.param.uOffset = str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- END FUNCTION


% --- Executes during object creation, after setting all properties.
function editText_MPCData_4_CreateFcn(hObject, eventdata, handles)

set(hObject,'BackgroundColor','white');

% --- END FUNCTION


function editText_MPCData_3_Callback(hObject, eventdata, handles)

handles.GUI.param.uScale = str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- END FUNCTION


% --- Executes during object creation, after setting all properties.
function editText_MPCData_3_CreateFcn(hObject, eventdata, handles)

set(hObject,'BackgroundColor','white');

% --- END FUNCTION


function editText_MPCData_2_Callback(hObject, eventdata, handles)

switch(get(handles.popupmenu_MPCdata,'Value'))
    case 4
        handles.GUI.param.u0 = str2num(get(hObject,'String'));
    case 5
        handles.GUI.param.udes = str2num(get(hObject,'String'));
    case 6
        handles.GUI.param.umax = str2num(get(hObject,'String'));
    case 7
        handles.GUI.param.xOffset = str2num(get(hObject,'String'));
end

guidata(hObject, handles);

% --- END FUNCTION


% --- Executes during object creation, after setting all properties.
function editText_MPCData_2_CreateFcn(hObject, eventdata, handles)

set(hObject,'BackgroundColor','white');

% --- END FUNCTION


function editText_MPCData_1_Callback(hObject, eventdata, handles)

switch(get(handles.popupmenu_MPCdata,'Value'))
    case 1
        handles.GUI.param.Thor = str2num(get(hObject,'String'));
    case 2
        handles.GUI.param.dt = str2num(get(hObject,'String'));
    case 3
        handles.GUI.param.Nhor = str2num(get(hObject,'String'));
    case 4
        handles.GUI.param.xk = str2num(get(hObject,'String'));
    case 5
        handles.GUI.param.xdes = str2num(get(hObject,'String'));
    case 6
        handles.GUI.param.umin = str2num(get(hObject,'String'));
    case 7
        handles.GUI.param.xScale = str2num(get(hObject,'String'));
    case 8
        handles.GUI.param.pSys  = str2num(get(hObject,'String'));
        handles.GUI.param.NpSys = length(handles.GUI.param.pSys);
    case 9
        handles.GUI.param.pCost  = str2num(get(hObject,'String'));
        handles.GUI.param.NpCost = length(handles.GUI.param.pCost);
    case 10
        handles.Tsim = str2num(get(hObject,'String'));
end

guidata(hObject, handles);

% --- END FUNCTION


% --- Executes during object creation, after setting all properties.
function editText_MPCData_1_CreateFcn(hObject, eventdata, handles)

set(hObject,'BackgroundColor','white');

% --- END FUNCTION


% --- Executes on selection change in popupmenu_MPCopt.
function popupmenu_MPCopt_Callback(hObject, eventdata, handles)

for i = 1:7
    if i == 1
        for j = 1:3
            eval(['set(handles.text_MPCopt_',num2str(i),num2str(j),',''Visible'',''off'')']);
            eval(['set(handles.editText_MPCopt_',num2str(i),num2str(j),',''Visible'',''off'')']);
        end
    else
        eval(['set(handles.text_MPCopt_',num2str(i),',''Visible'',''off'')']);
        eval(['set(handles.editText_MPCopt_',num2str(i),',''Visible'',''off'')']);
    end
end

tmp = get(handles.editText_MPCopt_11,'Position');
pos.editText_MPCopt_11 = tmp(1:2);
tmp = get(handles.editText_MPCopt_12,'Position');
pos.editText_MPCopt_12 = tmp(1:2);
tmp = get(handles.editText_MPCopt_13,'Position');
pos.editText_MPCopt_13 = tmp(1:2);

switch(get(hObject,'Value'))
    case 1  % Line search options
        set(handles.text_MPCopt_11,'String','adaptive');
        set(handles.text_MPCopt_12,'String','explicit1');
        set(handles.text_MPCopt_13,'String','explicit2');
        set(handles.editText_MPCopt_11,'Value',0);
        set(handles.editText_MPCopt_12,'Value',0);
        set(handles.editText_MPCopt_13,'Value',0);
        set(handles.editText_MPCopt_11,'Style','radiobutton');
        set(handles.editText_MPCopt_12,'Style','radiobutton');
        set(handles.editText_MPCopt_13,'Style','radiobutton');
        set(handles.editText_MPCopt_11,'Position',[pos.editText_MPCopt_11,handles.radiobutton.size]);
        set(handles.editText_MPCopt_12,'Position',[pos.editText_MPCopt_12,handles.radiobutton.size]);
        set(handles.editText_MPCopt_13,'Position',[pos.editText_MPCopt_13,handles.radiobutton.size]);
        set(handles.editText_MPCopt_11,'String','');
        set(handles.editText_MPCopt_12,'String','');
        set(handles.editText_MPCopt_13,'String','');
        set(handles.editText_MPCopt_11,'BackgroundColor',0.8*[1,1,1]);
        set(handles.editText_MPCopt_12,'BackgroundColor',0.8*[1,1,1]);
        set(handles.editText_MPCopt_13,'BackgroundColor',0.8*[1,1,1]);
        if strcmp(handles.GUI.opt.LineSearchType,'adaptive')
            set(handles.editText_MPCopt_11,'Value',1);
            set(handles.editText_MPCopt_12,'Value',0);
            set(handles.editText_MPCopt_13,'Value',0);
            %
            set(handles.editText_IdxLineSearch,'Style','edit');
            set(handles.editText_IdxLineSearch,'BackgroundColor','white');
            %
            set(handles.text_MPCopt_2,'String','init value');
            set(handles.editText_MPCopt_2,'Style','edit');
            set(handles.editText_MPCopt_2,'BackgroundColor',[1,1,1]);
            set(handles.editText_MPCopt_2,'String',mat2str(handles.GUI.opt.LineSearchInit));
            %
            set(handles.text_MPCopt_3,'String','min value');
            set(handles.editText_MPCopt_3,'String',mat2str(handles.GUI.opt.LineSearchMin));
            %
            set(handles.text_MPCopt_4,'String','max value');
            set(handles.editText_MPCopt_4,'String',mat2str(handles.GUI.opt.LineSearchMax));
            %
            set(handles.text_MPCopt_5,'String','interval factor');
            set(handles.editText_MPCopt_5,'String',mat2str(handles.GUI.opt.LineSearchIntervalFactor));
            %
            set(handles.text_MPCopt_6,'String','adapt. factor');
            set(handles.editText_MPCopt_6,'String',mat2str(handles.GUI.opt.LineSearchAdaptFactor));
            %
            set(handles.text_MPCopt_7,'String','interval tol.');
            set(handles.editText_MPCopt_7,'String',mat2str(handles.GUI.opt.LineSearchIntervalTol));
            %
            for i = 2:7
                eval(['set(handles.text_MPCopt_',num2str(i),',''Visible'',''on'')']);
                eval(['set(handles.editText_MPCopt_',num2str(i),',''Visible'',''on'')']);
            end
        elseif strcmp(handles.GUI.opt.LineSearchType,'explicit1')
            set(handles.editText_MPCopt_11,'Value',0);
            set(handles.editText_MPCopt_12,'Value',1);
            set(handles.editText_MPCopt_13,'Value',0);
            %
            set(handles.editText_IdxLineSearch,'Style','text');
            set(handles.editText_IdxLineSearch,'BackgroundColor',0.5*0.8*[1,1,1]);
            %
            set(handles.text_MPCopt_2,'String','init value');
            set(handles.editText_MPCopt_2,'String',mat2str(handles.GUI.opt.LineSearchInit));
            %
            set(handles.text_MPCopt_3,'String','min value');
            set(handles.editText_MPCopt_3,'String',mat2str(handles.GUI.opt.LineSearchMin));
            %
            set(handles.text_MPCopt_4,'String','max value');
            set(handles.editText_MPCopt_4,'String',mat2str(handles.GUI.opt.LineSearchMax));
            %
            for i = 2:4
                eval(['set(handles.text_MPCopt_',num2str(i),',''Visible'',''on'')']);
                eval(['set(handles.editText_MPCopt_',num2str(i),',''Visible'',''on'')']);
            end
        else
            set(handles.editText_MPCopt_11,'Value',0);
            set(handles.editText_MPCopt_12,'Value',0);
            set(handles.editText_MPCopt_13,'Value',1);
            %
            set(handles.editText_IdxLineSearch,'Style','text');
            set(handles.editText_IdxLineSearch,'BackgroundColor',0.5*0.8*[1,1,1]);
            %
            set(handles.text_MPCopt_2,'String','init value');
            set(handles.editText_MPCopt_2,'String',mat2str(handles.GUI.opt.LineSearchInit));
            %
            set(handles.text_MPCopt_3,'String','min value');
            set(handles.editText_MPCopt_3,'String',mat2str(handles.GUI.opt.LineSearchMin));
            %
            set(handles.text_MPCopt_4,'String','max value');
            set(handles.editText_MPCopt_4,'String',mat2str(handles.GUI.opt.LineSearchMax));
            %
            for i = 2:4
                eval(['set(handles.text_MPCopt_',num2str(i),',''Visible'',''on'')']);
                eval(['set(handles.editText_MPCopt_',num2str(i),',''Visible'',''on'')']);
            end
        end
        for j = 1:3
            eval(['set(handles.text_MPCopt_1',num2str(j),',''Visible'',''on'')']);
            eval(['set(handles.editText_MPCopt_1',num2str(j),',''Visible'',''on'')']);
        end
    case 2  % Maximum number of iterations
        set(handles.text_MPCopt_11,'String','value');
        set(handles.editText_MPCopt_11,'Style','edit');        
        set(handles.editText_MPCopt_11,'Position',[pos.editText_MPCopt_11,handles.textfield.size]);
        set(handles.editText_MPCopt_11,'BackgroundColor',[1,1,1]);
        set(handles.editText_MPCopt_11,'String',mat2str(handles.GUI.opt.MaxIter));
        set(handles.text_MPCopt_11,'Visible','on');
        set(handles.editText_MPCopt_11,'Visible','on');
    case 3  % Shift problem
        set(handles.text_MPCopt_11,'String','on');
        set(handles.text_MPCopt_12,'String','off');
        set(handles.editText_MPCopt_11,'String','');
        set(handles.editText_MPCopt_12,'String','');
        set(handles.editText_MPCopt_11,'Value',0);
        set(handles.editText_MPCopt_12,'Value',0);
        set(handles.editText_MPCopt_11,'Style','radiobutton');  
        set(handles.editText_MPCopt_12,'Style','radiobutton');
        set(handles.editText_MPCopt_11,'Position',[pos.editText_MPCopt_11,handles.radiobutton.size]);
        set(handles.editText_MPCopt_12,'Position',[pos.editText_MPCopt_12,handles.radiobutton.size]);
        set(handles.editText_MPCopt_11,'BackgroundColor',0.8*[1,1,1]);      
        set(handles.editText_MPCopt_12,'BackgroundColor',0.8*[1,1,1]);
        if strcmp(handles.GUI.opt.ShiftControl,'on')
            set(handles.editText_MPCopt_11,'Value',1);
            set(handles.editText_MPCopt_12,'Value',0);
        else
            set(handles.editText_MPCopt_11,'Value',0);
            set(handles.editText_MPCopt_12,'Value',1);
        end
        set(handles.text_MPCopt_11,'Visible','on');
        set(handles.text_MPCopt_12,'Visible','on');
        set(handles.editText_MPCopt_11,'Visible','on');
        set(handles.editText_MPCopt_12,'Visible','on');
    case 4  % Scale problem
        set(handles.text_MPCopt_11,'String','on');
        set(handles.text_MPCopt_12,'String','off');
        set(handles.editText_MPCopt_11,'String','');
        set(handles.editText_MPCopt_12,'String','');
        set(handles.editText_MPCopt_11,'Value',0);
        set(handles.editText_MPCopt_12,'Value',0);
        set(handles.editText_MPCopt_11,'Style','radiobutton');
        set(handles.editText_MPCopt_12,'Style','radiobutton');
        set(handles.editText_MPCopt_11,'Position',[pos.editText_MPCopt_11,handles.radiobutton.size]);
        set(handles.editText_MPCopt_12,'Position',[pos.editText_MPCopt_12,handles.radiobutton.size]);
        set(handles.editText_MPCopt_11,'BackgroundColor',0.8*[1,1,1]);
        set(handles.editText_MPCopt_12,'BackgroundColor',0.8*[1,1,1]);
        if strcmp(handles.GUI.opt.ScaleProblem,'on')
            set(handles.editText_MPCopt_11,'Value',1);
            set(handles.editText_MPCopt_12,'Value',0);
        else
            set(handles.editText_MPCopt_11,'Value',0);
            set(handles.editText_MPCopt_12,'Value',1);
        end
        set(handles.text_MPCopt_11,'Visible','on');
        set(handles.text_MPCopt_12,'Visible','on');
        set(handles.editText_MPCopt_11,'Visible','on');        
        set(handles.editText_MPCopt_12,'Visible','on');
    case 5  % Cost integration method
        set(handles.text_MPCopt_11,'String','trapezodial');
        set(handles.text_MPCopt_12,'String','simpson');
        set(handles.editText_MPCopt_11,'String','');
        set(handles.editText_MPCopt_12,'String','');
        set(handles.editText_MPCopt_11,'Value',0);
        set(handles.editText_MPCopt_12,'Value',0);
        set(handles.editText_MPCopt_11,'Style','radiobutton');
        set(handles.editText_MPCopt_12,'Style','radiobutton');
        set(handles.editText_MPCopt_11,'BackgroundColor',0.8*[1,1,1]);
        set(handles.editText_MPCopt_12,'BackgroundColor',0.8*[1,1,1]);
        set(handles.editText_MPCopt_11,'Position',[pos.editText_MPCopt_11,handles.radiobutton.size]);
        set(handles.editText_MPCopt_12,'Position',[pos.editText_MPCopt_12,handles.radiobutton.size]);
        if strcmp(handles.GUI.opt.CostIntegrator,'trapezodial')
            set(handles.editText_MPCopt_11,'Value',1);
            set(handles.editText_MPCopt_12,'Value',0);
        else
            set(handles.editText_MPCopt_11,'Value',0);
            set(handles.editText_MPCopt_12,'Value',1);
        end
        set(handles.text_MPCopt_11,'Visible','on');
        set(handles.text_MPCopt_12,'Visible','on');
        set(handles.editText_MPCopt_11,'Visible','on');
        set(handles.editText_MPCopt_12,'Visible','on');
    case 6  % Integrator
        set(handles.text_MPCopt_11,'String','type');
        set(handles.editText_MPCopt_11,'Value',0);
        set(handles.editText_MPCopt_11,'Style','popupmenu');
        set(handles.editText_MPCopt_11,'String',{'euler','modeuler','heun','ruku45'});
        set(handles.editText_MPCopt_11,'BackgroundColor',[1,1,1]);
        set(handles.editText_MPCopt_11,'Position',[pos.editText_MPCopt_11,handles.textfield.size]);
        switch (handles.GUI.opt.Integrator)
            case 'euler'
                set(handles.editText_MPCopt_11,'Value',1);
            case 'modeuler'
                set(handles.editText_MPCopt_11,'Value',2);
            case 'heun'
                set(handles.editText_MPCopt_11,'Value',3);
            case 'ruku45'
                set(handles.editText_MPCopt_11,'Value',4);
                set(handles.text_MPCopt_2,'String','Relative Tol.');
                set(handles.editText_MPCopt_2,'Style','edit');
                set(handles.editText_MPCopt_2,'String',mat2str(handles.GUI.opt.IntegratorRelTol));
                set(handles.editText_MPCopt_2,'BackgroundColor',[1,1,1]);
                set(handles.text_MPCopt_3,'String','Absolute Tol.');
                set(handles.editText_MPCopt_3,'String',mat2str(handles.GUI.opt.IntegratorAbsTol));
                for j = 2:3
                    eval(['set(handles.text_MPCopt_',num2str(j),',''Visible'',''on'')']);
                    eval(['set(handles.editText_MPCopt_',num2str(j),',''Visible'',''on'')']);
                end
            otherwise
                error('Error MPC options: Integratortype'); 
        end
        set(handles.text_MPCopt_11,'Visible','on');
        set(handles.editText_MPCopt_11,'Visible','on');
    case 7  % JacobianX
        set(handles.text_MPCopt_11,'String','sysjacxadj');
        set(handles.text_MPCopt_12,'String','sysjacx');
        set(handles.editText_MPCopt_11,'String','');
        set(handles.editText_MPCopt_12,'String','');
        set(handles.editText_MPCopt_11,'Value',0);
        set(handles.editText_MPCopt_12,'Value',0);
        set(handles.editText_MPCopt_11,'Style','radiobutton');
        set(handles.editText_MPCopt_12,'Style','radiobutton');
        set(handles.editText_MPCopt_11,'Position',[pos.editText_MPCopt_11,handles.radiobutton.size]);
        set(handles.editText_MPCopt_12,'Position',[pos.editText_MPCopt_12,handles.radiobutton.size]);
        set(handles.editText_MPCopt_11,'BackgroundColor',0.8*[1,1,1]);
        set(handles.editText_MPCopt_12,'BackgroundColor',0.8*[1,1,1]);
        if strcmp(handles.GUI.opt.JacobianX,'sysjacxadj')
            set(handles.editText_MPCopt_11,'Value',1);
            set(handles.editText_MPCopt_12,'Value',0);
        else
            set(handles.editText_MPCopt_11,'Value',0);
            set(handles.editText_MPCopt_12,'Value',1);
        end
        set(handles.text_MPCopt_11,'Visible','on');
        set(handles.text_MPCopt_12,'Visible','on');
        set(handles.editText_MPCopt_11,'Visible','on');
        set(handles.editText_MPCopt_12,'Visible','on');
    case 8  % JacobianU
        set(handles.text_MPCopt_11,'String','sysjacuadj');
        set(handles.text_MPCopt_12,'String','sysjacu');
        set(handles.editText_MPCopt_11,'String','');
        set(handles.editText_MPCopt_12,'String','');
        set(handles.editText_MPCopt_11,'Value',0);
        set(handles.editText_MPCopt_12,'Value',0);
        set(handles.editText_MPCopt_11,'Style','radiobutton');
        set(handles.editText_MPCopt_12,'Style','radiobutton');
        set(handles.editText_MPCopt_11,'Position',[pos.editText_MPCopt_11,handles.radiobutton.size]);
        set(handles.editText_MPCopt_12,'Position',[pos.editText_MPCopt_12,handles.radiobutton.size]);
        set(handles.editText_MPCopt_11,'BackgroundColor',0.8*[1,1,1]);
        set(handles.editText_MPCopt_12,'BackgroundColor',0.8*[1,1,1]);
        if strcmp(handles.GUI.opt.JacobianU,'sysjacuadj')
            set(handles.editText_MPCopt_11,'Value',1);
            set(handles.editText_MPCopt_12,'Value',0);
        else
            set(handles.editText_MPCopt_11,'Value',0);
            set(handles.editText_MPCopt_12,'Value',1);
        end
        set(handles.text_MPCopt_11,'Visible','on');
        set(handles.editText_MPCopt_11,'Visible','on');
        set(handles.text_MPCopt_12,'Visible','on');
        set(handles.editText_MPCopt_12,'Visible','on');
    case 9  % Integral cost
        set(handles.text_MPCopt_11,'String','on');
        set(handles.text_MPCopt_12,'String','off');
        set(handles.editText_MPCopt_11,'String','');
        set(handles.editText_MPCopt_12,'String','');
        set(handles.editText_MPCopt_11,'Value',0);
        set(handles.editText_MPCopt_12,'Value',0);
        set(handles.editText_MPCopt_11,'Style','radiobutton');
        set(handles.editText_MPCopt_12,'Style','radiobutton');
        set(handles.editText_MPCopt_11,'Position',[pos.editText_MPCopt_11,handles.radiobutton.size]);
        set(handles.editText_MPCopt_12,'Position',[pos.editText_MPCopt_12,handles.radiobutton.size]);
        set(handles.editText_MPCopt_11,'BackgroundColor',0.8*[1,1,1]);
        set(handles.editText_MPCopt_12,'BackgroundColor',0.8*[1,1,1]);
        if strcmp(handles.GUI.opt.IntegralCost,'on')
            set(handles.editText_MPCopt_11,'Value',1);
            set(handles.editText_MPCopt_12,'Value',0);
        else
            set(handles.editText_MPCopt_11,'Value',0);
            set(handles.editText_MPCopt_12,'Value',1);
        end
        set(handles.text_MPCopt_11,'Visible','on');
        set(handles.text_MPCopt_12,'Visible','on');
        set(handles.editText_MPCopt_11,'Visible','on');
        set(handles.editText_MPCopt_12,'Visible','on');
    case 10 % Final cost
        set(handles.text_MPCopt_11,'String','on');
        set(handles.text_MPCopt_12,'String','off');
        set(handles.editText_MPCopt_11,'String','');
        set(handles.editText_MPCopt_12,'String','');
        set(handles.editText_MPCopt_11,'Value',0);
        set(handles.editText_MPCopt_12,'Value',0);
        set(handles.editText_MPCopt_11,'Style','radiobutton');
        set(handles.editText_MPCopt_12,'Style','radiobutton');
        set(handles.editText_MPCopt_11,'Position',[pos.editText_MPCopt_11,handles.radiobutton.size]);
        set(handles.editText_MPCopt_12,'Position',[pos.editText_MPCopt_12,handles.radiobutton.size]);
        set(handles.editText_MPCopt_11,'BackgroundColor',0.8*[1,1,1]);
        set(handles.editText_MPCopt_12,'BackgroundColor',0.8*[1,1,1]);
        if strcmp(handles.GUI.opt.FinalCost,'on')
            set(handles.editText_MPCopt_11,'Value',1);
            set(handles.editText_MPCopt_12,'Value',0);
        else
            set(handles.editText_MPCopt_11,'Value',0);
            set(handles.editText_MPCopt_12,'Value',1);
        end
        set(handles.text_MPCopt_11,'Visible','on');
        set(handles.text_MPCopt_12,'Visible','on');
        set(handles.editText_MPCopt_11,'Visible','on');
        set(handles.editText_MPCopt_12,'Visible','on');
end

% --- END FUNCTION


% --- Executes during object creation, after setting all properties.
function popupmenu_MPCopt_CreateFcn(hObject, eventdata, handles)

set(hObject,'BackgroundColor','white');

% --- END FUNCTION


function editText_MPCopt_11_Callback(hObject, eventdata, handles)

switch (get(handles.popupmenu_MPCopt,'Value'))
    case 1
        set(handles.editText_MPCopt_11,'Value',1);
        set(handles.editText_MPCopt_12,'Value',0);
        set(handles.editText_MPCopt_13,'Value',0);
        handles.GUI.opt.LineSearchType = 'adaptive';
        popupmenu_MPCopt_Callback(handles.popupmenu_MPCopt, eventdata, handles);
        set(handles.editText_IdxLineSearch,'Style','edit');
        set(handles.editText_IdxLineSearch,'BackgroundColor','white');
    case 2
        handles.GUI.opt.MaxIter = str2num(get(hObject,'String'));
    case 3
        set(handles.editText_MPCopt_11,'Value',1);
        set(handles.editText_MPCopt_12,'Value',0);
        handles.GUI.opt.ShiftControl = 'on';
    case 4
        set(handles.editText_MPCopt_11,'Value',1);
        set(handles.editText_MPCopt_12,'Value',0);
        handles.GUI.opt.ScaleProblem = 'on';
    case 5
        set(handles.editText_MPCopt_11,'Value',1);
        set(handles.editText_MPCopt_12,'Value',0);
        handles.GUI.opt.CostIntegrator = 'trapezodial';
    case 6
        switch (get(hObject,'Value'))
            case 1
                handles.GUI.opt.Integrator = 'euler';
            case 2
                handles.GUI.opt.Integrator = 'modeuler';
            case 3
                handles.GUI.opt.Integrator = 'heun';
            case 4
                handles.GUI.opt.Integrator = 'ruku45';
        end
        popupmenu_MPCopt_Callback(handles.popupmenu_MPCopt, eventdata, handles);
    case 7
        set(handles.editText_MPCopt_11,'Value',1);
        set(handles.editText_MPCopt_12,'Value',0);
        handles.GUI.opt.JacobianX = 'sysjacxadj';
    case 8
        set(handles.editText_MPCopt_11,'Value',1);
        set(handles.editText_MPCopt_12,'Value',0);
        handles.GUI.opt.JacobianU = 'sysjacuadj';
    case 9
        set(handles.editText_MPCopt_11,'Value',1);
        set(handles.editText_MPCopt_12,'Value',0);
        handles.GUI.opt.IntegralCost = 'on';
    case 10
        set(handles.editText_MPCopt_11,'Value',1);
        set(handles.editText_MPCopt_12,'Value',0);
        handles.GUI.opt.FinalCost = 'on';
end
guidata(hObject, handles);

% --- END FUNCTION


% --- Executes during object creation, after setting all properties.
function editText_MPCopt_11_CreateFcn(hObject, eventdata, handles)

set(hObject,'BackgroundColor',0.8*[1,1,1]);

% --- END FUNCTION


% --- Executes on button press in editText_MPCopt_12.
function editText_MPCopt_12_Callback(hObject, eventdata, handles)

switch (get(handles.popupmenu_MPCopt,'Value'))
    case 1
        set(handles.editText_MPCopt_11,'Value',0);
        set(handles.editText_MPCopt_12,'Value',1);
        set(handles.editText_MPCopt_13,'Value',0);
        handles.GUI.opt.LineSearchType = 'explicit1';
        popupmenu_MPCopt_Callback(handles.popupmenu_MPCopt, eventdata, handles);
        set(handles.editText_IdxLineSearch,'Style','text');
        set(handles.editText_IdxLineSearch,'BackgroundColor',0.5*0.8*[1,1,1]);
    case 3
        set(handles.editText_MPCopt_11,'Value',0);
        set(handles.editText_MPCopt_12,'Value',1);
        handles.GUI.opt.ShiftControl = 'off';
    case 4
        set(handles.editText_MPCopt_11,'Value',0);
        set(handles.editText_MPCopt_12,'Value',1);
        handles.GUI.opt.ScaleProblem = 'off';
    case 5
        set(handles.editText_MPCopt_11,'Value',0);
        set(handles.editText_MPCopt_12,'Value',1);
        handles.GUI.opt.CostIntegrator = 'simpson';
    case 7
        set(handles.editText_MPCopt_11,'Value',0);
        set(handles.editText_MPCopt_12,'Value',1);
        handles.GUI.opt.JacobianX = 'sysjacx';
    case 8
        set(handles.editText_MPCopt_11,'Value',0);
        set(handles.editText_MPCopt_12,'Value',1);
        handles.GUI.opt.JacobianU = 'sysjacu';
    case 9
        set(handles.editText_MPCopt_11,'Value',0);
        set(handles.editText_MPCopt_12,'Value',1);
        handles.GUI.opt.IntegralCost = 'off';
    case 10
        set(handles.editText_MPCopt_11,'Value',0);
        set(handles.editText_MPCopt_12,'Value',1);
        handles.GUI.opt.FinalCost = 'off';
end
guidata(hObject, handles);

% --- END FUNCTION


% --- Executes during object creation, after setting all properties.
function editText_MPCopt_12_CreateFcn(hObject, eventdata, handles)

set(hObject,'BackgroundColor',0.8*[1,1,1]);

% --- END FUNCTION


% --- Executes on button press in editText_MPCopt_13.
function editText_MPCopt_13_Callback(hObject, eventdata, handles)

set(handles.editText_MPCopt_11,'Value',0);
set(handles.editText_MPCopt_12,'Value',0);
set(handles.editText_MPCopt_13,'Value',1);
handles.GUI.opt.LineSearchType = 'explicit2';
popupmenu_MPCopt_Callback(handles.popupmenu_MPCopt, eventdata, handles);
set(handles.editText_IdxLineSearch,'Style','text');
set(handles.editText_IdxLineSearch,'BackgroundColor',0.5*0.8*[1,1,1]);
        
guidata(hObject, handles);

% --- END FUNCTION


% --- Executes during object creation, after setting all properties.
function editText_MPCopt_13_CreateFcn(hObject, eventdata, handles)

set(hObject,'BackgroundColor',0.8*[1,1,1]);

% --- END FUNCTION


function editText_MPCopt_2_Callback(hObject, eventdata, handles)

switch (get(handles.popupmenu_MPCopt,'Value'))
    case 1
        handles.GUI.opt.LineSearchInit = str2num(get(hObject,'String'));
    case 6
        handles.GUI.opt.IntegratorRelTol = str2num(get(hObject,'String'));
end
guidata(hObject, handles);

% --- END FUNCTION


% --- Executes during object creation, after setting all properties.
function editText_MPCopt_2_CreateFcn(hObject, eventdata, handles)

set(hObject,'BackgroundColor','white');

% --- END FUNCTION


function editText_MPCopt_3_Callback(hObject, eventdata, handles)

switch (get(handles.popupmenu_MPCopt,'Value'))
    case 1
        handles.GUI.opt.LineSearchMin = str2num(get(hObject,'String'));
    case 6
        handles.GUI.opt.IntegratorAbsTol = str2num(get(hObject,'String'));
end
guidata(hObject, handles);

% --- END FUNCTION


% --- Executes during object creation, after setting all properties.
function editText_MPCopt_3_CreateFcn(hObject, eventdata, handles)

set(hObject,'BackgroundColor','white');

% --- END FUNCTION


function editText_MPCopt_4_Callback(hObject, eventdata, handles)

handles.GUI.opt.LineSearchMax = str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- END FUNCTION


% --- Executes during object creation, after setting all properties.
function editText_MPCopt_4_CreateFcn(hObject, eventdata, handles)

set(hObject,'BackgroundColor','white');

% --- END FUNCTION


function editText_MPCopt_5_Callback(hObject, eventdata, handles)

handles.GUI.opt.LineSearchIntervalFactor = str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- END FUNCTION


% --- Executes during object creation, after setting all properties.
function editText_MPCopt_5_CreateFcn(hObject, eventdata, handles)

set(hObject,'BackgroundColor','white');

% --- END FUNCTION


function editText_MPCopt_6_Callback(hObject, eventdata, handles)

handles.GUI.opt.LineSearchAdaptFactor = str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- END FUNCTION


% --- Executes during object creation, after setting all properties.
function editText_MPCopt_6_CreateFcn(hObject, eventdata, handles)

set(hObject,'BackgroundColor','white');

% --- END FUNCTION


function editText_MPCopt_7_Callback(hObject, eventdata, handles)

handles.GUI.opt.LineSearchIntervalTol = str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- END FUNCTION


% --- Executes during object creation, after setting all properties.
function editText_MPCopt_7_CreateFcn(hObject, eventdata, handles)

set(hObject,'BackgroundColor','white');

% --- END FUNCTION


function editText_IdxStates_Callback(hObject, eventdata, handles)

handles.plot.idxStates = get(hObject,'String');
guidata(hObject, handles);

% --- END FUNCTION


% --- Executes during object creation, after setting all properties.
function editText_IdxStates_CreateFcn(hObject, eventdata, handles)

set(hObject,'BackgroundColor','white');

% --- END FUNCTION



function editText_IdxControls_Callback(hObject, eventdata, handles)

handles.plot.idxControls = get(hObject,'String');
guidata(hObject, handles);

% --- END FUNCTION


% --- Executes during object creation, after setting all properties.
function editText_IdxControls_CreateFcn(hObject, eventdata, handles)

set(hObject,'BackgroundColor','white');

% --- END FUNCTION


function editText_IdxLineSearch_Callback(hObject, eventdata, handles)

handles.plot.idxLineSearch = get(hObject,'String');
idx = str2num(handles.plot.idxLineSearch);
for i = 1:length(idx)
  if idx(i) > handles.GUI.opt.MaxIter
    error('Selected indices does not correspond to number of iterations MaxIter.');
  end
end
guidata(hObject, handles);

% --- END FUNCTION


% --- Executes during object creation, after setting all properties.
function editText_IdxLineSearch_CreateFcn(hObject, eventdata, handles)

set(hObject,'BackgroundColor','white');

% --- END FUNCTION


% --- Executes on button press in checkbox_RealTraj.
function checkbox_RealTraj_Callback(hObject, eventdata, handles)

global figRealTraj;
global handleCheckboxRealTraj;

handleCheckboxRealTraj = hObject;

if ~get(hObject,'Value')
    set(hObject,'Value',0);
    if ishandle(figRealTraj)
        delete(figRealTraj);
        figRealTraj = [];
    end
else
    set(hObject,'Value',1);
    figRealTraj = createFigureRealTraj(handles);
end
guidata(hObject, handles);

% --- END FUNCTION


function figRealTraj = createFigureRealTraj(handles)

figWidth  = 560;
figHeight = 420; 

figRealTraj = figure('Visible','on','Units','pixels',...
                     'Name','Closed-loop trajectories',...
                     'CloseRequestFcn',@closeFigureRealTraj);

set(0,'Units','pixels');
set(handles.grampc_GUI,'Units','pixels');         

ScreenSize = get(0,'ScreenSize');
GUISize    = get(handles.grampc_GUI,'OuterPosition');
% FigSize    = get(figRealTraj,'OuterPosition');

xFig = GUISize(1)+GUISize(3);
yFig = ScreenSize(4)-(ScreenSize(4)-GUISize(2)-GUISize(4))-figHeight;
set(figRealTraj,'OuterPosition',[xFig,yFig,figWidth,figHeight]);

% --- END FUNCTION


function closeFigureRealTraj(varargin)

% global figRealTraj;
global handleCheckboxRealTraj;

if ishandle(handleCheckboxRealTraj)
  set(handleCheckboxRealTraj,'Value',0);
end
% if ishandle(figRealTraj)
%   delete(figRealTraj);
%   figRealTraj = [];
% end
closereq;

% --- END FUNCTION


% --- Executes on button press in checkbox_PredTraj.
function checkbox_PredTraj_Callback(hObject, eventdata, handles)

global figPredTraj;
global handleCheckboxPredTraj;

handleCheckboxPredTraj = hObject;

if ~get(hObject,'Value')
    set(hObject,'Value',0);
    if ishandle(figPredTraj)
        delete(figPredTraj);
        figPredTraj = [];
    end
else
    set(hObject,'Value',1);
    figPredTraj = createFigurePredTraj(handles);
end
guidata(hObject, handles);

% --- END FUNCTION


function figPredTraj = createFigurePredTraj(handles)

figWidth  = 560;
figHeight = 420; 

figPredTraj = figure('Visible','on','Units','pixels',...
                     'Name','Predicted trajectories',...
                     'CloseRequestFcn',@closeFigurePredTraj);

set(0,'Units','pixels');
set(handles.grampc_GUI,'Units','pixels'); 

ScreenSize = get(0,'ScreenSize');
GUISize    = get(handles.grampc_GUI,'OuterPosition');
% FigSize    = get(figPredTraj,'OuterPosition');

xFig = GUISize(1)+GUISize(3);
yFig = ScreenSize(4)-(ScreenSize(4)-GUISize(2)-GUISize(4))-2*figHeight;
set(figPredTraj,'OuterPosition',[xFig,yFig,figWidth,figHeight]);

% --- END FUNCTION


function closeFigurePredTraj(varargin)

% global figPredTraj;
global handleCheckboxPredTraj;

if ishandle(handleCheckboxPredTraj)
  set(handleCheckboxPredTraj,'Value',0);
end
% if ishandle(figPredTraj)
%     delete(figPredTraj);
%     figPredTraj = [];
% end
closereq;

% --- END FUNCTION


function editText_PlotSamplerate_Callback(hObject, eventdata, handles)

handles.plot.iSample = str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- END FUNCTION


% --- Executes during object creation, after setting all properties.
function editText_PlotSamplerate_CreateFcn(hObject, eventdata, handles)

set(hObject,'BackgroundColor','white');

% --- END FUNCTION


% --- Executes on selection change in popupmenu_MPCdata.
function popupmenu_MPCdata_Callback(hObject, eventdata, handles)

for i = 1:4
    eval(['set(handles.text_MPCData_',num2str(i),',''Visible'',''off'')']);
    eval(['set(handles.editText_MPCData_',num2str(i),',''Visible'',''off'')']);
end

switch(get(hObject,'Value'))
    case 1
        set(handles.text_MPCData_1,'String','Thor');
        set(handles.editText_MPCData_1,'String',num2str(handles.GUI.param.Thor));
        set(handles.text_MPCData_1,'Visible','on');
        set(handles.editText_MPCData_1,'Visible','on');
    case 2
        set(handles.text_MPCData_1,'String','dt');
        set(handles.editText_MPCData_1,'String',num2str(handles.GUI.param.dt));
        set(handles.text_MPCData_1,'Visible','on');
        set(handles.editText_MPCData_1,'Visible','on');
    case 3
        set(handles.text_MPCData_1,'String','Nhor');
        set(handles.editText_MPCData_1,'String',num2str(handles.GUI.param.Nhor));
        set(handles.text_MPCData_1,'Visible','on');
        set(handles.editText_MPCData_1,'Visible','on');
    case 4
        set(handles.text_MPCData_1,'String','x0');
        set(handles.text_MPCData_2,'String','u0');
        set(handles.editText_MPCData_1,'String',mat2str(handles.GUI.param.xk',3));
        set(handles.editText_MPCData_2,'String',mat2str(handles.GUI.param.u0',3));
        set(handles.text_MPCData_1,'Visible','on');
        set(handles.editText_MPCData_1,'Visible','on');
        set(handles.text_MPCData_2,'Visible','on');
        set(handles.editText_MPCData_2,'Visible','on');
    case 5
        set(handles.text_MPCData_1,'String','xdes');
        set(handles.text_MPCData_2,'String','udes');
        set(handles.editText_MPCData_1,'String',mat2str(handles.GUI.param.xdes',3));
        set(handles.editText_MPCData_2,'String',mat2str(handles.GUI.param.udes',3));
        set(handles.text_MPCData_1,'Visible','on');
        set(handles.editText_MPCData_1,'Visible','on');
        set(handles.text_MPCData_2,'Visible','on');
        set(handles.editText_MPCData_2,'Visible','on');
    case 6
        set(handles.text_MPCData_1,'String','umin');
        set(handles.text_MPCData_2,'String','umax');
        set(handles.editText_MPCData_1,'String',mat2str(handles.GUI.param.umin',3));
        set(handles.editText_MPCData_2,'String',mat2str(handles.GUI.param.umax',3));
        set(handles.text_MPCData_1,'Visible','on');
        set(handles.editText_MPCData_1,'Visible','on');
        set(handles.text_MPCData_2,'Visible','on');
        set(handles.editText_MPCData_2,'Visible','on');
    case 7
        set(handles.text_MPCData_1,'String','xScale');
        set(handles.text_MPCData_2,'String','xOffset');
        set(handles.text_MPCData_3,'String','uScale');
        set(handles.text_MPCData_4,'String','uOffset');
        set(handles.editText_MPCData_1,'String',mat2str(handles.GUI.param.xScale',3));
        set(handles.editText_MPCData_2,'String',mat2str(handles.GUI.param.xOffset',3));
        set(handles.editText_MPCData_3,'String',mat2str(handles.GUI.param.uScale',3));
        set(handles.editText_MPCData_4,'String',mat2str(handles.GUI.param.uOffset',3));
        set(handles.text_MPCData_1,'Visible','on');
        set(handles.editText_MPCData_1,'Visible','on');
        set(handles.text_MPCData_2,'Visible','on');
        set(handles.editText_MPCData_2,'Visible','on');
        set(handles.text_MPCData_3,'Visible','on');
        set(handles.editText_MPCData_3,'Visible','on');
        set(handles.text_MPCData_4,'Visible','on');
        set(handles.editText_MPCData_4,'Visible','on');
    case 8
        set(handles.text_MPCData_1,'String','pSys');
        set(handles.editText_MPCData_1,'String',mat2str(handles.GUI.param.pSys',3));
        set(handles.text_MPCData_1,'Visible','on');
        set(handles.editText_MPCData_1,'Visible','on');
    case 9
        set(handles.text_MPCData_1,'String','pCost');
        set(handles.editText_MPCData_1,'String',mat2str(handles.GUI.param.pCost',3));
        set(handles.text_MPCData_1,'Visible','on');
        set(handles.editText_MPCData_1,'Visible','on');
    case 10
        set(handles.text_MPCData_1,'String','Tsim');
        set(handles.editText_MPCData_1,'String',mat2str(handles.Tsim,3));
        set(handles.text_MPCData_1,'Visible','on');
        set(handles.editText_MPCData_1,'Visible','on');
end

% --- END FUNCTION


% --- Executes on button press in pushbutton_LoadData.
function pushbutton_LoadData_Callback(hObject, eventdata, handles)

handles.LOAD.grampcParamOpt.path     = [];
handles.LOAD.grampcParamOpt.filename = [];

[filename, path] = uigetfile('*.m', 'Select a MPC data function');
handles.LOAD.grampcParamOpt.path     = path;
handles.LOAD.grampcParamOpt.filename = filename;
guidata(hObject,handles); 

% --- END FUNCTION


% --- Executes on button press in pushbutton_LoadUserfct.
function pushbutton_LoadUserfct_Callback(hObject, eventdata, handles)

handles.LOAD.probfct.path     = [];
handles.LOAD.probfct.filename = [];

[filename, path] = uigetfile('*.c', 'Select a MPC probfct');
handles.LOAD.probfct.path     = path;
handles.LOAD.probfct.filename = filename;
guidata(hObject,handles);     

% --- END FUNCTION


% --- Executes on button press in pushbutton_Compile.
function pushbutton_Compile_Callback(hObject, eventdata, handles)

global CHEADERPATH;
global BINPATH;
global grampcGUI;

curPath = cd;
cd(handles.GUIpath);

if ~isfield(handles,'LOAD')
    commandwindow;
    error('Load a valid probfct and parameters/options.');
end

if ~isfield(handles.LOAD,'probfct')
  commandwindow;
  error('Load a valid probfct.');
end

if ~isfield(handles.LOAD,'grampcParamOpt')
  commandwindow;
  error('Load valid parameters/options.');
end

% delete previous mex files
file = {'grampc_init_Cmex','grampc_run_Cmex','grampc_setopt_Cmex','grampc_setparam_Cmex'};
for j = 1:length(file)
    warning('off');
        delete([file{j},'.*']);
    warning('on');
end

% compilation
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
       'grampc_run_Cmex',...
       'grampc_setopt_Cmex',...
       'grampc_setparam_Cmex'};

PROBFCT = [handles.LOAD.probfct.path,handles.LOAD.probfct.filename];
     
try
  disp('Building probfct ...');
  eval(['mex -c ',PROBFCT,' -I',CHEADERPATH]);
  if ~strcmp(handles.LOAD.probfct.filename,'probfct.c')
    movefile([handles.LOAD.probfct.filename(1:end-1),OBJEXT],['probfct.',OBJEXT]);
  end
catch %#ok<*CTCH>
  commandwindow;
  cd(curPath);
  error([PROBFCT,': Invalid probfct.']);
end

OBJ = '';
for i = 1:length(ObjFiles)
  OBJ = [OBJ,' ',BINPATH,'/',ObjFiles{i}]; %#ok<*AGROW>
end

for i = 1:length(EXE)
  disp(['Building ',EXE{i},' ...']);
  try
    eval(['mex -output ',EXE{i},' ',BINPATH,'/',EXE{i},'.',OBJEXT,OBJ,' probfct.',OBJEXT,' -I',CHEADERPATH]);
  catch 
    commandwindow;
    delete(['probfct.',OBJEXT]);
    cd(curPath);
    error([EXE{i},': Error during building process.']);
  end
end
disp('Building process successively finished.');

% delete userfunction
delete(['probfct.',OBJEXT]);

% update GUI data for current problem
grampc = grampc_init_Cmex(); %#ok<*ASGLU>;
handles.GUI.param = grampc.param;
grampcGUI.param   = grampc.param;
handles.GUI.opt   = grampc.opt;
grampcGUI.opt     = grampc.opt;
grampcGUI.rws     = grampc.rws;

try
    cd(handles.LOAD.grampcParamOpt.path);
    eval(['[grampc,handles.Tsim] = ',handles.LOAD.grampcParamOpt.filename(1:end-2),';']);
catch
    commandwindow;
    cd(curPath);
    error(['Provided data function ',handles.LOAD.grampcParamOpt.filename(1:end-2),' is not valid.']);
end

% Update parameters
if ~isempty(grampc.param)
    paramFields = fieldnames(grampc.param);
    for i = 1:length(paramFields)
        switch(paramFields{i})
          case 'Nx'
            handles.GUI.param.Nx = grampc.param.Nx;
          case 'Nu'
            handles.GUI.param.Nu = grampc.param.Nu;
          case 'xk'
            handles.GUI.param.xk = grampc.param.xk;
          case 'u0'
            handles.GUI.param.u0 = grampc.param.u0;
          case 'xdes'
            handles.GUI.param.xdes = grampc.param.xdes;
          case 'udes'
            handles.GUI.param.udes = grampc.param.udes;
          case 'Thor'
            handles.GUI.param.Thor = grampc.param.Thor;
          case 'dt'
            handles.GUI.param.dt = grampc.param.dt;
          case 'tk'
            handles.GUI.param.tk = grampc.param.tk;
          case 'Nhor'
            handles.GUI.param.Nhor = grampc.param.Nhor;
          case 'pCost'
            handles.GUI.param.pCost = grampc.param.pCost;
          case 'pSys'
            handles.GUI.param.pSys = grampc.param.pSys;
          case 'umax'
            handles.GUI.param.umax = grampc.param.umax;
          case 'umin'
            handles.GUI.param.umin = grampc.param.umin;
          case 'xScale'
            handles.GUI.param.xScale = grampc.param.xScale;
          case 'xOffset'
            handles.GUI.param.xOffset = grampc.param.xOffset;
          case 'uScale'
            handles.GUI.param.uScale = grampc.param.uScale;
          case 'uOffset'
            handles.GUI.param.uOffset = grampc.param.uOffset;
          case 'NpSys'
            handles.GUI.param.NpSys = grampc.param.NpSys;
          case 'NpCost'
            handles.GUI.param.NpCost = grampc.param.NpCost;
          otherwise
            commandwindow;
            cd(curPath);
            error([paramFields{i},': Invalid MPC parameter.']);
        end
    end
end

% Update MPC options
if ~isempty(grampc.opt)
    optFields = fieldnames(grampc.opt);
    for i = 1:length(optFields)
        switch(optFields{i})
            case 'MaxIter'
                handles.GUI.opt.MaxIter = grampc.opt.MaxIter;
            case 'ShiftControl'
                handles.GUI.opt.ShiftControl = grampc.opt.ShiftControl;
            case 'ScaleProblem'
                handles.GUI.opt.ScaleProblem = grampc.opt.ScaleProblem;
            case 'CostIntegrator'
                handles.GUI.opt.CostIntegrator = grampc.opt.CostIntegrator;
            case 'Integrator'
                handles.GUI.opt.Integrator = grampc.opt.Integrator;
            case 'IntegratorRelTol'
                handles.GUI.opt.IntegratorRelTol = grampc.opt.IntegratorRelTol;
            case 'IntegratorAbsTol'
                handles.GUI.opt.IntegratorAbsTol = grampc.opt.IntegratorAbsTol;
            case 'LineSearchType'
                handles.GUI.opt.LineSearchType = grampc.opt.LineSearchType;
            case 'LineSearchMax'
                handles.GUI.opt.LineSearchMax = grampc.opt.LineSearchMax;
            case 'LineSearchMin'
                handles.GUI.opt.LineSearchMin = grampc.opt.LineSearchMin;
            case 'LineSearchInit'
                handles.GUI.opt.LineSearchInit = grampc.opt.LineSearchInit;
            case 'LineSearchIntervalFactor'
                handles.GUI.opt.LineSearchIntervalFactor = grampc.opt.LineSearchIntervalFactor;
            case 'LineSearchAdaptFactor'
                handles.GUI.opt.LineSearchAdaptFactor = grampc.opt.LineSearchAdaptFactor;
            case 'LineSearchIntervalTol'
                handles.GUI.opt.LineSearchIntervalTol = grampc.opt.LineSearchIntervalTol;
            case 'JacobianX'
                handles.GUI.opt.JacobianX = grampc.opt.JacobianX;
            case 'JacobianU'
                handles.GUI.opt.JacobianU = grampc.opt.JacobianU;
            case 'IntegralCost'
                handles.GUI.opt.IntegralCost = grampc.opt.IntegralCost;
            case 'FinalCost'
                handles.GUI.opt.FinalCost = grampc.opt.FinalCost;
            otherwise
                commandwindow;
                cd(curPath);
                error([optFields{i},': Invalid option name.']);              
        end
    end
end

% update GUI
popupmenu_MPCdata_Callback(handles.popupmenu_MPCdata, eventdata, handles);
popupmenu_MPCopt_Callback(handles.popupmenu_MPCopt, eventdata, handles);

cd(curPath);

guidata(hObject,handles);

% --- END FUNCTION


% --- Executes on button press in pushbutton_Clean.
function pushbutton_Clean_Callback(hObject, eventdata, handles) %#ok<*INUSL,*DEFNU>

curPath = cd;
eval(['cd ',handles.GUIpath]);

% delete previous mex files
file    = {'grampc_init_Cmex','grampc_run_Cmex','grampc_setopt_Cmex','grampc_setparam_Cmex'};
for j = 1:length(file)
    warning('off');  %#ok<*WNOFF>
    delete([file{j},'.*']);
    warning('on'); %#ok<*WNON>
end

cd(curPath);

% --- END FUNCTION


% --- Executes on button press in pushbutton_Start.
function pushbutton_Start_Callback(hObject, eventdata, handles)

global grampcGUI;

if ~grampcGUI.START
    
  disp('MPC started ...');

  grampcGUI.START = 1;
  grampcGUI.STOP  = 0;
  grampcGUI.STEP  = 0;

  if ~grampcGUI.INIT
    grampc_init_opt_param(handles);
  end
  
  grampc_run(handles);
  
end

% --- END FUNCTION


% --- Executes on button press in pushbutton_Stop.
function pushbutton_Stop_Callback(hObject, eventdata, handles)

global grampcGUI;

disp('MPC stopped.');
grampcGUI.STOP  = 1;
grampcGUI.START = 0;
grampcGUI.STEP  = 0;
grampcGUI.INIT  = 0;

% --- END FUNCTION


% --- Executes on button press in pushbutton_Step.
function pushbutton_Step_Callback(hObject, eventdata, handles)

global grampcGUI;

if ~grampcGUI.STEP
    
  disp('MPC step ...');
  
  grampcGUI.STEP  = 1;
  grampcGUI.STOP  = 0;

  if ~grampcGUI.INIT
    grampc_init_opt_param(handles)
  end
  
  if ~grampcGUI.START 
    grampcGUI.START = 1;
    grampc_run(handles);
  end 

end

% --- END FUNCTION


% --- Executes on button press in pushbutton_Exit.
function pushbutton_Exit_Callback(hObject, eventdata, handles)

global grampcGUI;

grampcGUI.START = 0;
grampcGUI.STOP  = 1;
grampcGUI.STEP  = 0;

pushbutton_Clean_Callback(handles.pushbutton_Exit, eventdata, handles);

rmpath(handles.GUIpath);

closeFigureRealTraj();
closeFigurePredTraj();
% delete(handles.grampc_GUI);
closereq;

% --- END FUNCTION


% --- Starting MPC
function grampc_run(handles)

  global grampcGUI;
  grampc.param = grampcGUI.param;
  grampc.opt   = grampcGUI.opt;
  grampc.rws   = grampcGUI.rws;
  t = 0.0;

  while (grampcGUI.vec.t(grampcGUI.igrad) < handles.Tsim) && grampcGUI.START
    tic;
    [xnext,unext,Jnext] = grampc_run_Cmex(grampc);
    grampc_setparam_Cmex(grampc,'xk',xnext);
    t = t + grampc.param.dt;
    grampc_setparam_Cmex(grampc,'tk',t);
    grampcGUI.vec.CPUtime(grampcGUI.igrad) = toc;
    grampcGUI.vec.x(:,grampcGUI.igrad+1)   = xnext;
    grampcGUI.vec.u(:,grampcGUI.igrad+1)   = unext;
    grampcGUI.vec.J(grampcGUI.igrad+1)     = Jnext;
    grampcGUI.param = grampc.param;
    grampcGUI.rws   = grampc.rws;
    if mod(grampcGUI.igrad,handles.plot.iSample) == 0 || grampcGUI.STEP
      plotResults(handles);
      drawnow;
    end
    grampcGUI.igrad = grampcGUI.igrad + 1;
    if grampcGUI.STOP
      break;
    end
    if grampcGUI.STEP 
      grampcGUI.START = 0;
      grampcGUI.STEP  = 0;
    end
  end

  if grampcGUI.vec.t(grampcGUI.igrad) >= handles.Tsim
    disp('MPC calculation finished.');
    grampcGUI.START = 0;
    grampcGUI.STEP  = 0;
    grampcGUI.INIT  = 0;
    grampcGUI.STOP  = 0;
    % save results to workspace
    assignin('base','grampcGUI',grampcGUI);
  end

% --- END FUNCTION


% --- Setting grampcGUI's options and parameters
function grampc_init_opt_param(handles)

  global grampcGUI;

  grampcGUI.INIT = 1;
  
  % initialization
  grampc = grampc_init_Cmex();

  % set options
  grampc_setopt_Cmex(grampc,'MaxIter',handles.GUI.opt.MaxIter);
  grampc_setopt_Cmex(grampc,'ShiftControl',handles.GUI.opt.ShiftControl);
  grampc_setopt_Cmex(grampc,'ScaleProblem',handles.GUI.opt.ScaleProblem);
  grampc_setopt_Cmex(grampc,'CostIntegrator',handles.GUI.opt.CostIntegrator);
  grampc_setopt_Cmex(grampc,'Integrator',handles.GUI.opt.Integrator);
  grampc_setopt_Cmex(grampc,'LineSearchType',handles.GUI.opt.LineSearchType);
  grampc_setopt_Cmex(grampc,'IntegratorRelTol',handles.GUI.opt.IntegratorRelTol);
  grampc_setopt_Cmex(grampc,'IntegratorAbsTol',handles.GUI.opt.IntegratorAbsTol);
  grampc_setopt_Cmex(grampc,'LineSearchMax',handles.GUI.opt.LineSearchMax);
  grampc_setopt_Cmex(grampc,'LineSearchMin',handles.GUI.opt.LineSearchMin);
  grampc_setopt_Cmex(grampc,'LineSearchInit',handles.GUI.opt.LineSearchInit);
  grampc_setopt_Cmex(grampc,'LineSearchIntervalFactor',handles.GUI.opt.LineSearchIntervalFactor);
  grampc_setopt_Cmex(grampc,'LineSearchAdaptFactor',handles.GUI.opt.LineSearchAdaptFactor);
  grampc_setopt_Cmex(grampc,'LineSearchIntervalTol',handles.GUI.opt.LineSearchIntervalTol);
  grampc_setopt_Cmex(grampc,'JacobianX',handles.GUI.opt.JacobianX);
  grampc_setopt_Cmex(grampc,'JacobianU',handles.GUI.opt.JacobianU);
  grampc_setopt_Cmex(grampc,'IntegralCost',handles.GUI.opt.IntegralCost);
  grampc_setopt_Cmex(grampc,'FinalCost',handles.GUI.opt.FinalCost);

  % set MPC data
  grampc_setparam_Cmex(grampc,'xk',handles.GUI.param.xk);
  grampc_setparam_Cmex(grampc,'u0',handles.GUI.param.u0);
  grampc_setparam_Cmex(grampc,'xdes',handles.GUI.param.xdes);
  grampc_setparam_Cmex(grampc,'udes',handles.GUI.param.udes);
  grampc_setparam_Cmex(grampc,'Thor',handles.GUI.param.Thor);
  grampc_setparam_Cmex(grampc,'dt',handles.GUI.param.dt);
  grampc_setparam_Cmex(grampc,'tk',handles.GUI.param.tk);
  grampc_setparam_Cmex(grampc,'NpCost',handles.GUI.param.NpCost);
  if ~isempty(handles.GUI.param.pCost)
    grampc_setparam_Cmex(grampc,'pCost',handles.GUI.param.pCost);
  end
  grampc_setparam_Cmex(grampc,'NpSys',handles.GUI.param.NpSys);
  if ~isempty(handles.GUI.param.pSys)
    grampc_setparam_Cmex(grampc,'pSys',handles.GUI.param.pSys);
  end
  grampc_setparam_Cmex(grampc,'umax',handles.GUI.param.umax);
  grampc_setparam_Cmex(grampc,'umin',handles.GUI.param.umin);
  grampc_setparam_Cmex(grampc,'xScale',handles.GUI.param.xScale);
  grampc_setparam_Cmex(grampc,'xOffset',handles.GUI.param.xOffset);
  grampc_setparam_Cmex(grampc,'uScale',handles.GUI.param.uScale);
  grampc_setparam_Cmex(grampc,'uOffset',handles.GUI.param.uOffset);
  grampc_setparam_Cmex(grampc,'Nhor',handles.GUI.param.Nhor);

  grampcGUI.param = grampc.param;
  grampcGUI.opt   = grampc.opt;
  grampcGUI.rws   = grampc.rws;
  
  % data for plot
  grampcGUI.igrad = 1;

  grampcGUI.vec.t       = 0:grampcGUI.param.dt:handles.Tsim;
  grampcGUI.vec.x       = nan*zeros(grampcGUI.param.Nx,length(grampcGUI.vec.t));
  grampcGUI.vec.u       = nan*zeros(grampcGUI.param.Nu,length(grampcGUI.vec.t));
  grampcGUI.vec.J       = nan*zeros(1,length(grampcGUI.vec.t));
  grampcGUI.vec.CPUtime = nan*zeros(1,length(grampcGUI.vec.t));

  grampcGUI.vec.x(:,1) = grampcGUI.param.xk;
  
% --- END FUNCTION


function [alpha,J] = cost_approximation(lsAdapt)

  j = lsAdapt(5:7);
  a = lsAdapt(1:3);

  c0 = (j(3)*a(1)*(a(1)-a(2))*a(2)+a(3)*(j(1)*a(2)*(a(2)-a(3))+j(2)*a(1)*(-a(1)+a(3))))/((a(1)-a(2))*(a(1)-a(3))*(a(2)-a(3)));
  c1 = (j(3)*(-(a(1)*a(1))+(a(2)*a(2)))+j(2)*((a(1)*a(1))-(a(3)*a(3)))+j(1)*(-(a(2)*a(2))+(a(3)*a(3))))/((a(1)-a(2))*(a(1)-a(3))*(a(2)-a(3)));
  c2 = ((j(1)-j(3))/(a(1)-a(3))+(-j(2)+j(3))/(a(2)-a(3)))/(a(1)-a(2));

  alpha = linspace(a(1),a(3),30);
  J     = c2.*alpha.^2 + c1.*alpha + c0;
  
% --- END FUNCTION


