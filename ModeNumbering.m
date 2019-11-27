function varargout = ModeNumbering(varargin)
% MODENUMBERING MATLAB code for ModeNumbering.fig
%      MODENUMBERING, by itself, creates a new MODENUMBERING or raises the existing
%      singleton*.
%
%      H = MODENUMBERING returns the handle to a new MODENUMBERING or the handle to
%      the existing singleton*.
%
%      MODENUMBERING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODENUMBERING.M with the given input arguments.
%
%      MODENUMBERING('Property','Value',...) creates a new MODENUMBERING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ModeNumbering_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ModeNumbering_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ModeNumbering

% Last Modified by GUIDE v2.5 05-Oct-2017 17:33:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ModeNumbering_OpeningFcn, ...
                   'gui_OutputFcn',  @ModeNumbering_OutputFcn, ...
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


% --- Executes just before ModeNumbering is made visible.
function ModeNumbering_OpeningFcn(hObject, eventData, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ModeNumbering (see VARARGIN)

set(handles.Pump_Cos_uitable, 'Data', varargin{1}) ;
set(handles.Pump_Sin_uitable, 'Data', varargin{2}) ;
set(handles.PumpTitle_et, 'String', 'Pump Modes') ;

set(handles.Signal_Cos_uitable, 'Data', varargin{3}) ;
set(handles.Signal_Sin_uitable, 'Data', varargin{4}) ;
set(handles.SignalTitle_et, 'String', 'Signal Modes') ;

% Choose default command line output for ModeNumbering
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ModeNumbering wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ModeNumbering_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function PumpTitle_et_Callback(hObject, eventdata, handles)
% hObject    handle to PumpTitle_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PumpTitle_et as text
%        str2double(get(hObject,'String')) returns contents of PumpTitle_et as a double


% --- Executes during object creation, after setting all properties.
function PumpTitle_et_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PumpTitle_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SignalTitle_et_Callback(hObject, eventdata, handles)
% hObject    handle to SignalTitle_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SignalTitle_et as text
%        str2double(get(hObject,'String')) returns contents of SignalTitle_et as a double


% --- Executes during object creation, after setting all properties.
function SignalTitle_et_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SignalTitle_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
