function varargout = MI_value(varargin)
%MI_VALUE M-file for MI_value.fig
%      MI_VALUE, by itself, creates a new MI_VALUE or raises the existing
%      singleton*.
%
%      H = MI_VALUE returns the handle to a new MI_VALUE or the handle to
%      the existing singleton*.
%
%      MI_VALUE('Property','Value',...) creates a new MI_VALUE using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to MI_value_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      MI_VALUE('CALLBACK') and MI_VALUE('CALLBACK',hObject,...) call the
%      local function named CALLBACK in MI_VALUE.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MI_value

% Last Modified by GUIDE v2.5 02-Dec-2013 14:25:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MI_value_OpeningFcn, ...
                   'gui_OutputFcn',  @MI_value_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before MI_value is made visible.
function MI_value_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
% first=set(handles.edit_MI,'String',0.1);
% d=str2double(first);


% e=set(handles.edit_MI, 'String', '0.1');
% d=str2double(e);
d=0;
set(handles.edit_MI,'String', '0');
axes(handles.axes_main);
solv_ode(d);

% Choose default command line output for MI_value
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MI_value wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MI_value_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_MI_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MI as text
%        str2double(get(hObject,'String')) returns contents of edit_MI as a double


m=str2double(get(hObject,'String'));
% if isempty(m)
%       d = 0.1;     % when not input a number ,the contents is 5;
%    end 
axes(handles.axes_main);
solv_ode(m);

% --- Executes during object creation, after setting all properties.
function edit_MI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in video.
function video_Callback(hObject, eventdata, handles)
% hObject    handle to video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
m=str2double(get(handles.edit_MI,'String'));
if m<2
for i=m:0.05:5
    figure(1);
    h2=solv_ode(i);
 n=num2str(i);
 nn= strcat('MI=', n);
legend(h2,nn);
    drawnow 

guidata(hObject, handles);
end
end
close 

return

