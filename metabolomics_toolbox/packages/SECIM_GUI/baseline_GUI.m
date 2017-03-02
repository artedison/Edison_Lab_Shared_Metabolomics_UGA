function varargout = baseline_GUI(varargin)
% BASELINE_GUI MATLAB code for baseline_GUI.fig
%      BASELINE_GUI, by itself, creates a new BASELINE_GUI or raises the existing
%      singleton*.
%
%      H = BASELINE_GUI returns the handle to a new BASELINE_GUI or the handle to
%      the existing singleton*.
%
%      BASELINE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BASELINE_GUI.M with the given input arguments.
%
%      BASELINE_GUI('Property','Value',...) creates a new BASELINE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before baseline_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to baseline_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help baseline_GUI

% Last Modified by GUIDE v2.5 18-Aug-2014 09:16:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @baseline_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @baseline_GUI_OutputFcn, ...
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


% --- Executes just before baseline_GUI is made visible.
function baseline_GUI_OpeningFcn(hObject, eventdata, handles, parameters, ppm, X)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to baseline_GUI (see VARARGIN)
handles.XX=X;
handles.ppm=ppm;
handles.parameters= parameters;
% Choose default command line output for baseline_GUI
handles.output = hObject;
handles.ileft=1;
handles.iright=1;
default_A=5e-9*size(X,2)^4;
default_S=-quantile(X(:),.25);
set(handles.edit1,'string', num2str(default_A));
set(handles.edit2,'string', num2str(default_S));
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes baseline_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = baseline_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in bl1.
function bl1_Callback(hObject, eventdata, handles)
% hObject    handle to bl1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
A=str2double(get(handles.edit1,'String'));%handles.edit1;
S=str2double(get(handles.edit2,'String'));
size(handles.XX)
X=handles.XX(handles.ileft,:);
ppm=handles.ppm;
% A=5e9;
% S=100;
bd=Baseline(X',A,1.25,S);
plot(handles.axes1,ppm, X , ppm, bd*0,'k') 
hold(handles.axes1,'on')
plot(handles.axes1,ppm,bd,'r')
hold(handles.axes1,'off')
plot(handles.axes2, ppm, X-bd', ppm, bd*0,'k')
% x=handles.x;
% Y=handles.Y;
% plot(handles.axes1,x, Y) 
% hold(handles.axes1,'on')
% plot(handles.axes1,x,[1,1,1],'r')
% hold(handles.axes1,'off')
% plot(handles.axes2, x,Y-1,x, Y*0,'r')
linkaxes([handles.axes1,handles.axes2]);
set(handles.sample1,'string',num2str(handles.ileft));
handles.ileft=handles.ileft+1;
guidata(hObject, handles);

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sample1_Callback(hObject, eventdata, handles)
% hObject    handle to sample1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sample1 as text
%        str2double(get(hObject,'String')) returns contents of sample1 as a double


% --- Executes during object creation, after setting all properties.
function sample1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sample1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in reset1.
function reset1_Callback(hObject, eventdata, handles)
% hObject    handle to reset1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.sample1,'string', 1);
handles.ileft=1;
guidata(hObject, handles);

% --- Executes on button press in bl2.
function bl2_Callback(hObject, eventdata, handles)
% hObject    handle to bl2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
A=str2double(get(handles.edit3,'String'));
S=str2double(get(handles.edit4,'String'));
X=handles.XX(handles.iright,:);
ppm=handles.ppm;
bd=Baseline(X',A,1.25,S);
plot(handles.axes3,ppm, X,ppm, bd*0,'k') 
hold(handles.axes3,'on')
plot(handles.axes3,ppm,bd,'r')
hold(handles.axes3,'off')
plot(handles.axes4, ppm,X-bd',ppm, bd*0,'k')
linkaxes([handles.axes3,handles.axes4]);
set(handles.sample2,'string',num2str(handles.iright));
handles.iright=handles.iright+1;
guidata(hObject, handles);


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sample2_Callback(hObject, eventdata, handles)
% hObject    handle to sample2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sample2 as text
%        str2double(get(hObject,'String')) returns contents of sample2 as a double


% --- Executes during object creation, after setting all properties.
function sample2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sample2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in reset2.
function reset2_Callback(hObject, eventdata, handles)
% hObject    handle to reset2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.sample2,'string', 1);
handles.iright=1;
guidata(hObject, handles);


% --- Executes on button press in finalize.
function finalize_Callback(hObject, eventdata, handles)
% hObject    handle to finalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.radiobutton1, 'value')==1;
    handles.parameters.T.bl_noise = get(handles.edit1, 'string');
    handles.parameters.T.bl_fitting = get(handles.edit2, 'string');
else
    handles.parameters.T.bl_noise = get(handles.edit3, 'string');
    handles.parameters.T.bl_fitting = get(handles.edit4, 'string');
end
