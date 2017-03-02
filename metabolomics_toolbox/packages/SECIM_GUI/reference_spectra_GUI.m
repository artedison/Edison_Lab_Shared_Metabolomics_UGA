function varargout = reference_spectra_GUI(varargin)
% REFERENCE_SPECTRA_GUI MATLAB code for reference_spectra_GUI.fig
%      REFERENCE_SPECTRA_GUI, by itself, creates a new REFERENCE_SPECTRA_GUI or raises the existing
%      singleton*.
%
%      H = REFERENCE_SPECTRA_GUI returns the handle to a new REFERENCE_SPECTRA_GUI or the handle to
%      the existing singleton*.
%
%      REFERENCE_SPECTRA_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REFERENCE_SPECTRA_GUI.M with the given input arguments.
%
%      REFERENCE_SPECTRA_GUI('Property','Value',...) creates a new REFERENCE_SPECTRA_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before reference_spectra_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to reference_spectra_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help reference_spectra_GUI

% Last Modified by GUIDE v2.5 30-Jul-2014 19:54:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @reference_spectra_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @reference_spectra_GUI_OutputFcn, ...
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


% --- Executes just before reference_spectra_GUI is made visible.
function reference_spectra_GUI_OpeningFcn(hObject, eventdata, handles, parameters, spectra)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to reference_spectra_GUI (see VARARGIN)
handles.parameters= parameters;
handles.spectra=spectra;
[X, ppm, Xtitles]=Setup1D(spectra);
handles.X=X;
handles.ppm=ppm;
plot(handles.axes1, ppm, X); 
%set(handles.axes1,'xlim', ppmrange)
set(handles.axes1,'xdir','rev')
% Choose default command line output for reference_spectra_GUI
handles.output = hObject;
handles.ileft=1;
handles.iright=1;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes reference_spectra_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = reference_spectra_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
reference_spectra(handles.ileft, handles.axes2, handles.spectra, [handles.edit1, handles.edit2]);
disp([handles.edit1, handles.edit2]);
% tmp= double(handles.i);
% disp(class(tmp))
% handles.i= tmp+1;
%[X,ppm, Xtitles]=Setup1D(spectra);
%plot(handles.axes2, ppm, X); 
%disp(handles.ileft);
set(handles.sampleleft,'string',num2str(handles.ileft));
handles.ileft= handles.ileft +1;
guidata(hObject, handles);

% --- Executes on button press in pushbutton2.
% function pushbutton2_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton2 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% reference_spectra(handles.iright, handles.axes3, handles.spectra, [handles.edit3,handles.edit4]);
% disp(handles.iright);
% set(handles.sampleright,'string',num2str(handles.iright));
% handles.iright= handles.iright +1;
% guidata(hObject, handles);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%if get(handles.radioleft, 'value')==1;
    handles.parameters.T.ref_thresh_1 = handles.edit1;
    handles.parameters.T.ref_thresh_2 = handles.edit2;
% else
%     handles.parameters.T.ref_thresh_1 = handles.edit3;
%     handles.parameters.T.ref_thresh_2 = handles.edit4;
% end


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
A=str2double(get(hObject,'String'));
handles.edit1=A;
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
guidata(hObject, handles);
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
A=str2double(get(hObject,'String'));
handles.edit2=A;
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
guidata(hObject, handles);
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



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
A=str2double(get(hObject,'String'));
handles.edit3=A;
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
guidata(hObject, handles);
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
A=str2double(get(hObject,'String'));
handles.edit4=A;
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
guidata(hObject, handles);
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



function sample_number_Callback(hObject, eventdata, handles)
% hObject    handle to sample_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sample_number as text
%        str2double(get(hObject,'String')) returns contents of sample_number as a double


% --- Executes during object creation, after setting all properties.
function sample_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sample_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in resetleft.
function resetleft_Callback(hObject, eventdata, handles)
% hObject    handle to resetleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ileft=1;
guidata(hObject, handles);


% --- Executes on button press in resetright.
function resetright_Callback(hObject, eventdata, handles)
% hObject    handle to resetright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.iright=1;
guidata(hObject, handles);



function sampleleft_Callback(hObject, eventdata, handles)
% hObject    handle to sampleleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'string',num2str(handles.ileft));
% Hints: get(hObject,'String') returns contents of sampleleft as text
%        str2double(get(hObject,'String')) returns contents of sampleleft as a double
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function sampleleft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sampleleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sampleright_Callback(hObject, eventdata, handles)
% hObject    handle to sampleright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'string',num2str(handles.iright));
% Hints: get(hObject,'String') returns contents of sampleright as text
%        str2double(get(hObject,'String')) returns contents of sampleright as a double
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function sampleright_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sampleright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
