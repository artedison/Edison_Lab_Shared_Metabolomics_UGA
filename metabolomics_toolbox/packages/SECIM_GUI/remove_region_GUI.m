function varargout = remove_region_GUI(varargin)
%REMOVE_REGION_GUI M-file for remove_region_GUI.fig
%      REMOVE_REGION_GUI, by itself, creates a new REMOVE_REGION_GUI or raises the existing
%      singleton*.
%
%      H = REMOVE_REGION_GUI returns the handle to a new REMOVE_REGION_GUI or the handle to
%      the existing singleton*.
%
%      REMOVE_REGION_GUI('Property','Value',...) creates a new REMOVE_REGION_GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to remove_region_GUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      REMOVE_REGION_GUI('CALLBACK') and REMOVE_REGION_GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in REMOVE_REGION_GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help remove_region_GUI

% Last Modified by GUIDE v2.5 29-Jul-2014 14:25:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @remove_region_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @remove_region_GUI_OutputFcn, ...
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


% --- Executes just before remove_region_GUI is made visible.
function remove_region_GUI_OpeningFcn(hObject, eventdata, handles, parameters, X, ppm)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
handles.X=X;
handles.ppm=ppm;
handles.parameters= parameters;
plot(handles.axes1,ppm, X);
% Choose default command line output for remove_region_GUI
handles.output = hObject;
handles.region='water_region';
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes remove_region_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = remove_region_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.leftppm1=str2double(get(hObject,'String'));
%handles.edit1=A;
%set(hObject,'string',get(hObject,'String'));
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
guidata(hObject, handles);

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
% A=str2double(get(hObject,'String'));
% handles.edit2=A;
% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
guidata(hObject, handles);

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
% A=str2double(get(hObject,'String'));
% handles.edit3=A;
% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
guidata(hObject, handles);

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
% A=str2double(get(hObject,'String'));
% handles.edit4=A;
% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
guidata(hObject, handles);

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


% --- Executes on button press in finalize.
function finalize_Callback(hObject, eventdata, handles)
% hObject    handle to finalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.radioleft, 'value')==1;
    switch handles.region
        case 'water_region'
            handles.parameters.T.water_region_1 = str2num(get(handles.edit1,'string'));
            handles.parameters.T.water_region_2 = str2num(get(handles.edit2,'string'));
        case 'other_region'
            handles.parameters.T.region_removed_1 =str2num(get(handles.edit1,'string'));
            handles.parameters.T.region_removed_2 = str2num(get(handles.edit2,'string'));
        case 'ends_region'
            handles.parameters.T.ends_1 = str2num(get(handles.edit1,'string'));
            handles.parameters.T.ends_2 = str2num(get(handles.edit2,'string'));    
    end
else
    switch handles.region
        case 'water_region'
            handles.parameters.T.water_region_1 = str2num(get(handles.edit3,'string'));
            handles.parameters.T.water_region_2 = str2num(get(handles.edit4,'string'));
        case 'other_region'
            handles.parameters.T.region_removed_1 = str2num(get(handles.edit3,'string'));
            handles.parameters.T.region_removed_2 = str2num(get(handles.edit4,'string'));
        case 'ends_region'
            handles.parameters.T.ends_1 = str2num(get(handles.edit3,'string'));
            handles.parameters.T.ends_2 = str2num(get(handles.edit4,'string'));    
    end
end


% --- Executes on button press in remove_left.
function remove_left_Callback(hObject, eventdata, handles)
% hObject    handle to remove_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
shift1=str2double(get(handles.edit1,'String'));%handles.edit1;
shift2=str2double(get(handles.edit2,'String'));%handles.edit2;
diff=shift2-shift1;
X=handles.X;
ppm=handles.ppm;
disp(handles.region)
if (isequal(handles.region,'water_region') | isequal(handles.region,'other_region') )%(handles.region=='water_region' | handles.region=='other_region')
    XR=remove_region(X, ppm, shift1, shift2);
    plot(handles.axes2, ppm, XR);
    set(handles.axes2,'xlim', [shift1-diff/3, shift2+diff/3]);
    set(handles.axes2,'xdir','rev');
else
    [XR, ppmR]=remove_ends(X, ppm, shift1, shift2);
    plot(handles.axes2, ppmR, XR);
    set(handles.axes2,'xlim', [shift1, shift2]);
    set(handles.axes2,'xdir','rev');
end

% --- Executes on button press in remove_right.
function remove_right_Callback(hObject, eventdata, handles)
% hObject    handle to remove_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
shift1=str2double(get(handles.edit3,'String'));
shift2=str2double(get(handles.edit4,'String'));
diff=shift2-shift1;
X=handles.X;
ppm=handles.ppm;
if (isequal(handles.region,'water_region') | isequal(handles.region,'other_region') )
    XR=remove_region(X, ppm, shift1, shift2);
    plot(handles.axes3, ppm, XR);
    set(handles.axes3,'xlim', [shift1-diff/3, shift2+diff/3]);
    set(handles.axes3,'xdir','rev');
else
    [XR, ppmR]=remove_ends(X, ppm, shift1, shift2);
    plot(handles.axes3, ppmR, XR);
    set(handles.axes3,'xlim', [shift1, shift2]);
    set(handles.axes3,'xdir','rev');
end


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in menu.
function menu_Callback(hObject, eventdata, handles)
% hObject    handle to menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Determine the selected data set.
str = get(hObject, 'String');
val = get(hObject,'Value');
disp(str{val});
switch str{val};
case 'water_region'
   handles.region = 'water_region';
   set(handles.edit1,'String', '4.4');
   set(handles.edit2,'String', '4.6');
   %handles.edit1=4.4;
   %handles.edit2=4.6;
case 'ends_region' 
   handles.region = 'ends_region';
   set(handles.edit1,'String', '-0.5');
   set(handles.edit2,'String', '11');
   %handles.edit1=-0.5;
   %handles.edit2=11;
case 'other_region' 
   handles.region = 'other_region';
end
% Save the handles structure.
guidata(hObject,handles)
% Hints: contents = cellstr(get(hObject,'String')) returns menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menu


% --- Executes during object creation, after setting all properties.
function menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
