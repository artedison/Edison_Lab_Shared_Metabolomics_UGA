function varargout = alignment_GUI(varargin)
% ALIGNMENT_GUI MATLAB code for alignment_GUI.fig
%      ALIGNMENT_GUI, by itself, creates a new ALIGNMENT_GUI or raises the existing
%      singleton*.
%
%      H = ALIGNMENT_GUI returns the handle to a new ALIGNMENT_GUI or the handle to
%      the existing singleton*.
%
%      ALIGNMENT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ALIGNMENT_GUI.M with the given input arguments.
%
%      ALIGNMENT_GUI('Property','Value',...) creates a new ALIGNMENT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before alignment_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to alignment_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help alignment_GUI

% Last Modified by GUIDE v2.5 06-Sep-2014 21:55:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @alignment_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @alignment_GUI_OutputFcn, ...
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


% --- Executes just before alignment_GUI is made visible.
function alignment_GUI_OpeningFcn(hObject, eventdata, handles, parameters, ppm, X)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to alignment_GUI (see VARARGIN)
handles.X=X;
handles.ppm=ppm;
handles.parameters= parameters;
% Choose default command line output for alignment_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes alignment_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = alignment_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%ali%disp(class(get(handles.method1,'String')))
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in represent1.
function represent1_Callback(hObject, eventdata, handles)
% hObject    handle to represent1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns represent1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from represent1


% --- Executes during object creation, after setting all properties.
function represent1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to represent1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in method1.
function method1_Callback(hObject, eventdata, handles)
% hObject    handle to method1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns method1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from method1


% --- Executes during object creation, after setting all properties.
function method1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to method1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in algo1.
function algo1_Callback(hObject, eventdata, handles)
% hObject    handle to algo1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
names = get(hObject, 'String');
switch names{get(hObject,'Value')};
    case 'Guided'
        set(handles.represent1,'String', cellstr(char('correlation','spearman')));
    case 'Star'
        set(handles.represent1,'String', cellstr(char('mean','median','max','min','var')));
end            
% Hints: contents = cellstr(get(hObject,'String')) returns algo1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from algo1
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function algo1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to algo1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in represent2.
function represent2_Callback(hObject, eventdata, handles)
% hObject    handle to represent2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns represent2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from represent2


% --- Executes during object creation, after setting all properties.
function represent2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to represent2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in method2.
function method2_Callback(hObject, eventdata, handles)
% hObject    handle to method2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns method2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from method2


% --- Executes during object creation, after setting all properties.
function method2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to method2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in algo2.
function algo2_Callback(hObject, eventdata, handles)
% hObject    handle to algo2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
names = get(hObject, 'String');
switch names{get(hObject,'Value')};
    case 'Guided'
        set(handles.represent2,'String', cellstr(char('correlation','spearman')));
    case 'Star'
        set(handles.represent2,'String', cellstr(char('mean','median','max','min','var')));
end            
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns algo2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from algo2


% --- Executes during object creation, after setting all properties.
function algo2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to algo2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in represent3.
function represent3_Callback(hObject, eventdata, handles)
% hObject    handle to represent3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns represent3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from represent3


% --- Executes during object creation, after setting all properties.
function represent3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to represent3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in method3.
function method3_Callback(hObject, eventdata, handles)
% hObject    handle to method3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns method3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from method3


% --- Executes during object creation, after setting all properties.
function method3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to method3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in algo3.
function algo3_Callback(hObject, eventdata, handles)
% hObject    handle to algo3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
names = get(hObject, 'String');
switch names{get(hObject,'Value')};
    case 'Guided'
        set(handles.represent3,'String', cellstr(char('correlation','spearman')));
    case 'Star'
        set(handles.represent3,'String', cellstr(char('mean','median','max','min','var')));
end            
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns algo3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from algo3


% --- Executes during object creation, after setting all properties.
function algo3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to algo3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in represent4.
function represent4_Callback(hObject, eventdata, handles)
% hObject    handle to represent4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns represent4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from represent4


% --- Executes during object creation, after setting all properties.
function represent4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to represent4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in method4.
function method4_Callback(hObject, eventdata, handles)
% hObject    handle to method4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns method4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from method4


% --- Executes during object creation, after setting all properties.
function method4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to method4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in algo4.
function algo4_Callback(hObject, eventdata, handles)
% hObject    handle to algo4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
names = get(hObject, 'String');
switch names{get(hObject,'Value')};
    case 'Guided'
        set(handles.represent4,'String', cellstr(char('correlation','spearman')));
    case 'Star'
        set(handles.represent4,'String', cellstr(char('mean','median','max','min','var')));
end            
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns algo4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from algo4


% --- Executes during object creation, after setting all properties.
function algo4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to algo4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in finalize.
function finalize_Callback(hObject, eventdata, handles)
% hObject    handle to finalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.radiobutton1, 'Value')==1
            disp('ali');
            algo = get(handles.algo1,'String');
            algo = algo{get(handles.algo1,'Value')};
            handles.parameters.T.align_function = algo;
            method = get(handles.method1,'String');
            method = method{get(handles.method1,'Value')};
            handles.parameters.T.align_method = method;
            represent = get(handles.represent1,'String');
            represent = represent{get(handles.represent1,'Value')};
            handles.parameters.T.align_represent = represent;
elseif get(handles.radiobutton2, 'Value')==1
            algo = get(handles.algo2,'String');
            algo = algo{get(handles.algo2,'Value')};
            handles.parameters.T.align_function = algo;
            method = get(handles.method2,'String');
            method = method{get(handles.method2,'Value')};
            handles.parameters.T.align_method = method;
            represent = get(handles.represent2,'String');
            represent = represent{get(handles.represent2,'Value')};
            handles.parameters.T.align_represent = represent;
elseif get(handles.radiobutton3, 'Value')==1
            algo = get(handles.algo3,'String');
            algo = algo{get(handles.algo3,'Value')};
            handles.parameters.T.align_function = algo;
            method = get(handles.method3,'String');
            method = method{get(handles.method3,'Value')};
            handles.parameters.T.align_method = method;
            represent = get(handles.represent3,'String');
            represent = represent{get(handles.represent3,'Value')};
            handles.parameters.T.align_represent = represent;
else
            algo = get(handles.algo4,'String');
            algo = algo{get(handles.algo4,'Value')};
            handles.parameters.T.align_function = algo;
            method = get(handles.method4,'String');
            method = method {get(handles.method4,'Value')};
            handles.parameters.T.align_method = method;
            represent = get(handles.represent4,'String');
            represent = represent{get(handles.represent4,'Value')};
            handles.parameters.T.align_represent = represent;
end


% --- Executes on button press in align4.
function align4_Callback(hObject, eventdata, handles)
% hObject    handle to align4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Algo = get(handles.algo4,'String');
represent=get(handles.represent4,'String');
represent=represent{get(handles.represent4,'Value')};
alignment_method = get(handles.method4,'String');
alignment_method = alignment_method{get(handles.method4,'Value')};
switch Algo{get(handles.algo4,'Value')};
    case 'Guided'
        Xaligned4 = guide_align1D(handles.X,handles.ppm,represent,alignment_method);
    case 'Star'
        Xaligned4 = star_align1D(handles.X,handles.ppm, represent,alignment_method);
end;
plot(handles.axes4, handles.ppm, Xaligned4);
set(handles.axes4,'xdir','rev');

% --- Executes on button press in align3.
function align3_Callback(hObject, eventdata, handles)
% hObject    handle to align3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Algo = get(handles.algo3,'String');
represent=get(handles.represent3,'String');
represent=represent{get(handles.represent3,'Value')};
alignment_method = get(handles.method3,'String');
alignment_method = alignment_method{get(handles.method3,'Value')};
switch Algo{get(handles.algo3,'Value')};
    case 'Guided'
        Xaligned3 = guide_align1D(handles.X,handles.ppm,represent,alignment_method);
    case 'Star'
        Xaligned3 = star_align1D(handles.X,handles.ppm, represent,alignment_method);
end;
plot(handles.axes3, handles.ppm, Xaligned3);
set(handles.axes3,'xdir','rev');

% --- Executes on button press in align2.
function align2_Callback(hObject, eventdata, handles)
% hObject    handle to align2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Algo = get(handles.algo2,'String');
represent=get(handles.represent2,'String');
represent=represent{get(handles.represent2,'Value')};
alignment_method = get(handles.method2,'String');
alignment_method = alignment_method{get(handles.method2,'Value')};
switch Algo{get(handles.algo2,'Value')};
    case 'Guided'
        Xaligned2 = guide_align1D(handles.X,handles.ppm,represent,alignment_method);
    case 'Star'
        Xaligned2 = star_align1D(handles.X,handles.ppm, represent,alignment_method);
end;
plot(handles.axes2, handles.ppm, Xaligned2);
set(handles.axes2,'xdir','rev');

% --- Executes on button press in align1.
function align1_Callback(hObject, eventdata, handles)
% hObject    handle to align1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Algo = get(handles.algo1,'String');
represent=get(handles.represent1,'String');
represent=represent{get(handles.represent1,'Value')};
alignment_method = get(handles.method1,'String');
alignment_method = alignment_method{get(handles.method1,'Value')};
switch Algo{get(handles.algo1,'Value')};
    case 'Guided'
        Xaligned1 = guide_align1D(handles.X,handles.ppm,represent,alignment_method);
    case 'Star'
        Xaligned1 = star_align1D(handles.X,handles.ppm, represent,alignment_method);
end;
plot(handles.axes1, handles.ppm, Xaligned1);
set(handles.axes1,'xdir','rev');
