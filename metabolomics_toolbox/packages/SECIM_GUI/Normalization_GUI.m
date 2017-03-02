function varargout = Normalization_GUI(varargin)
% NORMALIZATION_GUI MATLAB code for Normalization_GUI.fig
%      NORMALIZATION_GUI, by itself, creates a new NORMALIZATION_GUI or raises the existing
%      singleton*.
%
%      H = NORMALIZATION_GUI returns the handle to a new NORMALIZATION_GUI or the handle to
%      the existing singleton*.
%
%      NORMALIZATION_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NORMALIZATION_GUI.M with the given input arguments.
%
%      NORMALIZATION_GUI('Property','Value',...) creates a new NORMALIZATION_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Normalization_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Normalization_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Normalization_GUI

% Last Modified by GUIDE v2.5 17-Oct-2014 12:03:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Normalization_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Normalization_GUI_OutputFcn, ...
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


% --- Executes just before Normalization_GUI is made visible.
function Normalization_GUI_OpeningFcn(hObject, eventdata, handles, parameters, ppm, X)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Normalization_GUI (see VARARGIN)
handles.X=X;
handles.ppm=ppm;
handles.parameters= parameters;
% Choose default command line output for Normalization_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Normalization_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Normalization_GUI_OutputFcn(hObject, eventdata, handles) 
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
            represent = get(handles.represent1,'String');
            represent = represent{get(handles.represent1,'Value')};
            handles.parameters.T.normalization_method = represent;
            if (isequal(represent,'Integral') | isequal(represent,'Intensity'))
                feature = get(handles.edit1,'String');
                handles.parameters.T.normalization_feature = feature;
            else
                handles.parameters.T.normalization_feature = 'NA';
            end
elseif get(handles.radiobutton2, 'Value')==1
            represent = get(handles.represent2,'String');
            represent = represent{get(handles.represent2,'Value')};
            handles.parameters.T.normalization_method = represent;
            if (isequal(represent,'Integral') | isequal(represent,'Intensity'))
                feature = get(handles.edit2,'String');
                handles.parameters.T.normalization_feature = feature;
            else
                handles.parameters.T.normalization_feature = 'NA';
            end
elseif get(handles.radiobutton3, 'Value')==1
            represent = get(handles.represent3,'String');
            represent = represent{get(handles.represent3,'Value')};
            handles.parameters.T.normalization_method = represent;
            if (isequal(represent,'Integral') | isequal(represent,'Intensity'))
                feature = get(handles.edit3,'String');
                handles.parameters.T.normalization_feature = feature;
            else
                handles.parameters.T.normalization_feature = 'NA';
            end
else
            represent = get(handles.represent4,'String');
            represent = represent{get(handles.represent4,'Value')};
            handles.parameters.T.normalization_method = represent;
            if (isequal(represent,'Integral') | isequal(represent,'Intensity'))
                feature = get(handles.edit4,'String');
                handles.parameters.T.normalization_feature = feature;
            else
                handles.parameters.T.normalization_feature = 'NA';
            end
end


% --- Executes on button press in align4.
function align4_Callback(hObject, eventdata, handles)
% hObject    handle to align4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
represent=get(handles.represent4,'String');
represent=represent{get(handles.represent4,'Value')};
if (isequal(represent,'Integral') | isequal(represent,'Intensity'))
    feature=str2num(get(handles.edit4,'string'));
    handles.Xnormalized4 = normalize(handles.X,handles.ppm,represent,feature);
else
    handles.Xnormalized4 = normalize(handles.X,handles.ppm,represent);
end
X=abs(handles.Xnormalized4);
F=X./repmat(median(X),[size(X,1),1]);
x=-4:.16:4;
n=histc(log(F)',x);
plot(handles.axes4,x,n);
set(handles.axes4,'xdir','rev');
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in align3.
function align3_Callback(hObject, eventdata, handles)
% hObject    handle to align3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
represent=get(handles.represent3,'String');
represent=represent{get(handles.represent3,'Value')};
if (isequal(represent,'Integral') | isequal(represent,'Intensity'))
    feature=str2num(get(handles.edit3,'string'));
    handles.Xnormalized3 = normalize(handles.X,handles.ppm,represent,feature);
else
    handles.Xnormalized3 = normalize(handles.X,handles.ppm,represent);
end
X=abs(handles.Xnormalized3);
F=X./repmat(median(X),[size(X,1),1]);
x=-4:.16:4;
n=histc(log(F)',x);
plot(handles.axes3,x,n);
set(handles.axes3,'xdir','rev');
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in align2.
function align2_Callback(hObject, eventdata, handles)
% hObject    handle to align2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
represent=get(handles.represent2,'String');
represent=represent{get(handles.represent2,'Value')};
if (isequal(represent,'Integral') | isequal(represent,'Intensity'))
    feature=str2num(get(handles.edit2,'string'));
    handles.Xnormalized2 = normalize(handles.X,handles.ppm,represent,feature);
else
    handles.Xnormalized2 = normalize(handles.X,handles.ppm,represent);
end
X=abs(handles.Xnormalized2);
F=X./repmat(median(X),[size(X,1),1]);
x=-4:.16:4;
n=histc(log(F)',x);
plot(handles.axes2,x,n);
set(handles.axes2,'xdir','rev');
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in align1.
function align1_Callback(hObject, eventdata, handles)
% hObject    handle to align1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
represent=get(handles.represent1,'String');
represent=represent{get(handles.represent1,'Value')};
if (isequal(represent,'Integral') | isequal(represent,'Intensity'))
    feature=str2num(get(handles.edit1,'string'));
    handles.Xnormalized1 = normalize(handles.X,handles.ppm,represent,feature);
else
    handles.Xnormalized1 = normalize(handles.X,handles.ppm,represent);
end
X=abs(handles.Xnormalized1);
F=X./repmat(median(X),[size(X,1),1]);
x=-4:.16:4;
n=histc(log(F)',x);
plot(handles.axes1,x,n);
set(handles.axes1,'xdir','rev');
% Update handles structure
guidata(hObject, handles);

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


% --- Executes on button press in boxplot4.
function boxplot4_Callback(hObject, eventdata, handles)
% hObject    handle to boxplot4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
X=abs(handles.Xnormalized4);
F=X./repmat(median(X),[size(X,1),1]);
x=-4:.16:4;
boxplot(handles.axes4,log(F)','plotstyle','compact','symbol',' ');
ylim(handles.axes4,[-4,4]);

% --- Executes on button press in boxplot3.
function boxplot3_Callback(hObject, eventdata, handles)
% hObject    handle to boxplot3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
X=abs(handles.Xnormalized3);
F=X./repmat(median(X),[size(X,1),1]);
x=-4:.16:4;
boxplot(handles.axes3,log(F)','plotstyle','compact','symbol',' ');
ylim(handles.axes3,[-4,4]);

% --- Executes on button press in boxplot2.
function boxplot2_Callback(hObject, eventdata, handles)
% hObject    handle to boxplot2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
X=abs(handles.Xnormalized2);
F=X./repmat(median(X),[size(X,1),1]);
x=-4:.16:4;
boxplot(handles.axes2,log(F)','plotstyle','compact','symbol',' ');
ylim(handles.axes2,[-4,4]);

% --- Executes on button press in boxplot1.
function boxplot1_Callback(hObject, eventdata, handles)
% hObject    handle to boxplot1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
X=abs(handles.Xnormalized1);
F=X./repmat(median(X),[size(X,1),1]);
x=-4:.16:4;
boxplot(handles.axes1,log(F)','plotstyle','compact','symbol',' ');
ylim(handles.axes1,[-4,4]);
