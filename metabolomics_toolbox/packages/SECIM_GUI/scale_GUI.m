function varargout = scale_GUI(varargin)
% SCALE_GUI MATLAB code for scale_GUI.fig
%      SCALE_GUI, by itself, creates a new SCALE_GUI or raises the existing
%      singleton*.
%
%      H = SCALE_GUI returns the handle to a new SCALE_GUI or the handle to
%      the existing singleton*.
%
%      SCALE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCALE_GUI.M with the given input arguments.
%
%      SCALE_GUI('Property','Value',...) creates a new SCALE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before scale_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to scale_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help scale_GUI

% Last Modified by GUIDE v2.5 22-Oct-2014 11:07:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @scale_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @scale_GUI_OutputFcn, ...
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


% --- Executes just before scale_GUI is made visible.
function scale_GUI_OpeningFcn(hObject, eventdata, handles, parameters, X)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to scale_GUI (see VARARGIN)
handles.X=X;
handles.parameters= parameters;
% Choose default command line output for scale_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes scale_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = scale_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


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


% --- Executes on button press in scale1.
function scale1_Callback(hObject, eventdata, handles)
% hObject    handle to scale1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
reset(handles.axes1);
represent=get(handles.method1,'String');
represent=represent{get(handles.method1,'Value')};
handles.Xscaled1 = scale(handles.X,represent);
X=handles.Xscaled1;
features=1:size(X,2);
variance_X=var(X);
[h2,k2]=sort(features);
cla(handles.axes1);
plot(handles.axes1,h2,variance_X(k2));
xlabel(handles.axes1,'Feature - chemical shift or m/z')
ylabel(handles.axes1,'Variance')
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in PCA1.
function PCA1_Callback(hObject, eventdata, handles)
% hObject    handle to PCA1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
reset(handles.axes1);
%reset(handles.axes2);
X=handles.Xscaled1;
model= nipalsPCA(X,8);
%isequal(handles.Xscaled1,handles.Xscaled4)
components= str2num(get(handles.edit1,'string'));
if exist('Y')==0
    Y=zeros(size(X,1),1);
    Ycolor=ones(size(X,1),1);
else
Ycolor=ceil(([(Y-mean(Y))/(2.01*max(abs(Y-mean(Y))))]+.5)*100);
end
cmap=jet(100);

if length(components)==1
    %figure, 
    for k=1:size(model.scores,1);
        hold(handles.axes1,'on')
        plot(handles.axes1,model.scores(k,components(1)),Y(k),'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',cmap(Ycolor(k),:),'MarkerSize',8)
    end
    xlabel(handles.axes1,['Component',num2str(components(1)),', ', num2str(model.variance(components(1))*100),' Percent of Variance'])
    hold(handles.axes1,'off')
    
elseif length(components)==2
    %figure, 
    cla(handles.axes1);
    for k=1:size(model.scores,1);
        hold(handles.axes1,'on')
        plot(handles.axes1,model.scores(k,components(1)),model.scores(k,components(2)),'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',cmap(Ycolor(k),:),'MarkerSize',8)
    end
    xlabel(handles.axes1,['Component',num2str(components(1)),', ', num2str(model.variance(components(1))*100),' Percent of Variance'])
    ylabel(handles.axes1,['Component',num2str(components(2)),', ', num2str(model.variance(components(2))*100),' Percent of Variance'])
    hold(handles.axes1,'off')
    
elseif length(components)==3 
    %reset(handles.axes2);
    reset(handles.axes1);
    axes(handles.axes1);
    for k=1:size(model.scores,1);
        hold(handles.axes1,'on')
        plot3(handles.axes1,model.scores(k,components(1)),model.scores(k,components(2)),model.scores(k,components(3)),'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',cmap(Ycolor(k),:),'MarkerSize',8)
    end
    view(3)
    xlabel(handles.axes1,['Component',num2str(components(1)),', ', num2str(model.variance(components(1))*100),' Percent of Variance'])
    ylabel(handles.axes1,['Component',num2str(components(2)),', ', num2str(model.variance(components(2))*100),' Percent of Variance'])
    zlabel(handles.axes1,['Component',num2str(components(3)),', ', num2str(model.variance(components(3))*100),' Percent of Variance'])
   
else
    error('components must be a vector with 1, 2, or 3 elements')
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


% --- Executes on button press in scale2.
function scale2_Callback(hObject, eventdata, handles)
% hObject    handle to scale2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%reset(handles.axes1);
reset(handles.axes2);
represent=get(handles.method2,'String');
represent=represent{get(handles.method2,'Value')};
handles.Xscaled2 = scale(handles.X,represent);
X=handles.Xscaled2;
features=1:size(X,2);
variance_X=var(X);
[h2,k2]=sort(features);
cla(handles.axes2);
plot(handles.axes2,h2,variance_X(k2));
xlabel(handles.axes2,'Feature - chemical shift or m/z')
ylabel(handles.axes2,'Variance')
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in PCA2.
function PCA2_Callback(hObject, eventdata, handles)
% hObject    handle to PCA2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
X=handles.Xscaled2;
model= nipalsPCA(X,8);
%isequal(handles.Xscaled1,handles.Xscaled4)
components= str2num(get(handles.edit2,'string'));
if exist('Y')==0
    Y=zeros(size(X,1),1);
    Ycolor=ones(size(X,1),1);
else
Ycolor=ceil(([(Y-mean(Y))/(2.01*max(abs(Y-mean(Y))))]+.5)*100);
end
cmap=jet(100);

if length(components)==1
    %figure, 
    for k=1:size(model.scores,1);
        hold(handles.axes2,'on')
        plot(handles.axes2,model.scores(k,components(1)),Y(k),'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',cmap(Ycolor(k),:),'MarkerSize',8)
    end
    xlabel(handles.axes2,['Component',num2str(components(1)),', ', num2str(model.variance(components(1))*100),' Percent of Variance'])
   
elseif length(components)==2
    %figure, 
    cla(handles.axes2);
    for k=1:size(model.scores,1);
        hold(handles.axes2,'on')
        plot(handles.axes2,model.scores(k,components(1)),model.scores(k,components(2)),'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',cmap(Ycolor(k),:),'MarkerSize',8)
    end
    xlabel(handles.axes2,['Component',num2str(components(1)),', ', num2str(model.variance(components(1))*100),' Percent of Variance'])
    ylabel(handles.axes2,['Component',num2str(components(2)),', ', num2str(model.variance(components(2))*100),' Percent of Variance'])
    hold(handles.axes2,'off')
    
elseif length(components)==3
    %figure, 
    reset(handles.axes2);
    axes(handles.axes2);
    cla(handles.axes2);
    for k=1:size(model.scores,1);
        hold(handles.axes2,'on')
        plot3(handles.axes2,model.scores(k,components(1)),model.scores(k,components(2)),model.scores(k,components(3)),'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',cmap(Ycolor(k),:),'MarkerSize',8)
    end
    view(3)
    xlabel(handles.axes2,['Component',num2str(components(1)),', ', num2str(model.variance(components(1))*100),' Percent of Variance'])
    ylabel(handles.axes2,['Component',num2str(components(2)),', ', num2str(model.variance(components(2))*100),' Percent of Variance'])
    zlabel(handles.axes2,['Component',num2str(components(3)),', ', num2str(model.variance(components(3))*100),' Percent of Variance'])
   
else
    error('components must be a vector with 1, 2, or 3 elements')
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


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.radiobutton1, 'Value')==1
            represent = get(handles.method1,'String');
            represent = represent{get(handles.method1,'Value')};
            handles.parameters.T.scaling_method = represent;
           
elseif get(handles.radiobutton2, 'Value')==1
            represent = get(handles.method2,'String');
            represent = represent{get(handles.method2,'Value')};
            handles.parameters.T.scaling_method = represent;
end
