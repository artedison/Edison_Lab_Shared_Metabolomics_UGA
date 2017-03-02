function varargout = ViewRawFCS(varargin)
% VIEWRAWFCS MATLAB code for ViewRawFCS.fig
%      VIEWRAWFCS, by itself, creates a new VIEWRAWFCS or raises the existing
%      singleton*.
%
%      H = VIEWRAWFCS returns the handle to a new VIEWRAWFCS or the handle to
%      the existing singleton*.
%
%      VIEWRAWFCS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEWRAWFCS.M with the given input arguments.
%
%      VIEWRAWFCS('Property','Value',...) creates a new VIEWRAWFCS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ViewRawFCS_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ViewRawFCS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ViewRawFCS

% Last Modified by GUIDE v2.5 28-Feb-2012 10:28:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ViewRawFCS_OpeningFcn, ...
                   'gui_OutputFcn',  @ViewRawFCS_OutputFcn, ...
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


% --- Executes just before ViewRawFCS is made visible.
function ViewRawFCS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ViewRawFCS (see VARARGIN)

% Choose default command line output for ViewRawFCS
handles.output = hObject;
handles.directoryname=[];
handles.fcs_filename = [];
handles.raw_data = [];
handles.transformed_data = [];
handles.marker_names = cell(0);
handles.transformation_option = 1; % 0 means no transformation, 1 means arcsinh, 2 means arcsinh followed by 0-mean 1-var
handles.arcsinh_cofactor = 5;
set(handles.radiobutton_1_arcsinh,'value',1);
set(handles.edit_arcsinh_cofactor,'string',5);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ViewRawFCS wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ViewRawFCS_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_fcs_filename_Callback(hObject, eventdata, handles)
% stop user from editing the working directory manually
set(handles.edit_fcs_filename,'String',handles.fcs_filename);


% --- Executes during object creation, after setting all properties.
function edit_fcs_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fcs_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_Browse.
function pushbutton_Browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('*.fcs', 'Pick one FCS file');
% check whether this button needs to do anything.
if isequal(filename,0) || isequal(pathname,0) % the user may click cancel in the selection dialog box
    fprintf('No fcs file selected!!\n\n');
    if isempty(handles.fcs_filename)
        fprintf('Current director is empty.\n\n');
    else
        fprintf('Current director is: %s\n\n',fullfile(handles.directoryname,handles.fcs_filename));
    end
    set(handles.edit_fcs_filename,'String',fullfile(handles.directoryname,handles.fcs_filename));
    return
end

% initialize analysis of this file
handles.directoryname = pathname;
handles.fcs_filename = filename;
set(handles.edit_fcs_filename,'String',fullfile(handles.directoryname,handles.fcs_filename));
[handles.raw_data, handles.marker_names] = readfcs(fullfile(pathname, filename));
handles.transformed_data = transform_data(handles.raw_data, handles.transformation_option, handles.arcsinh_cofactor);
guidata(hObject, handles);
% % biaxial plots
set(handles.Popup_x,'string',handles.marker_names,'value',1);
set(handles.Popup_y,'string',handles.marker_names,'value',2);
update_baxial_plots(handles);



% --- Executes on button press in radiobutton_0_no_transformation.
function radiobutton_0_no_transformation_Callback(hObject, eventdata, handles)
if handles.transformation_option ==0
    set(handles.radiobutton_0_no_transformation,'value',1);
    return
end
handles.transformation_option = 0;
handles.transformed_data = transform_data(handles.raw_data, handles.transformation_option, handles.arcsinh_cofactor);
guidata(hObject, handles);
update_baxial_plots(handles);


% --- Executes on button press in radiobutton_1_arcsinh.
function radiobutton_1_arcsinh_Callback(hObject, eventdata, handles)
if handles.transformation_option ==1
    set(handles.radiobutton_1_arcsinh,'value',1);
    return
end
handles.transformation_option = 1;
handles.transformed_data = transform_data(handles.raw_data, handles.transformation_option, handles.arcsinh_cofactor);
guidata(hObject, handles);
update_baxial_plots(handles);


% --- Executes on button press in radiobutton_2_arcsinh_norm.
function radiobutton_2_arcsinh_norm_Callback(hObject, eventdata, handles)
if handles.transformation_option ==2
    set(handles.radiobutton_2_arcsinh_norm,'value',1);
    return
end
handles.transformation_option = 2;
handles.transformed_data = transform_data(handles.raw_data, handles.transformation_option, handles.arcsinh_cofactor);
guidata(hObject, handles);
update_baxial_plots(handles);



function edit_arcsinh_cofactor_Callback(hObject, eventdata, handles)

tmp = get(handles.edit_arcsinh_cofactor,'string');
new_cofactor = str2num(tmp);
if isempty(str2num(tmp)) || new_cofactor<=0
	set(handles.edit_arcsinh_cofactor,'value',handles.arcsinh_cofactor,'string',num2str(handles.arcsinh_cofactor));
else
	handles.arcsinh_cofactor = new_cofactor;
    guidata(hObject,handles); 
    if handles.transformation_option==1 || handles.transformation_option==2
        handles.transformed_data = transform_data(handles.raw_data, handles.transformation_option, handles.arcsinh_cofactor);
        guidata(hObject,handles); 
        update_baxial_plots(handles);
    end
end


% --- Executes during object creation, after setting all properties.
function edit_arcsinh_cofactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_arcsinh_cofactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function [transformed_data] = transform_data(raw_data, transformation_option, arcsinh_cofactor)

fprintf('Data transformation options in SPADE parameters ... ');
switch transformation_option
    case 0, 
        fprintf('No transofrmation performed\n');
        transformed_data = raw_data;
    case 1, 
        fprintf(['arcsinh transformation with cofactor ',num2str(arcsinh_cofactor),'\n']);
        transformed_data = flow_arcsinh(raw_data,arcsinh_cofactor); 
    case 2, 
        fprintf(['arcsinh transformation with cofactor ',num2str(arcsinh_cofactor),', followed by 0-mean-1-var normalization \n']);
        transformed_data = SPADE_per_gene_normalization(flow_arcsinh(raw_data,arcsinh_cofactor)); 
    otherwise, 1;
end


% --- Executes on selection change in Popup_y.
function Popup_y_Callback(hObject, eventdata, handles)
update_baxial_plots(handles);


% --- Executes during object creation, after setting all properties.
function Popup_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Popup_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Popup_x.
function Popup_x_Callback(hObject, eventdata, handles)
update_baxial_plots(handles);


% --- Executes during object creation, after setting all properties.
function Popup_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Popup_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% my own function to draw the two biaxial plots
function update_baxial_plots(handles)
tmp = get(handles.Popup_x,'string');
marker_x = tmp(get(handles.Popup_x,'value'));
tmp = get(handles.Popup_y,'string');
marker_y = tmp(get(handles.Popup_y,'value'));
x = handles.transformed_data(ismember(handles.marker_names,marker_x),:);
y = handles.transformed_data(ismember(handles.marker_names,marker_y),:);
ind = ((~isnan(x)) & (~isnan(y)));
x = x(ind);
y = y(ind);
axes(handles.Axes_scatter);
plot(x,y,'g.');
axis([min(x)-(max(x)-min(x))*0.05,max(x)+(max(x)-min(x))*0.05,min(y)-(max(y)-min(y))*0.05,max(y)+(max(y)-min(y))*0.05])
axis_lim = axis;
axes(handles.Axes_contour);
SPADE_contour2D(x,y);
axis(axis_lim);



function SPADE_contour2D(x,y,contour_levels,outlier_percentage)
% FlowJo_contour2D(x,y, contour_levels, outlier_percentage)
% outlier_percentage = 1, 2, 5, 10, ...
% example:
%       x = randn(1,50000) + [zeros(1,20000), ones(1,30000)]; 
%       y = randn(1,50000)+ [3*ones(1,20000),zeros(1,30000)];
%       outlier_percentage = 5;
%       FlowJo_contour2D(x,y, [0:15:60,70:10:90,93:3:99])

if exist('contour_levels')~=1
    contour_levels = [0:15:60,70:10:90,93:3:99];
end

if exist('outlier_percentage')~=1
    outlier_percentage=5;
end

if length(x)>100000
    ind = randsample(1:length(x),100000);
    x_downsample = x(ind); 
    y_downsample = y(ind);
else
    x_downsample = x;      
    y_downsample = y;
end
display('computing ksdensity grids for each variable ...')
[F,XI,hx] = ksdensity(x_downsample); 
[F,YI,hy] = ksdensity(y_downsample); 
Z = zeros(length(XI),length(YI));

display('computing the gaussian kernel density of each point on the grid ...')
dist_XI_x = exp(-(XI(:)*ones(1,length(x_downsample)) - ones(length(XI),1)*(x_downsample(:)')).^2/(2*hx^2));
dist_YI_y = exp(-(YI(:)*ones(1,length(y_downsample)) - ones(length(YI),1)*(y_downsample(:)')).^2/(2*hy^2));
Z = (dist_XI_x*dist_YI_y')';
Z = Z./(sqrt(2*pi*hx^2*2*pi*hy)*length(x_downsample));


Y = sort(Z(:));
[dummy,ind] = min(abs(cumsum(Y)-sum(Y)*outlier_percentage/100));
outer_ring = Y(ind);

display('drawing the contour plot ...')
contour(XI,YI,Z,prctile(Y(ind:end),contour_levels));
axis_tmp = axis;

return