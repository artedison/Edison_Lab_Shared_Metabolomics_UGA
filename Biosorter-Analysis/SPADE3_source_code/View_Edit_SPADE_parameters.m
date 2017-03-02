function varargout = View_Edit_SPADE_parameters(varargin)
% VIEW_EDIT_SPADE_PARAMETERS M-file for View_Edit_SPADE_parameters.fig
%      VIEW_EDIT_SPADE_PARAMETERS, by itself, creates a new VIEW_EDIT_SPADE_PARAMETERS or raises the existing
%      singleton*.
%
%      H = VIEW_EDIT_SPADE_PARAMETERS returns the handle to a new VIEW_EDIT_SPADE_PARAMETERS or the handle to
%      the existing singleton*.
%
%      VIEW_EDIT_SPADE_PARAMETERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEW_EDIT_SPADE_PARAMETERS.M with the given input arguments.
%
%      VIEW_EDIT_SPADE_PARAMETERS('Property','Value',...) creates a new VIEW_EDIT_SPADE_PARAMETERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before View_Edit_SPADE_parameters_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to View_Edit_SPADE_parameters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help View_Edit_SPADE_parameters

% Last Modified by GUIDE v2.5 07-Mar-2013 23:59:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @View_Edit_SPADE_parameters_OpeningFcn, ...
                   'gui_OutputFcn',  @View_Edit_SPADE_parameters_OutputFcn, ...
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


% --- Executes just before View_Edit_SPADE_parameters is made visible.
function View_Edit_SPADE_parameters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to View_Edit_SPADE_parameters (see VARARGIN)

% Choose default command line output for View_Edit_SPADE_parameters
handles.output = hObject;
handles.mother_window_handles = varargin{1};
guidata(hObject, handles);
% initialize the listbox of used and not-used markers
used_markers_str = handles.mother_window_handles.used_markers;
[C,I] = setdiff(handles.mother_window_handles.all_overlapping_markers, used_markers_str);
not_used_markers_str = handles.mother_window_handles.all_overlapping_markers(sort(I));
set(handles.listbox_overlapping_markers_not_used,'string',not_used_markers_str);
set(handles.listbox_overlapping_markers_used,'string',used_markers_str);
% initialize the part about compensation
if ~isfield(handles.mother_window_handles,'apply_compensation')
    handles.mother_window_handles.apply_compensation=0;
    guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
end
if handles.mother_window_handles.apply_compensation==0 % no compensation
    set(handles.radiobutton_compensation_ignore,'value',1)
    set(handles.radiobutton_compensation_apply,'value',0)
else % yes compensation
    set(handles.radiobutton_compensation_ignore,'value',0)
    set(handles.radiobutton_compensation_apply,'value',1)
end
% initialize the transformation options
set(handles.edit_arcsinh_cofactor,'value',handles.mother_window_handles.arcsinh_cofactor,'string',num2str(handles.mother_window_handles.arcsinh_cofactor));
switch handles.mother_window_handles.transformation_option
    case 0    
        set(handles.radiobutton_0_no_transformation,'value',1);
    case 1    
        set(handles.radiobutton_1_arcsinh,'value',1);
    case 2    
        set(handles.radiobutton_2_arcsinh_norm,'value',1);
    otherwise
        set(handles.radiobutton_1_arcsinh,'value',1);
end
% density estimation parameters
set(handles.edit_kernel_width_factor,'string', num2str(handles.mother_window_handles.kernel_width_factor));
set(handles.edit_density_estimation_optimization_factor,'string', num2str(handles.mother_window_handles.density_estimation_optimization_factor));
% downsampling OD and TD
set(handles.edit_outlier_density,'string',num2str(handles.mother_window_handles.outlier_density));
set(handles.edit_target_density,'string',num2str(handles.mother_window_handles.target_density));
set(handles.edit_target_cell_number, 'string', num2str(handles.mother_window_handles.target_cell_number));
switch handles.mother_window_handles.target_density_mode
    case 1    
        set(handles.radiobutton_1_TD_percentile,'value',1);
    case 2    
        set(handles.radiobutton_2_TD_cell_number,'value',1);
    otherwise
        set(handles.radiobutton_2_TD_cell_number,'value',1);
end
% number of desired clusters
set(handles.edit_desired_number_of_clusters,'string',num2str(handles.mother_window_handles.number_of_desired_clusters));
set(handles.edit_max_allowable_events,'string',num2str(handles.mother_window_handles.max_allowable_events));
% clustering algorithm
switch handles.mother_window_handles.clustering_algorithm
    case 'kmeans'    
        set(handles.radiobutton_4_kmeans,'value',1);
        set(handles.radiobutton_5_agglomerative,'value',0);
    case 'agglomerative'    
        set(handles.radiobutton_4_kmeans,'value',0);
        set(handles.radiobutton_5_agglomerative,'value',1);
    otherwise
        set(handles.radiobutton_4_kmeans,'value',1);
        set(handles.radiobutton_5_agglomerative,'value',0);
        handles.mother_window_handles.clustering_algorithm = 'kmeans';
        guidata(hObject,handles); 
        guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
end
% files used or not used to build the SPADE tree
files_used = handles.mother_window_handles.file_used_to_build_SPADE_tree;
[C,I] = setdiff(handles.mother_window_handles.file_annot, files_used);
files_not_used = handles.mother_window_handles.file_annot(sort(I));
set(handles.listbox_files_not_used,'string',files_not_used);
set(handles.listbox_files_used,'string',files_used);

% UIWAIT makes View_Edit_SPADE_parameters wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = View_Edit_SPADE_parameters_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --- Executes on selection change in listbox_overlapping_markers_not_used.
function listbox_overlapping_markers_not_used_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_overlapping_markers_not_used (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_overlapping_markers_not_used contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_overlapping_markers_not_used


% --- Executes during object creation, after setting all properties.
function listbox_overlapping_markers_not_used_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_overlapping_markers_not_used (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_overlapping_markers_used.
function listbox_overlapping_markers_used_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_overlapping_markers_used (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_overlapping_markers_used contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_overlapping_markers_used


% --- Executes during object creation, after setting all properties.
function listbox_overlapping_markers_used_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_overlapping_markers_used (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in button_add_marker.
function button_add_marker_Callback(hObject, eventdata, handles)

not_used_markers_str = get(handles.listbox_overlapping_markers_not_used,'String');
used_markers_str = get(handles.listbox_overlapping_markers_used,'String');
ind = get(handles.listbox_overlapping_markers_not_used,'value');
if length(not_used_markers_str)==0 || length(not_used_markers_str)<ind
    return
end 
if ind == length(not_used_markers_str)
    set(handles.listbox_overlapping_markers_not_used,'value',ind-1);
end
tmp = not_used_markers_str(ind);  % this is the marker to be moved over
not_used_markers_str = not_used_markers_str(setdiff(1:end,ind)); 
[C,IA,IB] = intersect(handles.mother_window_handles.all_overlapping_markers, [used_markers_str(:);tmp]);
used_markers_str = handles.mother_window_handles.all_overlapping_markers(sort(IA));

set(handles.listbox_overlapping_markers_not_used,'String',not_used_markers_str);
if get(handles.listbox_overlapping_markers_used,'value')==0
    set(handles.listbox_overlapping_markers_used,'value',1);
end
set(handles.listbox_overlapping_markers_used,'String',used_markers_str);
handles.mother_window_handles.used_markers = used_markers_str;
guidata(hObject,handles); 
guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
% update the parameter file 
parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
used_markers = handles.mother_window_handles.used_markers;
save(parameter_filename,'used_markers','-append');
% update the .mat files
delete_LD_PooledData_Result_mat_files(handles)



% --- Executes on button press in button_remove_marker.
function button_remove_marker_Callback(hObject, eventdata, handles)
not_used_markers_str = get(handles.listbox_overlapping_markers_not_used,'String');
used_markers_str = get(handles.listbox_overlapping_markers_used,'String');
ind = get(handles.listbox_overlapping_markers_used,'value');
if length(used_markers_str)==0 || length(used_markers_str)<ind
    return
end 
if ind == length(used_markers_str)
    set(handles.listbox_overlapping_markers_used,'value',ind-1);
end
tmp = used_markers_str(ind);  % this is the marker to be moved over
used_markers_str = used_markers_str(setdiff(1:end,ind)); 
[C,IA,IB] = intersect(handles.mother_window_handles.all_overlapping_markers, [not_used_markers_str(:);tmp]);
not_used_markers_str = handles.mother_window_handles.all_overlapping_markers(sort(IA));

set(handles.listbox_overlapping_markers_used,'String',used_markers_str);
if get(handles.listbox_overlapping_markers_not_used,'value')==0
    set(handles.listbox_overlapping_markers_not_used,'value',1);
end
set(handles.listbox_overlapping_markers_not_used,'String',not_used_markers_str);
handles.mother_window_handles.used_markers = used_markers_str;
guidata(hObject,handles); 
guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
% update the parameter file 
parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
used_markers = handles.mother_window_handles.used_markers;
save(parameter_filename,'used_markers','-append');
% update the .mat files
delete_LD_PooledData_Result_mat_files(handles)





function edit_arcsinh_cofactor_Callback(hObject, eventdata, handles)
% hObject    handle to edit_arcsinh_cofactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = get(handles.edit_arcsinh_cofactor,'string');
new_cofactor = str2num(tmp);
if isempty(str2num(tmp)) || new_cofactor<=0
	set(handles.edit_arcsinh_cofactor,'value',handles.mother_window_handles.arcsinh_cofactor,'string',num2str(handles.mother_window_handles.arcsinh_cofactor));
else
	handles.mother_window_handles.arcsinh_cofactor = new_cofactor;
    guidata(hObject,handles); 
    guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
    % update the parameter file 
    parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
    arcsinh_cofactor = handles.mother_window_handles.arcsinh_cofactor;
    save(parameter_filename,'arcsinh_cofactor','-append');
    % update the .mat files
    delete_LD_PooledData_Result_mat_files(handles)
end


% --- Executes during object creation, after setting all properties.
function edit_arcsinh_cofactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_arcsinh_cofactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in radiobutton_0_no_transformation.
function radiobutton_0_no_transformation_Callback(hObject, eventdata, handles)
% set(handles.radiobutton_0_no_transformation,'value',1);
if handles.mother_window_handles.transformation_option==0
    set(handles.radiobutton_0_no_transformation,'value',1);
    return
end
handles.mother_window_handles.transformation_option=0;
guidata(hObject,handles); 
guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
% update the parameter file 
parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
transformation_option = handles.mother_window_handles.transformation_option;
save(parameter_filename,'transformation_option','-append');
% update the .mat files
delete_LD_PooledData_Result_mat_files(handles)



% --- Executes on button press in radiobutton_1_arcsinh.
function radiobutton_1_arcsinh_Callback(hObject, eventdata, handles)
% set(handles.radiobutton_1_arcsinh,'value',1);
if handles.mother_window_handles.transformation_option==1
    set(handles.radiobutton_1_arcsinh,'value',1);
    return
end
handles.mother_window_handles.transformation_option=1;
guidata(hObject,handles); 
guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
% update the parameter file 
parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
transformation_option = handles.mother_window_handles.transformation_option;
save(parameter_filename,'transformation_option','-append');
% update the .mat files
delete_LD_PooledData_Result_mat_files(handles)




% --- Executes on button press in radiobutton_2_arcsinh_norm.
function radiobutton_2_arcsinh_norm_Callback(hObject, eventdata, handles)
% set(handles.radiobutton_2_arcsinh_norm,'value',1);
if handles.mother_window_handles.transformation_option==2
    set(handles.radiobutton_2_arcsinh_norm,'value',1);
    return
end
handles.mother_window_handles.transformation_option=2;
guidata(hObject,handles); 
guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
% update the parameter file 
parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
transformation_option = handles.mother_window_handles.transformation_option;
save(parameter_filename,'transformation_option','-append');
% update the .mat files
delete_LD_PooledData_Result_mat_files(handles)




function edit_kernel_width_factor_Callback(hObject, eventdata, handles)
% hObject    handle to edit_kernel_width_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = get(handles.edit_kernel_width_factor,'string');
new_factor = str2num(tmp);
if isempty(str2num(tmp)) || new_factor<=0 
	set(handles.edit_kernel_width_factor,'string',num2str(handles.mother_window_handles.kernel_width_factor));
else
	handles.mother_window_handles.kernel_width_factor = new_factor;
    guidata(hObject,handles); 
    guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
    % update the parameter file 
    parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
    kernel_width_factor = handles.mother_window_handles.kernel_width_factor;
    save(parameter_filename,'kernel_width_factor','-append');
    % update the .mat files
    delete_LD_PooledData_Result_mat_files(handles)
end


% --- Executes during object creation, after setting all properties.
function edit_kernel_width_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_kernel_width_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_density_estimation_optimization_factor_Callback(hObject, eventdata, handles)
% hObject    handle to edit_density_estimation_optimization_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = get(handles.edit_density_estimation_optimization_factor,'string');
new_factor = str2num(tmp);
if isempty(str2num(tmp)) || new_factor<=0 
	set(handles.edit_density_estimation_optimization_factor,'string',num2str(handles.mother_window_handles.density_estimation_optimization_factor));
else
	handles.mother_window_handles.density_estimation_optimization_factor = new_factor;
    guidata(hObject,handles); 
    guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
    % update the parameter file 
    parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
    density_estimation_optimization_factor = handles.mother_window_handles.density_estimation_optimization_factor;
    save(parameter_filename,'density_estimation_optimization_factor','-append');
    % update the .mat files
    delete_LD_PooledData_Result_mat_files(handles)
end

% --- Executes during object creation, after setting all properties.
function edit_density_estimation_optimization_factor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function edit_outlier_density_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_outlier_density_Callback(hObject, eventdata, handles)
tmp = get(handles.edit_outlier_density,'string');
new_outlier_density = str2num(tmp);
if isempty(str2num(tmp)) || new_outlier_density<0 || new_outlier_density>100
	set(handles.edit_outlier_density,'string',num2str(handles.mother_window_handles.outlier_density));
else
	handles.mother_window_handles.outlier_density = new_outlier_density;
    guidata(hObject,handles); 
    guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
    % update the parameter file 
    parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
    outlier_density = handles.mother_window_handles.outlier_density;
    save(parameter_filename,'outlier_density','-append');
    % updae .mat files
    delete_PooledData_Result_mat_files(handles)
end




% --- Executes during object creation, after setting all properties.
function edit_target_density_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_target_density_Callback(hObject, eventdata, handles)
tmp = get(handles.edit_target_density,'string');
new_target_density = str2num(tmp);
if isempty(str2num(tmp)) || new_target_density<=0 || new_target_density>100
	set(handles.edit_target_density,'string',num2str(handles.mother_window_handles.target_density));
else
	handles.mother_window_handles.target_density = new_target_density;
    guidata(hObject,handles); 
    guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
    % update the parameter file 
    parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
    target_density = handles.mother_window_handles.target_density;
    save(parameter_filename,'target_density','-append');
    % updae .mat files
    delete_PooledData_Result_mat_files(handles)
end



% --- Executes during object creation, after setting all properties.
function edit_target_cell_number_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_target_cell_number_Callback(hObject, eventdata, handles)
tmp = get(handles.edit_target_cell_number,'string');
new_target_cell_number = str2num(tmp);
if isempty(str2num(tmp)) || new_target_cell_number<=0 || new_target_cell_number~=round(new_target_cell_number)
	set(handles.edit_target_cell_number,'string',num2str(handles.mother_window_handles.target_cell_number));
else
	handles.mother_window_handles.target_cell_number = new_target_cell_number;
    guidata(hObject,handles); 
    guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
    % update the parameter file 
    parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
    target_cell_number = handles.mother_window_handles.target_cell_number;
    save(parameter_filename,'target_cell_number','-append');
    % updae .mat files
    delete_PooledData_Result_mat_files(handles)
end


% --- Executes on button press in radiobutton_1_TD_percentile.
function radiobutton_1_TD_percentile_Callback(hObject, eventdata, handles)
if handles.mother_window_handles.target_density_mode==1
    set(handles.radiobutton_1_TD_percentile,'value',1);
    return
end
handles.mother_window_handles.target_density_mode=1;
guidata(hObject,handles); 
guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
% update the parameter file 
parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
target_density_mode = handles.mother_window_handles.target_density_mode;
save(parameter_filename,'target_density_mode','-append');
% updae .mat files
delete_PooledData_Result_mat_files(handles)




% --- Executes on button press in radiobutton_2_TD_cell_number.
function radiobutton_2_TD_cell_number_Callback(hObject, eventdata, handles)
if handles.mother_window_handles.target_density_mode==2
    set(handles.radiobutton_2_TD_cell_number,'value',1);
    return
end
handles.mother_window_handles.target_density_mode=2;
guidata(hObject,handles); 
guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
% update the parameter file 
parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
target_density_mode = handles.mother_window_handles.target_density_mode;
save(parameter_filename,'target_density_mode','-append');
% updae .mat files
delete_PooledData_Result_mat_files(handles)




% --- Executes during object creation, after setting all properties.
function edit_desired_number_of_clusters_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_desired_number_of_clusters_Callback(hObject, eventdata, handles)
tmp = get(handles.edit_desired_number_of_clusters,'string');
new_desired_number_of_clusters = str2num(tmp);
if isempty(str2num(tmp)) || new_desired_number_of_clusters<=0 || new_desired_number_of_clusters~=round(new_desired_number_of_clusters)
	set(handles.edit_desired_number_of_clusters,'string',num2str(handles.mother_window_handles.number_of_desired_clusters));
else
	handles.mother_window_handles.number_of_desired_clusters = new_desired_number_of_clusters;
    guidata(hObject,handles); 
    guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
    % update the parameter file 
    parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
    number_of_desired_clusters = handles.mother_window_handles.number_of_desired_clusters;
    save(parameter_filename,'number_of_desired_clusters','-append');
    % updae .mat files
    delete_Result_mat_files(handles)
end



% --- Executes on selection change in listbox_files_not_used.
function listbox_files_not_used_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_files_not_used (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_files_not_used contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_files_not_used


% --- Executes during object creation, after setting all properties.
function listbox_files_not_used_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_files_not_used (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_files_used.
function listbox_files_used_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_files_used (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_files_used contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_files_used


% --- Executes during object creation, after setting all properties.
function listbox_files_used_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_files_used (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_add_file.
function button_add_file_Callback(hObject, eventdata, handles)

not_used_files_str = get(handles.listbox_files_not_used,'String');
used_files_str = get(handles.listbox_files_used,'String');
ind = get(handles.listbox_files_not_used,'value');
if length(not_used_files_str)==0 || length(not_used_files_str)<ind
    return
end 
if ind == length(not_used_files_str)
    set(handles.listbox_files_not_used,'value',ind-1);
end
tmp = not_used_files_str(ind);  % this is the flie to be moved over
not_used_files_str = not_used_files_str(setdiff(1:end,ind)); 
[C,IA,IB] = intersect(handles.mother_window_handles.file_annot, [used_files_str(:);tmp]);
used_files_str = handles.mother_window_handles.file_annot(sort(IA));

set(handles.listbox_files_not_used,'String',not_used_files_str);
if get(handles.listbox_files_used,'value')==0
    set(handles.listbox_files_used,'value',1);
end
set(handles.listbox_files_used,'String',used_files_str);
handles.mother_window_handles.file_used_to_build_SPADE_tree = used_files_str;
guidata(hObject,handles); 
guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
% update the parameter file 
parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
file_used_to_build_SPADE_tree = handles.mother_window_handles.file_used_to_build_SPADE_tree;
save(parameter_filename,'file_used_to_build_SPADE_tree','-append');
% updae .mat files
delete_PooledData_Result_mat_files(handles)


% --- Executes on button press in button_remove_file.
function button_remove_file_Callback(hObject, eventdata, handles)
used_files_str = get(handles.listbox_files_used,'String');
not_used_files_str = get(handles.listbox_files_not_used,'String');
ind = get(handles.listbox_files_used,'value');
if length(used_files_str)==0 || length(used_files_str)<ind
    return
end 
if ind == length(used_files_str)
    set(handles.listbox_files_used,'value',ind-1);
end
tmp = used_files_str(ind);  % this is the flie to be moved over
used_files_str = used_files_str(setdiff(1:end,ind)); 
[C,IA,IB] = intersect(handles.mother_window_handles.file_annot, [not_used_files_str(:);tmp]);
not_used_files_str = handles.mother_window_handles.file_annot(sort(IA));

set(handles.listbox_files_used,'String',used_files_str);
if get(handles.listbox_files_not_used,'value')==0
    set(handles.listbox_files_not_used,'value',1);
end
set(handles.listbox_files_not_used,'String',not_used_files_str);
handles.mother_window_handles.file_used_to_build_SPADE_tree = used_files_str;
guidata(hObject,handles); 
guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
% update the parameter file 
parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
file_used_to_build_SPADE_tree = handles.mother_window_handles.file_used_to_build_SPADE_tree;
save(parameter_filename,'file_used_to_build_SPADE_tree','-append');
% updae .mat files
delete_PooledData_Result_mat_files(handles)





function edit_max_allowable_events_Callback(hObject, eventdata, handles)
tmp = get(handles.edit_max_allowable_events,'string');
new_max_allowable_events = str2num(tmp);
if isempty(str2num(tmp)) || new_max_allowable_events<=0 || new_max_allowable_events~=round(new_max_allowable_events)
	set(handles.edit_max_allowable_events,'string',num2str(handles.mother_window_handles.max_allowable_events));
else
	handles.mother_window_handles.max_allowable_events = new_max_allowable_events;
    guidata(hObject,handles); 
    guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
    % update the parameter file 
    parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
    max_allowable_events = handles.mother_window_handles.max_allowable_events;
    save(parameter_filename,'max_allowable_events','-append');
    % updae .mat files
    delete_PooledData_Result_mat_files(handles)
end


% --- Executes during object creation, after setting all properties.
function edit_max_allowable_events_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% edit_arcsinh_cofactor_Callback(handles.edit_arcsinh_cofactor, [], handles);
% edit_kernel_width_factor_Callback(handles.edit_kernel_width_factor, [], handles);
% edit_density_estimation_optimization_factor_Callback(handles.edit_density_estimation_optimization_factor, [], handles);
% edit_outlier_density_Callback(handles.edit_outlier_density, [], handles);
% edit_target_density_Callback(handles.edit_target_density, [], handles);
% edit_target_cell_number_Callback(handles.edit_target_cell_number, [], handles);
% edit_desired_number_of_clusters_Callback(handles.edit_desired_number_of_clusters, [], handles);
% edit_max_allowable_events_Callback(handles.edit_max_allowable_events, [], handles);
% delete(hObject);


% --- Executes on button press in button_close_parameter_window.
function button_close_parameter_window_Callback(hObject, eventdata, handles)
delete(handles.figure1);
% figure1_CloseRequestFcn(handles.figure1, [], handles);




% parameter update induced change of files. 
function delete_LD_PooledData_Result_mat_files(handles)
for i=1:length(handles.mother_window_handles.all_fcs_filenames)
    fcs_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.all_fcs_filenames{i});
    mat_filename = [fcs_filename(1:end-3),'mat'];
    if exist(mat_filename)==2 % if the .mat file that stores the downsampling info already exist
        delete(mat_filename);
    end
end
pooled_downsampled_filename = fullfile(handles.mother_window_handles.directoryname, handles.mother_window_handles.pooled_downsampled_filename);
if exist(pooled_downsampled_filename)==2
    delete(pooled_downsampled_filename);
end
cluster_mst_upsample_filename = fullfile(handles.mother_window_handles.directoryname, handles.mother_window_handles.cluster_mst_upsample_filename);
if exist(cluster_mst_upsample_filename)==2
    delete(cluster_mst_upsample_filename);
end



% parameter update induced change of files. 
function delete_PooledData_Result_mat_files(handles)
pooled_downsampled_filename = fullfile(handles.mother_window_handles.directoryname, handles.mother_window_handles.pooled_downsampled_filename);
if exist(pooled_downsampled_filename)==2
    delete(pooled_downsampled_filename);
end
cluster_mst_upsample_filename = fullfile(handles.mother_window_handles.directoryname, handles.mother_window_handles.cluster_mst_upsample_filename);
if exist(cluster_mst_upsample_filename)==2
    delete(cluster_mst_upsample_filename);
end


% parameter update induced change of files. 
function delete_Result_mat_files(handles)
cluster_mst_upsample_filename = fullfile(handles.mother_window_handles.directoryname, handles.mother_window_handles.cluster_mst_upsample_filename);
if exist(cluster_mst_upsample_filename)==2
    delete(cluster_mst_upsample_filename);
end


% --- Executes on button press in radiobutton_4_kmeans.
function radiobutton_4_kmeans_Callback(hObject, eventdata, handles)
if isequal(handles.mother_window_handles.clustering_algorithm,'kmeans')
    set(handles.radiobutton_4_kmeans,'value',1);
    set(handles.radiobutton_5_agglomerative,'value',0);
    return
end
set(handles.radiobutton_4_kmeans,'value',1);
set(handles.radiobutton_5_agglomerative,'value',0);
handles.mother_window_handles.clustering_algorithm='kmeans';
guidata(hObject,handles); 
guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
% update the parameter file 
parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
clustering_algorithm = handles.mother_window_handles.clustering_algorithm;
save(parameter_filename,'clustering_algorithm','-append');
% updae .mat files
delete_Result_mat_files(handles)


% --- Executes on button press in radiobutton_5_agglomerative.
function radiobutton_5_agglomerative_Callback(hObject, eventdata, handles)
if isequal(handles.mother_window_handles.clustering_algorithm,'agglomerative')
    set(handles.radiobutton_4_kmeans,'value',0);
    set(handles.radiobutton_5_agglomerative,'value',1);
    return
end
set(handles.radiobutton_4_kmeans,'value',0);
set(handles.radiobutton_5_agglomerative,'value',1);
handles.mother_window_handles.clustering_algorithm='agglomerative';
guidata(hObject,handles); 
guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
% update the parameter file 
parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
clustering_algorithm = handles.mother_window_handles.clustering_algorithm;
save(parameter_filename,'clustering_algorithm','-append');
% updae .mat files
delete_Result_mat_files(handles)


% --- Executes when selected object is changed in uipanel8.
function uipanel8_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel8 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if get(handles.radiobutton_compensation_ignore,'value')==1 && get(handles.radiobutton_compensation_apply,'value')==0
    apply_compensation = 0; % no compensation is applied
elseif get(handles.radiobutton_compensation_ignore,'value')==0 && get(handles.radiobutton_compensation_apply,'value')==1
    apply_compensation = 1; % use compensated data
else
    error('these two radiobuttons cannot be of the same status')
end
if handles.mother_window_handles.apply_compensation==apply_compensation  % if the selection is the same as that in the mother window, do nothing
    return
end
handles.mother_window_handles.apply_compensation = apply_compensation;
guidata(hObject,handles); 
guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
% update the parameter file 
parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
apply_compensation = handles.mother_window_handles.apply_compensation;
save(parameter_filename,'apply_compensation','-append');
% update the .mat files
delete_LD_PooledData_Result_mat_files(handles)
if exist(fullfile(handles.mother_window_handles.directoryname,'check_loaded_data'),'dir')
    delete(fullfile(handles.mother_window_handles.directoryname,'check_loaded_data','*.*'));
    rmdir(fullfile(handles.mother_window_handles.directoryname,'check_loaded_data'));
end