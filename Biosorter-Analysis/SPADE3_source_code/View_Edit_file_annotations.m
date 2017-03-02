function varargout = View_Edit_file_annotations(varargin)
% VIEW_EDIT_FILE_ANNOTATIONS MATLAB code for View_Edit_file_annotations.fig
%      VIEW_EDIT_FILE_ANNOTATIONS, by itself, creates a new VIEW_EDIT_FILE_ANNOTATIONS or raises the existing
%      singleton*.
%
%      H = VIEW_EDIT_FILE_ANNOTATIONS returns the handle to a new VIEW_EDIT_FILE_ANNOTATIONS or the handle to
%      the existing singleton*.
%
%      VIEW_EDIT_FILE_ANNOTATIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEW_EDIT_FILE_ANNOTATIONS.M with the given input arguments.
%
%      VIEW_EDIT_FILE_ANNOTATIONS('Property','Value',...) creates a new VIEW_EDIT_FILE_ANNOTATIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before View_Edit_file_annotations_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to View_Edit_file_annotations_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help View_Edit_file_annotations

% Last Modified by GUIDE v2.5 18-May-2011 10:43:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @View_Edit_file_annotations_OpeningFcn, ...
                   'gui_OutputFcn',  @View_Edit_file_annotations_OutputFcn, ...
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


% --- Executes just before View_Edit_file_annotations is made visible.
function View_Edit_file_annotations_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to View_Edit_file_annotations (see VARARGIN)

% Choose default command line output for View_Edit_file_annotations
handles.output = hObject;
handles.all_fcs_filenames = varargin{1};
handles.file_annot = varargin{2};
handles.mother_window_handles = varargin{3};
handles.table_data = [handles.all_fcs_filenames,handles.file_annot];
set(handles.uitable,'data',handles.table_data);
set(handles.uitable,'columnEditable',[false true]);
set(handles.uitable,'columnname',[{'filename'},{'file short annotation'}]);
set(handles.uitable,'ColumnWidth',[{330},{100}]);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes View_Edit_file_annotations wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = View_Edit_file_annotations_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_close.
function button_close_Callback(hObject, eventdata, handles)
% hObject    handle to button_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.table_data = get(handles.uitable,'data');
handles.mother_window_handles.file_annot = handles.table_data(:,2);
parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
file_annot = handles.mother_window_handles.file_annot;
save(parameter_filename,'file_annot','-append');
delete(handles.figure1);


% --- Executes when entered data in editable cell(s) in uitable.
function uitable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.table_data = get(handles.uitable,'data');
r = eventdata.Indices(1);
handles.file_annot(r) = handles.table_data(r,2);
handles.mother_window_handles.file_annot(r) = handles.table_data(r,2);
guidata(hObject,handles);
guidata(handles.mother_window_handles.button_browse_directory, handles.mother_window_handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% handles.table_data = get(handles.uitable,'data');
% handles.mother_window_handles.file_annot = handles.table_data(:,2);
% parameter_filename = fullfile(handles.mother_window_handles.directoryname,handles.mother_window_handles.parameter_filename);
% file_annot = handles.mother_window_handles.file_annot;
% save(parameter_filename,'file_annot','-append');
% delete(hObject);
