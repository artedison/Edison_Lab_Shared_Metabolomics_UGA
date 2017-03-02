function varargout = View_Edit_file_annotation(varargin)
% VIEW_EDIT_FILE_ANNOTATION M-file for View_Edit_file_annotation.fig
%      VIEW_EDIT_FILE_ANNOTATION, by itself, creates a new VIEW_EDIT_FILE_ANNOTATION or raises the existing
%      singleton*.
%
%      H = VIEW_EDIT_FILE_ANNOTATION returns the handle to a new VIEW_EDIT_FILE_ANNOTATION or the handle to
%      the existing singleton*.
%
%      VIEW_EDIT_FILE_ANNOTATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEW_EDIT_FILE_ANNOTATION.M with the given input arguments.
%
%      VIEW_EDIT_FILE_ANNOTATION('Property','Value',...) creates a new VIEW_EDIT_FILE_ANNOTATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before View_Edit_file_annotation_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to View_Edit_file_annotation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help View_Edit_file_annotation

% Last Modified by GUIDE v2.5 03-May-2011 16:24:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @View_Edit_file_annotation_OpeningFcn, ...
                   'gui_OutputFcn',  @View_Edit_file_annotation_OutputFcn, ...
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


% --- Executes just before View_Edit_file_annotation is made visible.
function View_Edit_file_annotation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to View_Edit_file_annotation (see VARARGIN)

% Choose default command line output for View_Edit_file_annotation
handles.output = hObject;
handles.all_fcs_filenames = varargin{1};
handles.file_annot = varargin{2};
handles.mother_window_handles = varargin{3};
% Update handles structure
% display filenames and the file_annot in the table below
set(handles.activex1,'Rows',length(handles.all_fcs_filenames)+1, 'Cols',3);
set(handles.activex1,'Row',0,'Col',1,'Text','filename');
set(handles.activex1,'Row',0,'Col',2,'Text','file short annotation');
for i=1:length(handles.all_fcs_filenames)
    set(handles.activex1,'Row',i,'Col',0,'Text',num2str(i));
    set(handles.activex1,'Row',i,'Col',1,'Text',handles.all_fcs_filenames{i});
    set(handles.activex1,'Row',i,'Col',2,'Text',handles.file_annot{i});
end 
set(handles.activex1,'ColWidth',0,500);
set(handles.activex1,'ColWidth',1,6000);
set(handles.activex1,'ColWidth',2,1500);
guidata(hObject, handles);

% UIWAIT makes View_Edit_file_annotation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = View_Edit_file_annotation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.figure1;


% --------------------------------------------------------------------
function activex1_AfterEdit(hObject, eventdata, handles)
% hObject    handle to activex1 (see GCBO)
% eventdata  structure with parameters passed to COM event listener
% handles    structure with handles and user data (see GUIDATA)
r = get(handles.activex1,'Row');
c = get(handles.activex1,'Col');
if c==1
    set(handles.activex1,'Row',r,'Col',c,'Text',handles.all_fcs_filenames{r});
end
if c==2
%     handles.file_annot{r} = get(handles.activex1,'Text');  % somehow that guidata() function does not work within this call back of activeX control%
% NOTE: I have tried to call another of my own function from here, but it does not work either, the same error message appears.
    handles.mother_window_handles.file_annot{r} = get(handles.activex1,'Text');
    guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
end




% --------------------------------------------------------------------
function activex1_KeyUpEdit(hObject, eventdata, handles)
% hObject    handle to activex1 (see GCBO)
% eventdata  structure with parameters passed to COM event listener
% handles    structure with handles and user data (see GUIDATA)
r = get(handles.activex1,'Row');
c = get(handles.activex1,'Col');
if c==1
    set(handles.activex1,'Row',r,'Col',c,'Text',handles.all_fcs_filenames{r});
end




% --- Executes on button press in button_close.
function button_close_Callback(hObject, eventdata, handles)
% hObject    handle to button_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure1_CloseRequestFcn(handles.figure1, [], handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
for i=1:length(handles.all_fcs_filenames)
    set(handles.activex1,'Row',i,'Col',2)
    handles.mother_window_handles.file_annot{i} = get(handles.activex1,'Text');
end 
% update the parameter file
parameter_filename = fullfile(handles.mother_window_handles.directoryname,'SPADE_parameters.mat');
file_annot = handles.mother_window_handles.file_annot;
save(parameter_filename,'file_annot', '-append');
% update this variable in the main window
guidata(handles.mother_window_handles.button_browse_directory,handles.mother_window_handles);
delete(hObject);


