function varargout = Alignment(varargin)
% ALIGNMENT MATLAB code for Alignment.fig
%      ALIGNMENT, by itself, creates a new ALIGNMENT or raises the existing
%      singleton*.
%
%      H = ALIGNMENT returns the handle to a new ALIGNMENT or the handle to
%      the existing singleton*.
%
%      ALIGNMENT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ALIGNMENT.M with the given input arguments.
%
%      ALIGNMENT('Property','Value',...) creates a new ALIGNMENT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Alignment_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Alignment_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Alignment

% Last Modified by GUIDE v2.5 01-May-2014 11:21:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Alignment_OpeningFcn, ...
                   'gui_OutputFcn',  @Alignment_OutputFcn, ...
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


% --- Executes just before Alignment is made visible.
function Alignment_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Alignment (see VARARGIN)

% Choose default command line output for Alignment
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Alignment wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Alignment_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in guide_method.
function guide_method_Callback(hObject, eventdata, handles)
% hObject    handle to guide_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns guide_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from guide_method


% --- Executes during object creation, after setting all properties.
function guide_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to guide_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in star_method.
function star_method_Callback(hObject, eventdata, handles)
% hObject    handle to star_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns star_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from star_method


% --- Executes during object creation, after setting all properties.
function star_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to star_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in cluster.
function cluster_Callback(hObject, eventdata, handles)
% hObject    handle to cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cluster contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cluster


% --- Executes during object creation, after setting all properties.
function cluster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in represent.
function represent_Callback(hObject, eventdata, handles)
% hObject    handle to represent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns represent contents as cell array
%        contents{get(hObject,'Value')} returns selected item from represent


% --- Executes during object creation, after setting all properties.
function represent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to represent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in align.
function align_Callback(hObject, eventdata, handles)
% hObject    handle to align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ppmr Xr Xa y;
star_methods=cellstr(get(handles.star_method,'String'));
star_methods_selected=star_methods{get(handles.star_method,'Value')};
rep_spectra=cellstr(get(handles.represent,'String'));
rep_spectra_selected=rep_spectra{get(handles.represent,'Value')};
Xa=star_align1D(Xr,ppmr,rep_spectra_selected,star_methods_selected,0.08,0.5);
compare1D(Xr,Xa,ppmr,y);