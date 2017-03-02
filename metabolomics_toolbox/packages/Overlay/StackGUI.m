 function varargout = StackGUI(varargin)
%% STACKGUI MATLAB code for StackGUI.fig
%      STACKGUI, by itself, creates a new STACKGUI or raises the existing
%      singleton*.
%
%      H = STACKGUI returns the handle to a new STACKGUI or the handle to
%      the existing singleton*.
%
%      STACKGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STACKGUI.M with the given input arguments.
%
%      STACKGUI('Property','Value',...) creates a new STACKGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before StackGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to StackGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%       
%       Go to Help > Documentation inside the STACKGUI program to see more information
%       on operating the program
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Author: Dawit Woldegiorgis
% Edison Lab, University of Florida

% Last Modified 05-Mar-2015 
% Graphics through GUIDE v2.5

%% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @StackGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @StackGUI_OutputFcn, ...
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


%% --- Executes just before StackGUI is made visible.
function StackGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)
% varargin   command line arguments to StackGUI (see VARARGIN)

% Instantiate handle variables
handles.base = [];          % base spectrum data
handles.overlay = [];       % overlay spectrum data
handles.paths = cell(1);    % location of .mat files (base = 1, overlay = 2:end)
handles.ptr = 0;            % index of item in queue (for overlay)
handles.tptr = 1;           % index of the next empty slot in table
handles.bthresh = -1;       % threshold of base spectra
handles.othresh = -1;       % threshold of overlay spectra (% of max peak)
handles.warning = 1;        % swtich to control warning (warning: on/off)
handles.baseplot = [];      % handle to contour plot of base spectra
handles.overlayplot = [];   % handle to contour plot of overlay spectra
handles.bpeakplot = [];     % handle to scatter plot of base peaks
handles.opeakplot = [];     % handle to scatter plot of overlay peaks
handles.start = 0;          % switch for starting process
handles.peaks = struct('bx',[],'by',[],'ox',[],'oy',[],'thresh',0.3);    % peak data

% Default values
handles.warningIcon = get(handles.disableWarnings,'CData');
handles.peakIcon = get(handles.peakFill,'CData');

% Choose default command line output for StackGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes StackGUI wait for user response (see UIRESUME)
% uiwait(handles.StackGUI);

%% --- Outputs from this function are returned to the command line.
function varargout = StackGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%% **************************** MENU **************************************
%% Menu Item -- File >
function file_Callback(hObject, eventdata, handles)
% hObject    handle to file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)

%% Menu Item -- File > Load >
function file_load_Callback(hObject, eventdata, handles)
% hObject    handle to file_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)

%% Menu Item -- File > Load > Base
function file_load_base_Callback(hObject, eventdata, handles)
% hObject    handle to file_load_base (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)
loadBase_ClickedCallback(hObject, eventdata, handles);

%% Menu Item -- File > Load > Overlay
function file_load_overlay_Callback(hObject, eventdata, handles)
% hObject    handle to file_load_overlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)
loadOverlay_ClickedCallback(hObject, eventdata, handles);

%% Menu Item -- File > Save
function file_save_Callback(hObject, eventdata, handles)
% hObject    handle to file_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)
save_ClickedCallback(hObject, eventdata, handles);

%% Menu Item -- File > Quit
function quit_Callback(hObject, eventdata, handles)
% hObject    handle to quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
quit = questdlg('Any unsaved data will be lost. Are you sure you want to quit?','Quit StackGUI','Yes','No','No');
if strcmp(quit,'Yes')
    close(handles.StackGUI);
else
    return;
end

%% Menu Item -- Help >
function help_Callback(hObject, eventdata, handles)
% hObject    handle to help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)

%% Menu Item -- Help > About
function help_about_Callback(hObject, eventdata, handles)
% hObject    handle to help_about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)
line1 = {'      STACKGUI (version 2.1.1)     '};
line2 = {'                                 '};
line3 = {'  Edison lab, University of Florida'};
line4 = {'            Copyright © 2015        '};
msg = [line1;line2;line3;line4;line2];
msgbox(msg,'About');

%% Menu Item -- Help > Documentation
function help_documentation_Callback(hObject, eventdata, handles)
% hObject    handle to help_documentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)


%% **************************** TOOLBAR ***********************************
%% Toolbar: Load base
function loadBase_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to loadBase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)

% Load sequence: code taken from loadallft.m
[file,path] = uigetfile('*.mat','Select Base Spectrum ...');

% break if user cancels load
if (~file) 
    return; 
end

% Update paths and base name
handles.paths(1) = {strcat(path,file)};
name = strsplit(file,'.');
set(handles.baseName,'String',name{1});

% Get base spectrm
handles.base = load(handles.paths{1});

handles.ptr = 1;                           % update pointer
updateSpectra(hObject, handles,'base');    % update spectra

%% ToolBar: Load Overlay
function loadOverlay_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to loadOverlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)

% Load sequence: code borrowed from loadallft.m
[files,path] = uigetfile('*.mat','Select Overlay Spectra','MultiSelect','on');

if (~iscell(files)) % single or no file selected
   if (~files)   % break if cancelled
       return;
   elseif (ischar(files))  % convert to cell
       files = {files};
   end
end

% Update paths of spectra
temp = cellfun(@(x) ([path x]),files,'UniformOutput',0);
if (isempty(handles.paths{1}) && handles.warning)
    msgbox('Please load the base spectrum before starting.','Base spectrum not loaded','warn');
end
handles.paths = [handles.paths,temp];

% Update names of metabolites in GUI -- extension of filename is removed
names = cellfun(@(x) x(1:end-4),files,'UniformOutput',0);
oldList = cellstr(get(handles.list,'String'));
if isempty(oldList{1}) % nothing loaded -- add as new
    newList = names';
else    % update list
    newList = [oldList;names'];
end
set(handles.list,'String',newList);

% Update handles structure
guidata(hObject, handles);

%% Toolbar: Save Table
function save_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)

% Get table from GUI table
columnName = cellstr(get(handles.table,'ColumnName'))';
data = get(handles.table,'Data');
data = [columnName;data];   % concatenate with columnName

% Prompt save name and location
[name, path] = uiputfile('.csv','Save as','Untitled');
path = [path name];

% Open and override any existing table
fid = fopen(path, 'w');   
for row = 1:size(data,1)
    % stop when encountering empty cells
    if isempty(data{row,1})
        break;
    end
    % write comma delimited csv file
    for col = 1:size(data,2)
        fprintf(fid,'%s,',data{row,col});
    end
    fprintf(fid,'\n');
end
% close file
fclose(fid);

%% Toolbar: Peak Mode > ON
function peakMode_OnCallback(hObject, eventdata, handles)
% hObject    handle to peakMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Warn if not started
if ~handles.start
    if handles.warning
        msgbox('Please START before continuing');
    end
    set(hObject,'State','off');
    return;
end

% Show peak plots of both spectra
set(handles.bpeakplot,'Visible','on');
set(handles.opeakplot,'Visible','on');

% Toggle on showBase and showOverlay
set(handles.showBase,'Value',1);
set(handles.showOverlay,'Value',1);

% Hide contour plots of both spectra
set(handles.baseplot,'Visible','off');
set(handles.overlayplot,'Visible','off');

%% Toolbar: Peak Mode > OFF
function peakMode_OffCallback(hObject, eventdata, handles)
% hObject    handle to peakMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Do nothing if not started
if ~handles.start
    return;
end

% Hide peak plots of both spectra
set(handles.bpeakplot,'Visible','off');
set(handles.opeakplot,'Visible','off');

% Toggle off showBase and showOverlay
set(handles.showBase,'Value',0);
set(handles.showOverlay,'Value',0);

% Show contour plots of both spectra
set(handles.baseplot,'Visible','on');
set(handles.overlayplot,'Visible','on');

%% Toolbar: Fill Peak > ON
function peakFill_OnCallback(hObject, eventdata, handles)
% hObject    handle to peakFill (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Do nothing if no peaks plotted
if isempty(handles.bpeakplot) && isempty(handles.opeakplot)
    set(hObject,'State','off');
    return;
end

% Filled Icon
filledIcon = handles.peakIcon;
fillSpace = filledIcon == 0.94;
filledIcon(fillSpace) = 0;
set(hObject,'CData',filledIcon);

% Set peak plot data point markers to filled
set(handles.bpeakplot,'MarkerFaceColor','flat');
set(handles.opeakplot,'MarkerFaceColor','flat');

% Update handles structure
guidata(hObject,handles);

%% Toolbar: Fill Peak > OFF
function peakFill_OffCallback(hObject, eventdata, handles)
% hObject    handle to peakFill (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Do nothing if no peaks plotted
if isempty(handles.bpeakplot) && isempty(handles.opeakplot)
    return;
end

% Restore original icon
set(hObject,'CData',handles.peakIcon);

% Set peak plot data point markers to unfilled
set(handles.bpeakplot,'MarkerFaceColor','none');
set(handles.opeakplot,'MarkerFaceColor','none');

% Update handles structure
guidata(hObject,handles);

%% Toolbar: Disable warnings > ON
function disableWarnings_OnCallback(hObject, eventdata, handles)
% hObject    handle to disableWarnings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% turn warnings off
handles.warning = 0;    
% Grey out icon
disableIcon = handles.warningIcon;        % RGB color data
disableIcon(:,:,1) = disableIcon(:,:,2);  % R = G
disableIcon(:,:,3) = disableIcon(:,:,2);  % B = G
set(hObject,'CData',disableIcon);

% Update handles structure
guidata(hObject, handles);

%% Toolbar: Disable Warnings > OFF
function disableWarnings_OffCallback(hObject, eventdata, handles)
% hObject    handle to disableWarnings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% turn warnings on
handles.warning = 1;    
% Restore default icon
set(hObject,'CData',handles.warningIcon);

% Update handles structure
guidata(hObject, handles);

%% **************************** MATCH PANEL *******************************
%% --- Executes on button press in yes.
function yes_Callback(hObject, eventdata, handles)
% hObject    handle to yes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)
addToTable(hObject, handles, 'Yes');
handles = guidata(hObject);     % load updated handles
next_Callback(hObject, eventdata, handles);

%% --- Executes on button press in maybe.
function maybe_Callback(hObject, eventdata, handles)
% hObject    handle to maybe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)
addToTable(hObject, handles, 'Maybe');
handles = guidata(hObject);     % load updated handles
next_Callback(hObject, eventdata, handles);

%% --- Executes on button press in no.
function no_Callback(hObject, eventdata, handles)
% hObject    handle to no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)
addToTable(hObject,handles,'No');
handles = guidata(hObject);     % load updated handles
next_Callback(hObject, eventdata, handles);


%% **************************** PEAKS PANEL *******************************
%% --- Executes during object creation, after setting all properties.
function peakThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to peakThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% --- Executres when peakThresh value changes
function peakThresh_Callback(hObject, eventdata, handles)
% hObject    handle to peakThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of peakThresh as text
%        str2double(get(hObject,'String')) returns contents of peakThresh as a double

newThresh = str2double(get(hObject,'String'));
handles.peaks.thresh = newThresh;

% Update peak ratio
if ~isempty(handles.bpeakplot) && ~isempty(handles.opeakplot)
    ratio = getPeakRatio(handles.peaks.bx,handles.peaks.by,handles.peaks.ox,handles.peaks.oy,newThresh);
    set(handles.peakRatio,'String',ratio);
end

% Update handles structure
guidata(hObject,handles);

%% --- Executes on button press in showBase.
function showBase_Callback(hObject, eventdata, handles)
% hObject    handle to showBase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showBase

% Warn if not started
if isempty(handles.baseplot)
    if handles.warning
        msgbox('No base spectrum found. Please load base before continuing');
    end
    set(hObject,'Value',0);     % Toggle off
    return;
end

% Get button toggle state
state = get(hObject,'Value');

% Show/Hide
if state == 1
    set(handles.bpeakplot,'Visible','on');
else
    set(handles.bpeakplot,'Visible','off');
end

%% --- Executes on button press in bpu.
function bpu_Callback(hObject, eventdata, handles)
% hObject    handle to bpu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check for base spectrum
if isempty(handles.bpeakplot)
    % Warn if no base loaded
    if handles.warning
        msgbox('No base spectrum found. Please load base before continuing');
    end
    return;
end

show = get(handles.bpeakplot,'Visible');

% warn if peaks are not shown in GUI
if strcmp(show,'Off') && handles.warning
    msgbox('Warning: You cannnot see the adjustments you are making. Please enable show base peaks to see changes');
end

size = get(handles.bpeakplot,'SizeData');   % get current size
set(handles.bpeakplot,'SizeData',size+2);   % increment size 

%% --- Executes on button press in bpd.
function bpd_Callback(hObject, eventdata, handles)
% hObject    handle to bpd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check for base spectrum
if isempty(handles.bpeakplot)
    if handles.warning
        msgbox('No base spectrum found. Please load base before continuing');
    end
    return;
end

show = get(handles.bpeakplot,'Visible');

% warn if peaks are not shown in GUI
if strcmp(show,'Off') && handles.warning
    msgbox('Warning: You cannnot see the adjustments you are making. Please enable show base peaks to see changes');
end

size = get(handles.bpeakplot,'SizeData');       % get current size
lowPeak = min(size);                            % get smallest peak size:
                                                % useful if size is a vector    
if lowPeak > 2
    set(handles.bpeakplot,'SizeData',size-2);   % decrement size
end


%% --- Executes on button press in showOverlay.
function showOverlay_Callback(hObject, eventdata, handles)
% hObject    handle to showOverlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showOverlay

% Check for overlay spectrum
if isempty(handles.overlayplot)
    % Warn if not started
    if handles.warning
        msgbox('No overlay sepectrum found. Please START before continuing');
    end
    set(hObject,'Value',0);     % Toggle off
    return;
end

% Get button toggle state
state = get(hObject,'Value');

% Show/Hide
if state == 1
    set(handles.opeakplot,'Visible','on');
else
    set(handles.opeakplot,'Visible','off');
end

%% --- Executes on button press in opu.
function opu_Callback(hObject, eventdata, handles)
% hObject    handle to opu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check for overlay spectrum
if isempty(handles.opeakplot)
    if handles.warning
        msgbox('No overlay spectrum found. Please START before continuing');
    end
    return;
end

show = get(handles.opeakplot,'Visible');

% warn if peaks are not shown in GUI
if strcmp(show,'Off') && handles.warning
    msgbox('Warning: You cannnot see the adjustments you are making. Please enable show overlay peaks to see changes');
end

size = get(handles.opeakplot,'SizeData');   % get current size
set(handles.opeakplot,'SizeData',size+2);   % increment size 

%% --- Executes on button press in opd.
function opd_Callback(hObject, eventdata, handles)
% hObject    handle to opd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check for overlay spectrum
if isempty(handles.opeakplot)
    if handles.warning
        msgbox('No overlay spectrum found. Please START before continuing');
    end
    return;
end

show = get(handles.opeakplot,'Visible');

% warn if peaks are not shown in GUI
if strcmp(show,'Off') && handles.warning
    msgbox('Warning: You cannnot see the adjustments you are making. Please enable show overlay peaks to see changes');
end

size = get(handles.opeakplot,'SizeData');       % get current size
lowPeak = min(size);                            % get smallest peak size:
                                                % useful if size is a vector    
if lowPeak > 2
    set(handles.opeakplot,'SizeData',size-2);   % decrement size
end

%% **************************** CONTROLS PANEL ****************************
%% --- Executes on button press in back
function back_Callback(hObject, eventdata, handles)
% hObject    handle to back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)

% start
if (~handles.start)
    start_Callback(hObject, eventdata, handles);
    return;
end

% step back in queue
if (handles.ptr > 2)    
    handles.ptr = handles.ptr-1;
    if get(handles.list,'Value') > 1
        newValue = get(handles.list,'Value') - 1;
        set(handles.list,'Value',newValue);
    end
else
    if (handles.warning)
        msgbox('You''re at the top of the list.');
    end
    return;
end

handles.othresh = -1;                      % reset threshold
updateSpectra(hObject,handles,'overlay');  % update spectra

%% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)

% start
if (~handles.start)
    start_Callback(hObject, eventdata, handles);
    return;
end

listValue = get(handles.list,'Value');
listItems = cellstr(get(handles.list,'String'));

% step forward in queue
if (handles.ptr < length(handles.paths))
    if (listValue < length(listItems)) && (handles.ptr > 1)
        set(handles.list,'Value',listValue+1);
    end
    handles.ptr = handles.ptr+1;   
else
    if (handles.warning)
        msgbox('Check complete! You''re at the end of the list.');
    end
    return;
end

handles.othresh = -1;                      % reset threshold
updateSpectra(hObject,handles,'overlay');  % update spectra


%% --- Executes on button press in start.
function start_Callback(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Update switch
if ~handles.start
    handles.start = 1;
else
    return;
end

% check if base is loaded
if (isempty(handles.paths{1}))
    if (handles.warning)
        msgbox('You need to load the base spectrum before starting','No base spectrum','warn');
    end
    return;
elseif length(handles.paths) == 1
    if (handles.warning)
        msgbox('You need to load the overlay spectra before starting', 'No overlay spectra','warn');
    end
    return;
end

% plot first overlay spectrum
if handles.ptr == 1
    handles.ptr = handles.ptr+1;
    updateSpectra(hObject,handles,'overlay');
end

%% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)

% Do nothing if base is not loaded
if isempty(handles.paths{1})
    return
end

reply = questdlg('Unsaved data will be lost. Are you sure you want to RESET?','Reset program','Yes','No','No');

if strcmp(reply,'No') || strcmp(reply,'Cancel')
    return;
end

% Delete plots
delete(handles.baseplot);               % delete base plot
delete(handles.overlayplot);            % delete overlay plot
delete(handles.bpeakplot);              % delete base peak plot
delete(handles.opeakplot);              % delete overlay peak plot

% Reset global variables
handles.ptr = 1;                        % reset pointer
handles.tptr = 1;                       % reset table pointer
handles.paths = cell(1);                % empty paths list
handles.bthresh = -1;                   % reset base thresh
handles.othresh = -1;                   % reset overlay thresh
handles.baseplot = [];                  % reset to empty
handles.overlayplot = [];               % reset to empty
handles.opeakplot = [];                 % reset to empty
handles.bpeakplot = [];                 % reset to emtpy
handles.start = 0;                      % reset start flag
handles.peaks.bx = [];                  % reset peak data
handles.peaks.by = [];
handles.peaks.ox = [];                  
handles.peaks.oy = [];
handles.peaks.thresh = 0.3;

% Reset GUI variables
set(handles.list,'Value',1);            % reset list pointer
set(handles.table,'Data',cell(1,3));    % empty table data
set(handles.list,'String',[]);          % emtpy list contents
set(handles.peakRatio,'String',[]);     % remove peak ratio
set(handles.baseName,'String',[]);      % remove base name
set(handles.peakThresh,'String',0.3);   % reset peak thresh

% Update handles structure
guidata(hObject,handles);

%% **************************** CONTROLS: BASE ****************************
%% --- Executes during object creation, after setting all properties.
function bThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% --- Threshold input for base
function bThresh_Callback(hObject, eventdata, handles)
% hObject    handle to bthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bthresh as text
%        str2double(get(hObject,'String')) returns contents of bthresh as a double
handles.bthresh = str2double(get(hObject,'String'));
if handles.bthresh < 0 || isnan(handles.bthresh)
    handles.bthresh = 0;
end

updateSpectra(hObject,handles,'base');


%% --- Executes on button press in base thresh up.
function btu_Callback(hObject, eventdata, handles)
% hObject    handle to btu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)

handles.bthresh = handles.bthresh + 10;
updateSpectra(hObject,handles,'base');

%% --- Executes on button press in btd.
function btd_Callback(hObject, eventdata, handles)
% hObject    handle to btd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)

if (handles.bthresh < 10)
    handles.bthresh = 0;
else
    handles.bthresh = handles.bthresh - 10;
end

updateSpectra(hObject,handles,'base');

%% **************************** CONTROLS: OVERLAY *************************
%% --- Executes during object creation, after setting all properties.
function oThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to othresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% --- Threshold input for overlay
function oThresh_Callback(hObject, eventdata, handles)
% hObject    handle to oThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)

% Hints: get(hObject,'String') returns contents of oThresh as text
%        str2double(get(hObject,'String')) returns contents of oThresh as a double
handles.othresh = str2double(get(hObject,'String'));
if handles.othresh < 0 || isnan(handles.othresh)
    handles.othresh = 0;
end

updateSpectra(hObject,handles,'overlay');

%% --- Executes on button press in overlay thresh up.
function otu_Callback(hObject, eventdata, handles)
% hObject    handle to otu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)

handles.othresh = handles.othresh + 10;
updateSpectra(hObject,handles,'overlay');

%% --- Executes on button press in overlay thresh down
function otd_Callback(hObject, eventdata, handles)
% hObject    handle to otd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)

if handles.othresh < 10
    handles.othresh = 0;
else
    handles.othresh = handles.othresh - 10;
end

updateSpectra(hObject,handles,'overlay');

%% **************************** METABOLITES PANEL ************************* 
%% --- Executes during object creation, after setting all properties.
function list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% --- Executes on selection change in list.
function list_Callback(hObject, eventdata, handles)
% hObject    handle to list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user table (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list

if (handles.ptr == get(hObject,'Value')+1 || isempty(get(hObject,'String')))  
    return;     % do nothing
else
    handles.start = 1;
    handles.ptr = get(hObject,'Value')+1;   % jump to selection
end
updateSpectra(hObject,handles,'overlay');

%% --- Executes on key press with focus on list and none of its controls.
function list_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to list (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

% BUGGY

listItems = cellstr(get(hObject,'String'));
listValue = get(hObject,'Value');
        
if strcmp(eventdata.Key,'delete') || strcmp(eventdata.Key,'d')
    if isempty(listItems{1})
        return;
    else
        handles.ptr = listValue+1;
        handles.paths(handles.ptr) = [];            % delete from path list
        listItems(listValue) = [];                  % delete from list
        set(handles.list,'String',listItems);       % update list
        if(listValue < 1)
            listValue = listValue - 1;
            set(handles.list,'Value',listValue);    % update value
        end
    end
end

%% **************************** HELPER FUCNTIONS ************************** 
%% --- Add data to table
function addToTable(hObject,handles,call)
% hObject   handle to caller
% handles   strucutre with handles and user table (see GUIDATA)
% exists    string to determine how item is added to table

if handles.ptr <= 1
    if handles.warning
        msgbox('There is nothing plotted yet. Please press START before trying to add to table.','Nothing plotted','help');
    end
    return;     % stop if nothing is plotted
end

Data = get(handles.table,'Data');                   % table data
listItems = cellstr(get(handles.list,'String'));    % list data
listValue = get(handles.list,'Value');              % current list item
name = listItems{listValue};
peakRatio = get(handles.peakRatio,'String');        % current peak ratio

if strcmp(call,'No') addNew = false;     % flag: update current data
else                 addNew = true;      % flag: add new data
end

% Check for duplicate and update
for i = 1:size(Data,1)
    if strcmp(Data(i,1),name)
        if strcmp(call,'No')
            Data(i,:) = [];                     % delete row
            handles.tptr = handles.tptr - 1;    % update pointer
            break;                              % exit loop
        end
        addNew = false;           % flag: update current data
        Data(i,2) = {call};       % update existence
        Data(i,3) = {peakRatio};  % update peak ratio
    end
end

% Insert new data
if addNew
    Data(handles.tptr,1) = {name};       % insert name
    Data(handles.tptr,2) = {call};       % insert existence
    Data(handles.tptr,3) = {peakRatio};  % insert peak ratio
    handles.tptr = handles.tptr + 1;     % update pointer
end

set(handles.table,'Data',Data);   % update table data

% Update handles structure
guidata(hObject,handles);


%% --- Update Spectra
function updateSpectra(hObject,handles,type)
% hObject   handle to caller
% handles   strucutre with handles and user table (see GUIDATA)
% type      determines which type of spectrum to update

% Determine plot visibility
showContour = get(handles.peakMode,'State');
if strcmp(showContour,'off')
    showContour = 1;
else
    showContour = 0;
end
showBasePeaks = get(handles.showBase,'Value');
showOverlayPeaks = get(handles.showOverlay,'Value');

switch (type)
    case 'base'
        % Plot contour
        delete(handles.baseplot);  % delete old base spectrum
        if handles.bthresh == -1   % plot with default threshold
            [handles.baseplot,handles.bthresh] = contourplot(handles.base.X,handles.base.XNoise,handles.base.ppm1,handles.base.ppm2,'blue');
        else
            handles.baseplot = contourplot(handles.base.X,handles.base.XNoise,handles.base.ppm1,handles.base.ppm2,'blue',handles.bthresh);
        end
        % Set contour visibility
        if ~showContour
            set(handles.baseplot,'Visible','off');
        end
        % Plot peak
        delete(handles.bpeakplot);  % delete old base peak plot
        [handles.bpeakplot,handles.peaks.bx,handles.peaks.by] = peakplot(handles.base,handles.bthresh,'cyan');
        % Set peak visiblity
        if ~showBasePeaks
            set(handles.bpeakplot,'Visible','off');
        end
        % update base threshold (GUI)
        set( handles.bThresh,'String',handles.bthresh);  
    case 'overlay'
        % Plot contour
        delete(handles.overlayplot);  % delete old overlay spectrum
        handles.overlay = load(handles.paths{handles.ptr}); % get overlay spectra
        if handles.othresh == -1
            % plot with default threshold
            [handles.overlayplot,handles.othresh] = contourplot(handles.overlay.X,handles.overlay.XNoise,handles.overlay.ppm1,handles.overlay.ppm2,'red');
        else
            % plot with specified threshold
            handles.overlayplot = contourplot(handles.overlay.X,handles.overlay.XNoise,handles.overlay.ppm1,handles.overlay.ppm2,'red',handles.othresh);
        end
        % Set contour visibility
        if ~showContour
            set(handles.overlayplot,'Visible','off');
        end
        % plot overlay peaks
        delete(handles.opeakplot);  % delete old overlay peak plot
        [handles.opeakplot,handles.peaks.ox,handles.peaks.oy] = peakplot(handles.overlay,handles.othresh,'magenta');
        % Set peak visibility
        if ~showOverlayPeaks
            set(handles.opeakplot,'Visible','off');
        end
        set(handles.oThresh,'String',handles.othresh);  % update overlay threshold (GUI)
end

% Check peak fill state
if strcmp(get(handles.peakFill,'State'),'on')
   peakFill_OnCallback(handles.peakFill,[],handles) 
end

% Set peak ratio
if ~isempty(handles.bpeakplot) && ~isempty(handles.opeakplot)
    ratio = getPeakRatio(handles.peaks.bx,handles.peaks.by,handles.peaks.ox,handles.peaks.oy,handles.peaks.thresh);
    set(handles.peakRatio,'String',ratio);
end

% Update handles structure
guidata(hObject, handles);

%% Calculate and set the peak ratio
function  [ratio] = getPeakRatio(bx,by,ox,oy,thresh)
% bx        x-Coordinates of peaks in base spectrum
% by        y-Coordinates of peaks in base spectrum
% ox        x-Coordinates of peaks in overlay spectrum
% oy        y-Coordinates of peaks in overlay spectrum
% thresh    maximum distance between peaks to consider them overlapping

overlap = 0;
% Determine overlap: Overlap occurs when the distance between the overlay
% peak in question and any other peak in base is less than the threshold
% specified

for i = 1:length(ox) % get peak from overlay
    for j = 1:length(bx) % compare to each peak in base
        dx = ox(i) - bx(j);
        dy = oy(i) - by(j);
        peakDistance = sqrt(dx^2 + dy^2);   % distance b/n peaks
        if peakDistance < thresh
            overlap = overlap + 1; % overlapping peaks found
            break;
        end
    end
end     

% Calculate ratio (overlapping peaks:total peaks in overlay)
ratio = overlap/length(ox);

%% --- Peak Plot
function [handle,x,y] = peakplot(spectrum,thresh,color)
% spectrum      struct containing data of spectrum
%               Fields:  X,XNoise,ppm1,ppm2,XTitles
%               **See CONTOURPLOT for details on fields
% thresh        threshold
% handle        handle to the scatter plot of peaks

data = abs(spectrum.X ./ spectrum.XNoise);    % positive peaks
[row,col] = getPeaks(data,thresh);
x = spectrum.ppm1(col);     % x-Coordinates of peaks
y = spectrum.ppm2(row);     % y-Coordinates of peaks

% Determine peak sizes
peaks = zeros(1,length(col));     % instantiate peaks
for i = 1:length(col)
    peaks(i) = data(row(i),col(i));     % get peak heights
end
maxPeak = max(peaks);           % get max peak height
size = (peaks./maxPeak)*36;     % size vector based on relative height
% plot peaks
hold on;
handle = scatter(x,y,size,color);


%% --- Get Peaks: Find row and column indices of peaks
function [row,col] = getPeaks(data,thresh)
% data      matrix containing data points
% thresh    threshold to determine minmum peak height
% row       row index of peaks in data
% col       column index of peaks in data

if ~exist('thresh','var')
    thresh = mean(data(data~=0))+5*std(data(data~=0));  % default thresh
    thresh = thresh * 0.125;    % minimum peak height
else
    thresh = thresh * 0.1250;   % minimum peak height
end

% Truncate data less than minimum peak height to zero
truncate = data < thresh;
data(truncate) = 0;          

% Add a "boundary" of zeros:
% Doing this will allow use of single algorithm for all data points.
% Note that only three data points need to be compared against for edges
% and two for corners of matrix. Furthermore, the direction of comparison
% changes for data points on different edges of the matrix. This procedure
% eliminates accounting for the variations.
zrow = zeros(1,size(data,2));
data = [zrow;data;zrow];
zcol = zeros(size(data,1),1);
data = [zcol data zcol];

% Iterate through data
nrows = size(data,1);    % number of rows
ncols = size(data,2);    % number of cols
row = [];
col = [];
% Note: first and last indices are "boundary" zeros. See above.
for i = 2:nrows-1
    for j = 2:ncols-1
        if data(i,j) > data(i,j-1) && data(i,j) > data(i,j+1)
            if data(i,j) > data(i+1,j) && data(i,j) > data(i-1,j)
                % data (i,j) is local maxima
                % Note that (i,j) == (i-1,j-1) in actual data, due to the
                % boundary of zeros
                row = [row i-1];  % add to list of row indices
                col = [col j-1];  % add to list of column indices
                j = j+1;    % skip next data point since this > next
            end
        end
    end
end

%% --- Contourplot: Based on stackplot.m with some modifications
function [handle,thresh] = contourplot(X,XNoise,ppm1,ppm2,color,thresh)
 % Plots stacked contour plot of X colored by Y (if known).
% 
% X            table matrix of spectra
% XNoise	   Calculated noise matrix
% ppm1         Chemical shift vector of F2
% ppm2         Chemical shift vector of F1
% thresh       Intensity threshold- can be scalar or vector; if vector
%              same number of elements as size(X,3)
% handle       handle to the contour plot figure
    
data=abs(X./XNoise);        % data matrix with noise reduction

if ~exist('thresh','var')
    thresh=mean(data(data~=0))+5*std(data(data~=0));    % default threshold
end

range = [-3:2/3:3];         % log range of contour levels
vector=(2.^range)*thresh;   % contour levels
hold on;
[~,handle] = contour(ppm1,ppm2,data,vector,'EdgeColor',color);
set(gca,'XDir','rev');      % reverse X-axis
set(gca,'YDir','rev');      % reverse Y-axis
