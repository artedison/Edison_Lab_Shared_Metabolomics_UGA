function varargout = MetaNMRQC(varargin)
% METANMRQC MATLAB code for MetaNMRQC.fig
%      METANMRQC, by itself, creates a new METANMRQC or raises the existing
%      singleton*.
%
%      H = METANMRQC returns the handle to a new METANMRQC or the handle to
%      the existing singleton*.
%
%      METANMRQC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in METANMRQC.M with the given input arguments.
%
%      METANMRQC('Property','Value',...) creates a new METANMRQC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MetaNMRQC_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MetaNMRQC_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MetaNMRQC

% Last Modified by GUIDE v2.5 01-May-2014 12:19:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MetaNMRQC_OpeningFcn, ...
                   'gui_OutputFcn',  @MetaNMRQC_OutputFcn, ...
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


% --- Executes just before MetaNMRQC is made visible.
function MetaNMRQC_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MetaNMRQC (see VARARGIN)

% Choose default command line output for MetaNMRQC
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MetaNMRQC wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MetaNMRQC_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadSpectra.
function loadSpectra_Callback(hObject, eventdata, handles)
% hObject    handle to loadSpectra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ppm X y S_id;
loadallft;
[X,ppm,Xtitle]=Setup1D(spectra);
y=dlmread('yVec.txt');
S_id=char(Xtitle);
handles.plotF=display1D(X,ppm,y);
guidata(hObject,handles);

% --- Executes on button press in showID.
function showID_Callback(hObject, eventdata, handles)
% hObject    handle to showID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(0,'CurrentFigure',handles.plotF);
dcm = datacursormode(gcf);
datacursormode on;
set(dcm,'UpdateFcn',@whichLineCallback);
guidata(hObject,handles);

% --- Executes on button press in correctBase.
function correctBase_Callback(hObject, eventdata, handles)
% hObject    handle to correctBase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ppmr Xr Xb y;
disp('It may take a few minutes, Please be patient...');
Xb=CorrectBl(Xr);
handles.plotF=compare1D(Xr,Xb,ppmr,y);
guidata(hObject,handles);

% --- Executes on button press in alignSpectra.
function alignSpectra_Callback(hObject, eventdata, handles)
% hObject    handle to alignSpectra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ppmr Xb Xa y;
Alignment;
guidata(hObject,handles);


% --- Executes on button press in removeReg.
function removeReg_Callback(hObject, eventdata, handles)
% hObject    handle to removeReg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ppm Xf Xr ppmr y;
[Xr1,ppmr]=remove_ends(Xf,ppm,0.01,11.0);
Xr=remove_region(Xr1,ppmr,4.4,4.89);
handles.plotF=display1D(Xr,ppmr,y);
guidata(hObject,handles);


% --- Executes on button press in correctRef.
function correctRef_Callback(hObject, eventdata, handles)
% hObject    handle to correctRef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ppm X y Xf;
 [~,k0]=min(abs(ppm));
 [~,k1]=min(abs(ppm+0.5));
for i=1:size(X,1)
 XDSS=X(i,k1:k0);
 [maxtab,~]=peakdet(XDSS,1e5);
 [~,idx]=max(maxtab(:,2));
 kdss=k0-k1-maxtab(idx,1);
 Xf(i,:)=circshift(X(i,:),[0,kdss]);
end
handles.plotF=compare1D(X,Xf,ppm,y);
guidata(hObject,handles);


function output_txt = whichLineCallback(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

%In a plot with a bunch of lines, figure out which line I'm clicking
global S_id;
persistent oldLH oldColor oldLW;
lineHandle=get(event_obj, 'Target');
oldColor=get(lineHandle,'Color');
oldLW=get(lineHandle,'LineWidth');
if ishandle(oldLH) 
    set(oldLH,'LineWidth',oldLW,'Color',oldColor);
end
parent=get(lineHandle,'parent');
axes=get(parent,'parent');

groupHandles=get(axes,'children');
children=[];
for i=1:length(groupHandles)
    children=[children; get(groupHandles(i),'children')];
end    

idx=find(sort(children)==lineHandle);

output_txt=['Spectrum: ', S_id(idx,1:2)];
set(lineHandle,'LineWidth',5,'Color','g');
oldLH=lineHandle;



% --- Executes on button press in calcP.
function calcP_Callback(hObject, eventdata, handles)
% hObject    handle to calcP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ppmr Xb y;
figure;
[H,Hn,Ib,Is]=opt_bucket(ppmr,Xb,.05,.5);
beta=0.05/size(Ib,1);
pos=find(y==0);
neg=find(y==1);
for k = 1:size(Ib,1)
    [h(k),p(k)]=ttest2(Hn(pos,k),Hn(neg,k),beta);
end
[p,sigs]=MWAS(Hn,y,'bonferroni');
manhattan(Hn,y,Ib(:,1),p,sigs,'bonferroni');


% --- Executes on button press in checkVar.
function checkVar_Callback(hObject, eventdata, handles)
% hObject    handle to checkVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ppmr Xb;
varcheck(Xb, ppmr);
