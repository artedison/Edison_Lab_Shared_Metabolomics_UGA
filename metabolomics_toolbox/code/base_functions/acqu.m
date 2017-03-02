function acqu(location,flag)
%% ACQU gets the data from the acqu file, in a SECIM folder hierarchy, and
% saves it to CSV file in methods.
% Input:
%   location    name/path of parent folder. Default is current folder
%               * Note: give it the raw_data folder, but works with any
%               folder having the same structure format
%   flag        OPTIONAL: Detemine organizaiton of data in workspace
% Output:
%   CSV file in the methods folder, named according to the child folder
%   from which data was extracted.

% instantiate location
if ~exist('location','var')
    location = pwd;
    home = location;
else
    home = pwd;
end

%instantiate flag
if ~exist('flag','var')
    flag = 0;
end

% go to location
try 
    cd(location);
    location = pwd;     % change location to path format
catch
    error(strcat('No "',location,'" directory found here.'));
end

parentPath = fileparts(location);           % get path to parent folder
parentSub = getFolders(parentPath);         % get subfolders in parent
subFolders = getFolders(location);          % get subfolders in target

found = false;               % flag to determine if "methods" folder exists
                                    
% search for index of methods in subFolders
for methods_index = 1:length(parentSub)
    name = parentSub{methods_index};    % get folder name
    % check if name contains 'methods' at the end
    if length(name) >= 7 
        if strcmpi(name(end-6:end),'methods') 
            found = true;
            break;
        end
    end
end

% prompt for new folder creation if no methods folder found
if ~found
    warning(strcat('No "methods" folder found in ',parentPath));
    create = input('Do you want to create new "methods" folder? (Y/N) ','s');
    if strcmpi(create,'Y')
        [~,parentName] = fileparts(parentPath);
        name = strcat(parentName,'_methods');
        mkdir(parentPath,name);   % create folder in parent
    else
        if ~strcmpi(create,'N')
            warning('Unkown command.');
        end
        error('Cannot contiue without required folder. Program will terminate.');
        cd(home);
        return;
    end
end

methodsPath = strcat(parentPath,'/',name);
    
% Extract data from each folder, if available
for index = 1:length(subFolders)
    folder = subFolders{index};             % folder name
    data = getData(folder);                 % extract acqu data
    if ~isempty(data)                       % does data exist in folder?
        save(data,folder,methodsPath,flag);    % save data to methods folder
    end
end

if exist('home','var')
    cd(home);       % return home
end

function [data] = getData(path)
%% GETDATA extracts the variables specified in "var" from the acqu/acqus
% file (1D/2D respectively).

% Modify this var list if you want to add/delete variables you want to get
% from the acqus file.
var = {'D1','DS','DW','DE','EXP','NS','NUC1','O1','PL','PULPROG','RG', ...
       'SFO1','SOLVENT','SW','TD','TD0','TE'};
var = [var,{'P1'}];    % DO NOT MODIFY: unless you want to discard P1!
file = 'acqus';

data = cell(0);              % inintialize data
rows = cell(0);              % initialize row name

home = pwd;                  % save caller's path
cd(path);                    % go to path
subFolders = getFolders;     % children folders
for f = 1:length(subFolders) % get data from all subFolders
    cd(subFolders{f});       % go to each subFolder
    files = getFolders(1);   % get list of files
    for i = 1:length(files)
        if strcmp(files{i},file)
            filePath = strcat(pwd,'/',file);
            data = [data;getValue(filePath,var)];   % concatenate w/ old data
            rows = [rows;subFolders(f)];            % concatenate w/ row names
        end
    end
    cd ..;                   % go back to parent
end
cd(home);                   % go back to caller's path

% Construct table
if ~isempty(data)
    data = [rows,data];     % concatenate row names with respective data value
    header = [{'run_order'},var];      % leave blank header for column with row names
    data = [header;data];   % concatenate header with data
end

function data = getValue(path,var)
% GETVALUE returns the values of the variables in the file specified

% Input:
%   path    location of file in string format
%   var     cell array with list of variables to be searched.
% Output:
%   data    cell array of values corresponding to each variable

% load file into cell array
fid = fopen(path,'r');
text = textscan(fid,'%s');
text = text{1};             

data = cell(size(var));     % instantiate data

% extract data
for i = 1:length(var)
    if i == length(var) % CASE OF 'P1' = 'P' in acqu
        varName = strcat('##$',var{i}(1),'=');
        % Get data from text
        for j = 1:length(text)
            if strcmp(varName,text{j})
                 data(i) = text(j+3);
            end
        end
    else  % other variables
        varName = strcat('##$',var{i},'=');
        % Get data from text
        for j = 1:length(text)
            if strcmp(varName,text{j})
                 data(i) = text(j+1);
            end
        end
    end
end

fclose(fid);

function D = getFolders(location,revertFlag)
% GETFOLDERS returns the list of non-system folders at specified location.
% If revertFlag is active (-1), then getFodlers returns the list of files
% at specified location.

% Input: 
%   location    (optional) path/name of folder
%   reverFlag   (optional) numeric value. See function overview above
% Output:
%   D           cell array of strings -- folder/file names
% -------
% Example: list = getFolders(pwd);  % folders in current directory
%          list = getFolders(1);    % files in current directory

% get folder contents
if exist('location','var') && isa(location,'char')
    currentFolder = pwd;    % save current location
    cd(location);           % go to custom location
    D = struct(dir);        % get folder contents
    cd(currentFolder);      % return to original location
else
    D = struct(dir);        % get contents of current directory
end

% check for revert flag
if ~exist('revertFlag','var')       % no second argument
    if exist('location','var') && isa(location,'numeric')
        revertFlag = location;      % first argument is revertFlag
    else
        revertFlag = 0;
    end
end
    
% Get directory names
i = 1;          
e = length(D);
while i <= e
    if revertFlag == 0
        % remove system folders and files from list
        if (strcmp(D(i).name(1),'.') || ~D(i).isdir)
            D(i) = [];
            e = e-1;
        else
            i = i+1;
        end
    else
        % remove all folders
        if D(i).isdir
            D(i) = [];
            e = e-1;
        else
            i = i+1;
        end
    end
end

temp = [];
for i = 1:length(D)
    temp = [temp;{D(i).name}];
end
D = temp;

function save(data,name,saveTo,flag)
%% SAVE writes data to a CSV file, at the specified path
% Input:
%   data        data to be written (cell format)
%   name        name of CSV file
%   saveTo      path where CSV file is saved
%   flag        determines how workspace variables are structured

home = pwd;                     % save caller's path
cd(saveTo);                     % go to save location
name = strcat(name,'.csv');     % add extension

% write first row
fid = fopen(name, 'w');         % open file and override any existing data
for row = 1:size(data,1)
    for col = 1:size(data,2)
        fprintf(fid,'%s,',data{row,col});     % print line
    end
    fprintf(fid,'\n');      % go to next line
end
fclose(fid);    % close file

cd(home);                       % go to caller's path

% SAVE data to workspace
if ~flag     % save to struct
    for i = 1:size(data,2)
        s.(genvarname(data{1,i})) = toNum(data(2:end,i));
    end
    assignin('base',name(1:end-4),s);
else        % save directly to workspace
    for i = 1:size(data,2)
        assignin('base',data{1,i},toNum(data(2:end,i)));
    end
end

function [nVector] = toNum(vector)
% Converts the cell vector to a numerical vector. If it can't it'll return
% the original vector
nVector = [];
try
    % covert to numerical vector
    for i = 1:length(vector)
        nVector = [nVector;str2num(vector{i})];
    end
catch
    % revert to original
    nVector = vector;
end
    
        