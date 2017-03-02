function readacqu(parent,filename)
%READACQU reads acqu data and creates for required parameters 
%  The function parses through the parent folder to identify the location
%  of the acqu file. Once identified, the funtion parses through the file
%  looking for parameters and their values. Data is collected from all acuq
%  files found in subfolders of parent folder and the result is written to
%  an excel table.
%
%   Input: parent = main folder name
%          filename(optional) = name of excel file to be created.
%                               default filename is 'Untitled.xlsx'
%
%   IMPORTANT: Function is applicable only to parent folder that is 3 levels
%   above the data. For example, consider the following structure:
%   DATA [DATA1, DATA2, DATA3]
%      where: DATA is parent folder and DATA(1-3) are its contents
%   Expanding the contents let DATA and all its subcontents be noted by
%   DATA[DATA1[1[acqu],2[acqu],3[acqu]], DATA2[1[acqu], 2[acqu]], DATA3[1[acqu]]]
%       thus one acqu file can be mapped as
%       DATA -> DATA1 -> 1 -> acqu, where 1 is 1st level, DATA1 is at 2nd
%       level and DATA is at 3rd level in the folder hierarchy
%
%   For the above example, calling readacqu('DATA') will create an excel 
%   spreadsheet contining 3 sheets - namely DATA1, DATA2, and DATA3 and 
%   the associated data
%       Note: The first sheet in excel will be an empty sheet, which must
%       be deleted manually. 


if nargin == 1
    filename = 'Untitled.xlsx';
end

S = subfolders(parent);
E = [];
last = size(S,2);

% Collect data from all subfolders
for i = 1:last
    newParent = strcat(parent,'/',S{i});
    E{i} = sheet(newParent);
    if (size(E{i},1) == 1)
        E(i) = [];
    end
end

warning ('off','all');  % silent mode
% Write to excel file
last = size(E,2);
for i = 1:last
    xlswrite(filename,E{i},S{i});
end
warning ('on','all'); % normal mode
end

%--------------------------------------------------------------------------

%=====================%
%   HELPER FUNCTIONS  %
%=====================%

%--------------------------------------------------------------------------

function [D] = sheet(parent)
% SHEET creates an excel sheet comprising the combined data extracted from
% each subfolder of parent. Returns a strucutre.

S = subfolders(parent);

% list of variables required
var = {'SW','RG','TD','PULPROG','NS','DS','PE','EXP','O1','SFO1','SOLVENT','DW','D1','TD0','PL','NUC1','TE'};
% name of data file
file = 'acqus';

% base folder
if isempty(S) || isequal(S,{'pdata'}) || isequal(S,{'pdata2'}) || isequal(S,{'pdata3'})
    % check if file exist
    files = cellstr(ls(parent));    % list of files in parent
    all = size(files,1);            % all files
    for i = 1:all
        if strcmp(files(i),file)    % search for file from list of files
            loc = strcat(parent,'/',file);
            D = extract(loc,var);  % return data extracted
            row = relative(parent);
            D = [row,D];
            return
        end
    end
end

% Parse through subfolders
numRow = size(S,2);
numCol = size(var,2)+1;
D = cell(numRow,numCol);  % combined data
for i = 1:numRow
    newParent = strcat(parent,'/',S{i});
    D(i,:) = sheet(newParent);     % add each data in separate rows
end
Col = [{''},var];   % column data
D = [Col;D];
end

%--------------------------------------------------------------------------

function [NAME] = relative(parent)
% RELATIVE returns the folder name, provided a path to the folder
last = size(parent,2);
pos = 1;
for i = 1:last
    index = last-i+1;
    if isstrprop(parent(index),'punct')
        pos = index+1;
        break;
    end
end
NAME = parent(pos:end);
end

%--------------------------------------------------------------------------

function [S] = subfolders(parent)
% SUBFOLDERS returns a cell containing all immediate subfolder names.
    D = dir(parent);
    sub = [D.isdir];
    S = {D(sub).name};
    i = 1;
    col = size(S,2);
    while i <= col
        if (strcmp(S{i},'.') || strcmp(S{i},'..'))
            S(i) = [];          % delete cell
            col = size(S,2);   % update size
        else
            i = i+1;
        end; 
    end
end

%--------------------------------------------------------------------------

function [S] = extract(file, var)
% EXTRACT takes the variables listed in var and searches for them in the
% file. When located, the data of the variables are stored in a cell array
% and returned.

% open file into a cell array
A = textscan(fopen(file,'r'),'%s');
A = A{1};

% create cell
S = cell(size(var));

% fill data
for i = 1:size(var,2)
    find = false;
    for j = 1:size(A,1)
        if strcmp(strcat('##$',var{i},'='),A{j})
            S(i) = A(j+1);
            find = true;
        end
    end
end

fclose('all');

end

%--------------------------------------------------------------------------
