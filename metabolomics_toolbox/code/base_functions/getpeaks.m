function [data] = getpeaks(output,varargin)
% Select picked peaks from .tab files and put them in COLMAR format:
% GETPEAKS extracts the 1H and 13C columns from mutliple files containing
% peak picked list, and 
% Input: file name of .csv
%        comma separated list of files containing peaks (.tab files)
% Output: 'output'.csv - contains list of C and H peaks

narginchk(2,inf);
data = [];
for k = 1:size(varargin,2)
    data = [data;extract(varargin{k})];
end
csvwrite(strcat(output,'.csv'),data);

end

function [val] = extract(file)
% Open .tab file, extract data and return required columns as a double
% array
val = [];

fid = fopen(file);
if (fid < 0)
    error('extract:badFile',strcat('File ''',file,''' does not exist'));
    return;
end

% read lines of file
i = 1;
while (~feof(fid)) 
    line{i,1} = fgetl(fid); 
    i = i+1; 
end

% delete lines 1-19
line = line(19:end);
val = [];
for i=1:size(line,1)
    cell = strsplit(line{i},' ');         
    cell = cell(~cellfun(@isempty,cell)); % remove empty cells
    val = [val;cellfun(@str2num,cell,'UniformOutput',false)];   % convert to number and concatentate with previous
end

val = val(:,6:7);   % grab 6th and 7th columns

end