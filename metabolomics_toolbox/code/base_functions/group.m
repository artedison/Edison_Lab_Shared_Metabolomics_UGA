% GROUP reads a specified Microsoft Excel file containing randomized sample 
% list, and imports the list into a single column variable Y. Groups are 
% separated by numerical identifiers; A = 0, B = 1, ... , Z = 25
%
% Dawit Woldegiorgis 10/7/2014
%--------------------------------------------------------------------------

function [Y] = group(filename);
% Input:
%   filename = name of excel document
% Output:
%   Y = single column double with numbers denoting groups


[~,txt] = xlsread(filename);    % load excel as text matrix
for i = 1:size(txt,1)           % find randomized sample list
    for j = 1:size(txt,2)
        if strcmp(txt{i,j},'Randomized Sample List')
            group_cell = txt(i+1:end,j);    % extract column of groups
        end
    end
end

for i = 1:size(group_cell,1)        % convert to number representation
    group = group_cell{i,1}(1);     % extract letter of group
    for j = 0:25
        if group == char(65+j)      % compare to alphabet iteratively
            Y(i,1) = j;             % set value of group
        end
    end
end