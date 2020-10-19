function [cells] = castToCells(arr)
%% addReasonableLegend

% Author: MTJ
% Version: 0.1
% Date: 2020
%
% Description
%
%   Takes an array of numericals and casts to cells (element-wise).
%
% Inputs:
%       
%       'arr'   double or logical array to be converted
%
% Output:
%
%       cells   cell array with contents of arr
%       
%
% Usage: 
%         
%         a = [10,25,50,100,250,500;10,25,50,100,250,500;10,25,50,100,250,500];
%         cells = castToCells(a);
%        
%%

    % Ensure that the pull list is a cell array
    if iscell(arr)
        cells = arr;
    else
        if or(isnumeric(arr),islogical(arr))
            cells = num2cell(arr);
            else 
                if or(ischar(arr),isstring(arr))
                    cells = {arr};
                end
        end
    end
end