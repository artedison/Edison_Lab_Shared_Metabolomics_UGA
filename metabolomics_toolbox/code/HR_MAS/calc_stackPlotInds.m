function [plotInds,plotIndsCat] = calc_stackPlotInds(data,numPoints,varargin)
% Calculates the inds needed for plotting one or more datasets together at
% a given resolution (e.g. number of spectra). This should make it easier
% to plot reasonable numbers of spectra using slow functions like
% stackSpectra().

% Takes something like 
%
% >> data
%
% ans = 
%
%   1×5 cell array
% 
%     {1×32768 double}    {1×32768 double}    {120×32768 double}    {974×32768 double}    {916×32768 double}
% 
%   and pulls row indices from each cell (sample group) to enable
%   plotting. These indices apply as follows:
%
%   matrix = [data];
%   matrix( [plotInds{:}] , :)
%
% This could also be worked out for timepoint - based (continuous)
% distribution of inds. 

% MTJ JAN2021


%% Parse optional args

%    optParams = {'maxInd','minInd'};
%    varargin =[ {'minInd'}  {0}  {'maxInd'}];
%    
%     if ~isempty(varargin)
%         isstr = find(cellfun(@ischar,varargin));
%         [a,b] = ismember(varargin(isstr),optParams);
%         for parm = isstr(a)
%             varargin(parm)
%         end
%         ismember(op
%         
%     end

%% Convert data to rows vect

    type = class(data);
    switch type
        case 'double'
            
            % Simply distribute on rows
            
                plotInds = {round(linspace(1,size(data,1),numPoints))};
                plotIndsCat = plotInds;
                
        case 'cell'
            
            % Get the size of each dataset
            
                rows = cellfun(@(x) size(x,1), data);
                cols = cellfun(@(x) size(x,2), data);
                
                
            % Calculate the inds lists
            
                cr = cumsum(rows);    
                
                if ~exist('maxInd','var')
                    maxInd = max(cr); 
                end
                
                indsLists = fillRegions([1,cr(1:end-1) + 1;...
                                         cr                  ]);
                    % Note: this operation is equivalent to:
                    %                   [     1 : cr(1)
                    %                   cr(1)+1 : cr(2)
                    %                   cr(2)+1 : cr(3)
                    %                   ...
                    %                   cr(end-1)+1 : cr(end)]
                                 
                 % Trim this indsLists to maxInd
                    
                    for i = 1:length(indsLists)
                        indsLists{i} = indsLists{i}(indsLists{i} < maxInd);
                    end
                    rows(cellfun(@isempty,indsLists)) = [];
                    indsLists(cellfun(@isempty,indsLists)) = [];
                    
            % Distribute numPoints across each range
                % If there are > 1 timepoint, make sure at least the first
                % and last points are represented
                
                numPoints_row = zeros(size(rows));
                numPoints_row(rows==1) = 1;

%                 numPoints_row = ones(size(rows));
%                 numPoints_row(rows>1) = 2;
%                 
            % How many points are left to distribute?
                remPoints = numPoints - sum(numPoints_row);     % number of points left to distribute over, accounting for single-point runs
                needPoints = rows>1;
                  
            % Divide these up to the groups (always include endpoints)
            
                ratios = rows(needPoints) / sum(rows(needPoints));
                numPoints_row(needPoints) = round(ratios * remPoints);
                numPoints_row(numPoints_row<1) = 1;
                numPoints_row(rows>1 & numPoints_row < 2) = 2;
                
                    % In the future, we could rank the remainders of these and assign remaining
                    % points based on highest remainders
                    
                plotInds = cell(size(rows));
                
            % Get inds for each group
                for i = 1:length(numPoints_row)
                    plotInds{i} = round(linspace(indsLists{i}(1),indsLists{i}(end),numPoints_row(i)));
                end      
                
        case 'struct'
            
        otherwise
            error('Input format is incorrect. ''data'' should be one matrix (double, m x n), or multiple matrices (cell or struct, each element m x n).')
    end
    
    plotIndsCat = [plotInds{:}]';
    
end
