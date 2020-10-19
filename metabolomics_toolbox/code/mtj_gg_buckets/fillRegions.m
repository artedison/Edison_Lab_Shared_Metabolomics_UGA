function [regionsCells,concatenatedInds,featureNumbers] = fillRegions(ROIs)

%% fillRegions

% Author: MTJ
% Version: 0.1
% Date: 2020
%
% Description:
%
%   Takes n Regions of Interest (ROIs) as ppm indices for their
%   boundaries, and returns the full index list between the boundaries in 
%   an n x 1 cell array. Default is 2 x n for ROIs. Also provides other
%   useful outputs. 
%
% Inputs:
%
%     ROIs                  2 x n 
%
% Output:
%       
%     regionsCells          n x 1 cell array containing expanded indices
%                           for each region/bucket
%     concatenatedInds      expanded indices, unlisted into a vector.
%                           Useful for linear indexing/vectorizing. 
%     featureNumbers        feature indices corresponding to
%                           concatenatedInds
%
% Usage: 
%         
%         regionsCells = fillRegions([1243;5236,124:756]);
%         [regionsCells,concatenatedInds,featureNumbers] = fillRegions(ROIs)
%                 
% MTJ 2020


    % Determine orientation of matrix 
        binDim = find((size(ROIs))==2,1);  % find the dimension that = 2 (if 2 x 2, first dim is preferred)
        
    % Convert each ROI to a cell
        cells = num2cell(ROIs,binDim);     % make cells according to proper dim
    
    % Apply fillRegion to each cell    
        regionsCells = cellfun(@(x) x(1):x(2),cells,'UniformOutput',0);

    % Generate vect
        concatenatedInds = [regionsCells{:}];
    
    % Generate feature numbers to correspond to concatenatedInds
    
        tmp = regionsCells;
        for i = 1:length(regionsCells)
            tmp{i}(1:end) = i;
        end
        
        featureNumbers = [tmp{:}];
        
end
