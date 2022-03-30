function [updatedBuckets] = combineBuckets(currentBuckets)
%% combineBuckets

% Author: MTJ
% Version: 0.1
% Date: 2020
%
% Description:
%
%   Allows the user to select a group of buckets on the current active
%   figure axis and modifies a provided bucket list to combine them into
%   one set of boundaries. Outer boundaries of the outermost buckets that
%   were touched by rectangle defined by a mouse click+drag are retained. 
%   * Intended for use with the refineBuckets() function in the opt_bucket
%   pipeline. 
%
% Inputs:
%
%       currentBuckets  n x 2 matrix of bucket boundaries corresponding to
%                       patch objects in current open figure
%
% Output:
%
%       updatedBuckets  n x 2 matrix of new buckets with merge result
%                       accounted for
%
% Usage: 
%         
%       From an open figure with buckets:
%       
%       modifiedBuckets = combineBuckets(currentBuckets);
%                 
% MTJ 2020

    % Check to make sure there are even buckets
    
        if isempty(currentBuckets)
            updatedBuckets = currentBuckets;
            return
        end

    % Select one or more buckets by drawing box
        combBuckets  = selectBuckets(currentBuckets);
        
        oldBuckets = currentBuckets(combBuckets,:);

        % Make a new one using the lowest and highest x val from all
        % oldbuckets
       
            currentBuckets(end+1,:) = [min(oldBuckets(:)),max(oldBuckets(:))];
        
        % Clean up the old buckets
        
            currentBuckets(combBuckets,:) = [];       
            updatedBuckets = currentBuckets;
    
end