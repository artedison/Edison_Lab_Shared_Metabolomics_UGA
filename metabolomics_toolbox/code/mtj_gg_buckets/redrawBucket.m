function [updatedBuckets,selBounds] = redrawBucket(currentBuckets)
%%
    % Author: MTJ
    % Version: 0.1
    % Tested on Matlab Version R2020a
    % Date: JUL2020
    %
    % Description:
    %
    %       Replace selected buckets with one new bucket. Behaves as a  
    %       combined delete and add function. Accessory function for 
    %       refineBuckets().
    %
    % Input:
    %
    %       currentBuckets: array of bucket bounds
    %
    % Output:
    %
    %       updatedBuckets: updated array of bucket bounds
    %       selBounds:      bounds selected by the user
    %
    % Log:
    %
    % Example run:
    %
    %       
    %

    % Select one or more buckets by drawing box
        
        [selectedBuckets,selBounds]  = selectBuckets(currentBuckets);
        
        % Replace selected buckets with one using selection bounds
        
            currentBuckets(end+1,:) = [selBounds];
        
        % Clean up the old buckets
        
            currentBuckets(selectedBuckets,:) = [];       
            updatedBuckets = currentBuckets;
        
end