function [bucketInds,xbounds]  = selectBuckets(currentBuckets)
%% selectBuckets

    % Author: MTJ
    % Version: 0.1
    % Tested on Matlab Version R2020a
    % Date: JUL2020
    %
    % Description:
    %
    %       Allows the user to select one or more buckets by drawing a 
    %       rectangle to touch or contain them. Accessory function for 
    %       refineBuckets() . Replaces selectBucket() . 
    %
    % Input:
    %       
    %       currentBuckets:     n x 2 array of n bucket bounds (ppm values)
    %
    % Output:
    %
    %       bucketInds:         indices of the bucket selected, in
    %                           currentBuckets
    %       xbounds:            ppm bounds of the region selected
    %
    % Log:
    %       
    %
    % Example run:
    %       
    %       [bucketInds,xbounds] = selectBuckets(currentBuckets)
    %       
    %
%%    
    % Draw the ROI
            
                [xbounds,~] = drawROI();
                
                xbounds = sort(xbounds);
                
            % Find buckets which overlap with those bounds
                    % Find all instances where a currentBucket bound falls
                    % within xbounds
                %bucketInds = find(any(currentBuckets < xbounds(2) & currentBuckets > xbounds(1), 2) );  
                   % If a currentBucket bound falls within xbounds, then an
                   % xbound must fall within a currentBucket. 
                   
                   % If both xbounds are within a bin, then  
                   
                   % At least one xbound is both
                        % greater than the lower bound AND
                        % less than the upper bound
                        
                    bucketInds = find(...
                        any(xbounds < currentBuckets,2) & any(xbounds > currentBuckets,2)...
                        |...
                        any(currentBuckets < xbounds(2) & currentBuckets > xbounds(1), 2)...
                        );
end
