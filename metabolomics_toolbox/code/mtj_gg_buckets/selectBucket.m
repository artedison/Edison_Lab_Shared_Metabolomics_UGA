function [bucketInd,breakPoint] = selectBucket(currentBuckets)
%% selectBucket

    % Author: MTJ
    % Version: 0.1
    % Tested on Matlab Version R2020a
    % Date: JUL2020
    %
    % Description:
    %
    %       Allows the user to click a region of interest (e.g. bucket) 
    %       to select it. Bucket should be drawn as a patch object on an 
    %       active figure. Previous accessory function for refineBuckets() . 
    %
    % Input:
    %       
    %       currentBuckets:     n x 2 array of n bucket bounds (ppm values)
    %
    % Output:
    %
    %       bucketInd:          index of the bucket selected in
    %                           currentBuckets
    %       breakPoint:         the ppm value of the point that was clicked
    %
    % Log:
    %       
    %
    % Example run:
    %       
    %       [bucketInd,breakPoint] = selectBucket(currentBuckets)
    %       
    %
%%
%             title('Click a bucket to select it') % provide instructions

    % Check to make sure there are even buckets
    
        if isempty(currentBuckets) % skip and return empties
            breakPoint = [];
            bucketInd = [];
            return
        end
        [breakPoint,~] = ginput(1);

        % Find the FIRST bucket in the list that matches the selection

            bucketInd = find(currentBuckets(:,1) < breakPoint & currentBuckets(:,2) > breakPoint,...
                1); % Find the FIRST bucket
            
end