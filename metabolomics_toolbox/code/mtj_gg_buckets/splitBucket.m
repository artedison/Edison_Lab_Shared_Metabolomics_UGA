    function updatedBuckets = splitBucket(currentBuckets)
%% splitBucket

    % Author: MTJ
    % Version: 0.1
    % Tested on Matlab Version R2020a
    % Date: JUL2020
    %
    % Description:
    %
    %       Allows the user to split an existing bucket by clicking at the 
    %       point where the bucket should be provided. Accessory function 
    %       for refineBuckets() . 
    %
    % Input:
    %       
    %       currentBuckets:     n x 2 array of n bucket bounds (ppm values)
    %
    % Output:
    %
    %       updatedBuckets:     updated currentBuckets in which the
    %                           selected bucket has been split into two
    %                           adjacent buckets. 
    %
    % Log:
    %       
    %
    % Example run:
    %       
    %       updatedBuckets = splitBucket(currentBuckets)
    %       
    %
%%    
            
            % Draw the ROI
            
                [bucketInd,breakPoint] = selectBucket(currentBuckets);
                
                    if ~isempty(bucketInd)
                                                
                        % Add the new bucket to the list

                            updatedBuckets = [currentBuckets;[currentBuckets(bucketInd,1),breakPoint]];    % lower segment goes from old lowerbound to breakpoint
                            updatedBuckets = [updatedBuckets;[breakPoint,currentBuckets(bucketInd,2)]];    % upper segment goes from breakpoint to upperbound
                                                                                   
                            updatedBuckets(bucketInd,:) = []; % remove old bucket
                    else
                        
                        updatedBuckets = currentBuckets;
                        
                    end
                    

    end
