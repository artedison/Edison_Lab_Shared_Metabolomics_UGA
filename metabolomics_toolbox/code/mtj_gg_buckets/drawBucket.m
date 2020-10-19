    function updatedBuckets = drawBucket(currentBuckets)
            
            % Draw the ROI
                [bounds,ys] = drawROI();
                
            % Add the new bucket to the list
            
                updatedBuckets = [currentBuckets;bounds];
        
    end
