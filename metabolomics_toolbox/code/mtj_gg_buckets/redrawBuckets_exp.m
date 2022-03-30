function [updatedBuckets,selBounds] = redrawBuckets_exp(currentBuckets)

%
% MTJ 2020

        
    % Select one or more buckets by drawing box

        [selectedBucketInds,selBounds]  = selectBuckets(currentBuckets); 
            
    if size(currentBuckets,1) > 1

            % See if the new bounds fall in old buckets

                lmb = selBounds(1) > currentBuckets(:,1) & selBounds(1) < currentBuckets(:,2);
                umb = selBounds(2) > currentBuckets(:,1) & selBounds(2) < currentBuckets(:,2);

                % Going to have issues if expanding into two buckets
                % If so, replace the closer bound with the new bound
                    currentBuckets(end+1,:) = [currentBuckets(lmb,1),selBounds(1)];
                    currentBuckets(end+1,:) = [selBounds(2),currentBuckets(umb,2)];

            % Replace selected buckets with one using selection bounds

                currentBuckets(end+1,:) = selBounds;

            % Clean up the old buckets

                currentBuckets(selectedBucketInds,:) = [];       
                updatedBuckets = currentBuckets;
    else 
        
        updatedBuckets = selBounds;
        
    end
end