function [hits,overlaps] = findOverlappingRegions(bounds,testset,currentppm)
%% findOverlappingRegions

    % This function allows you to find the intersections (overlaps) between
    % two lists of intervals (bounds and testset) on a linear axis
    % (currentppm). hits returns the indices of the testset which partially
    % or completely contain a member of bounds. 
    
    % Initialize Params for testing:
        %currentppm = currentppm;
        %bounds = sort(bounds);
        hits = cell(size(bounds,2),1);
        overlaps = cell(size(bounds,2),1);
        for i = 1:size(testset,2)
            testset(:,i) = sort(testset(:,i));
        end
        for i = 1:size(bounds,2)
            bounds(:,i) = sort(bounds(:,i));
        end
        [~,inds] = sort(sum(testset,1));
        testset = testset(:,inds);
        
        for j = 1:size(bounds,2)                                        % go through all the query features
            % Find overlapping intervals
                containsLB = (double(testset(1,:) < bounds(1,j)) + double(testset(2,:) > bounds(1,j)) > 1);
                containsUB = (double(testset(1,:) < bounds(2,j)) + double(testset(2,:) > bounds(2,j)) > 1);
                intersections = double(containsLB)+double(containsUB);
                %plot(intersections)                                    % 2 indicates containment, 1 indicates overlap
                hits{j} = find(intersections > 0);
                overlaps{j} = cell(1,length(hits{j}));
            % Find the region of overlap between two intervals
                for i = 1:length(hits{j})                               % go through the hits and get the overlap for each
                    overlaps{j}{i} = intersect(matchPPMs(testset(1,hits{j}(i)),currentppm):matchPPMs(testset(2,hits{j}(i)),currentppm), matchPPMs(bounds(1,j),currentppm):matchPPMs(bounds(2,j),currentppm));  
                end
        end
end
