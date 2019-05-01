function [intensities,ridgeTimes,toDelete] = ridge_tracing_pruneRidge(intensities,ridgeTimes)
    % Replace any duplicate times with the higher value:
        fullTimes = unique(ridgeTimes);
        repTimePts = find(hist(ridgeTimes,fullTimes)>1); % these find the times in fullTimes that 
        if ~isempty(repTimePts)   
            toDelete = [];
            for k = 1:length(repTimePts)
                 currentPts = find(ridgeTimes==fullTimes(repTimePts(k)) );
                 [~,b] = max(intensities( currentPts ));
                 toDelete = [toDelete;currentPts(find([1:length(currentPts)]~=b))];
            end
            intensities(toDelete) = [];      
            ridgeTimes(toDelete) = [];
        else
            toDelete = [];
            %fprintf('\n\n\t\tAll points have unique times. Returning ''[]'' for the deleted indices\n\n')  
        end    
end