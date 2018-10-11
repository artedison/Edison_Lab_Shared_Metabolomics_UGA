function [newRidge_rowInds,newRidge_times,newRidge_colInds,newRidge_PPMs,newRidge_intensities,gapInds] = ridge_tracing_fillGaps(matrix,currentppm,window,windowPPMs,ppms,ridgeInds_row,ridgeTimes,fullTimes,ridgeIntensities,plotFigs)

    %% Find the gaps (not including the ends):    
        fullTimesClipped = fullTimes(min(ridgeInds_row):max(ridgeInds_row)); % exclude the ends 
        [gapTimes,gapIndsClipped] = setdiff(fullTimesClipped,ridgeTimes);  
        gapInds = matchPPMs(gapTimes,fullTimes)'; % these inds apply to fullTimes
        
    %% Compute the positions of new points  
        newPoints = zeros(length(gapInds),4);
        if ~isempty(gapInds)              
            for k = 1:length(gapInds)
                %% Linear interpolation
                    % Get the two closest points in ridgeTimes:
                        %[~,inds] = sort(   abs(  ridgeTimes-
                        %fullTimes(gapInds(k)) )   ); % In terms of time. Better to use index
                            %inds = inds(1:2);  % does order matter here? No.                                    
                        inds(1) = max(find(ridgeInds_row < gapInds(k))); % get next lower index in ridgeInds
                        inds(2) = min(find(ridgeInds_row > gapInds(k))); % get next higher index in ridgeInds
                        closestPoints = [ppms(inds);ridgeTimes(inds)']';
                    % Calculate the line between them
                        [p,S] = polyfit(closestPoints(:,2),closestPoints(:,1),1);
                        newPoints(k,3) = gapTimes(k);                           % this gets the gap times
                        newPoints(k,4) = gapInds(k);                            % this gets the gap indices on fullTimes                                                             
                        newPoints(k,1) = polyval(p,newPoints(k,3));             % this gets the ppm value of the new point, given the time
                        newPoints(k,2) = matchPPMs(newPoints(k,1),currentppm);  % this calculates the ppm index on currentppm for the new point
            end            
        else            
            gapInds = 0;
            %fprintf('\n\n\t\tThere were no gaps. Returning ''0'' for the added indices\n\n')  
        end    
        %% Add the new points to the ridge
            
            newRidge_PPMs = [ppms,newPoints(:,1)'];
            newRidge_colInds = [matchPPMs(ppms,currentppm),newPoints(:,2)'];
            newRidge_times = [ridgeTimes;newPoints(:,3)]';
            newRidge_rowInds = [ridgeInds_row';newPoints(:,4)]; 
                        
            % This one is for matrix values
%             newRidge_intensities = matrix(sub2ind(size(matrix),newRidge_rowInds,newRidge_colInds));                              
            % This one gives the smoothed values, provided the window and
            % windowppms
                newRidge_intensities = window(      sub2ind(    size(window),   newRidge_rowInds',  matchPPMs(newRidge_PPMs,windowPPMs)  )       );
                    % in the case of no new points (no gaps), this should
                    % be the same as ridgeIntensities:
                        %isequal(newRidge_intensities,ridgeIntensities)
        %% Sort the ridge according to timepoint
            [newRidge_times,sortInds] = sort(newRidge_times);
            newRidge_colInds = newRidge_colInds(sortInds);
            newRidge_PPMs = newRidge_PPMs(sortInds);
            newRidge_rowInds = newRidge_rowInds(sortInds);
            newRidge_intensities = newRidge_intensities(sortInds);
            %*** NEED TO ADJUST GAPINDS
                % gapIndsClipped contains the row indices of the ridges, IF
                % in order. Thus, no sorting is necessary if all other
                % calculations are correct.
                
        %% Plot the new ridge points in a different color
            % Autoset the ROI:
%                 ROI = [min(currentppm(newRidge_colInds))-0.3,max(currentppm(newRidge_colInds))+0.3];
                 ROI = [min(windowPPMs),max(windowPPMs)];
                 ROIinds = matchPPMs(ROI(1),currentppm):matchPPMs(ROI(2),currentppm);
    if strcmp(plotFigs,'plotFigs')
                figure, hold on
%                     surf(currentppm(ROIinds),1:size(matrix,1),matrix(:,ROIinds),'FaceColor','Interp');%,hold on
                    surf(windowPPMs,fullTimes,window,'FaceColor','Interp');%,hold on
                    shading interp
                    xlabel('trajectory')
                    zlabel('Scaled Intensity')
                    ylabel('Time (h)')
                    title('Gap-Filled Ridge on Smoothed Data')
                    set(gca,'xdir','reverse')
                    % Plot the initial ridge points in blue:
                      %scatter3(ppms, ridgeInds_row,ridgeIntensities)
                      %plot3(currentppm(newRidge_colInds), newRidge_rowInds,newRidge_intensities,'r')
                      scatter3(currentppm(newRidge_colInds), fullTimes(newRidge_rowInds),newRidge_intensities,'k')
                      %scatter3(ppms, fullTimes(ridgeInds_row),ridgeIntensities,'k')
                    % Plot the new ridge points in red:
                      %scatter3(newPoints(:,1), newPoints(:,3),newPtInts,'r')
                      % gapIndsClipped contains the row indices of the ridges, IF in order:  
                        scatter3(currentppm(newRidge_colInds(gapIndsClipped)), fullTimes(newRidge_rowInds(gapIndsClipped)),newRidge_intensities(gapIndsClipped),'r')
                      %scatter3(newRidge_PPMs(gapInds), newRidge_rowInds(gapInds), newRidge_intensities(gapInds),'r')
                hold off 
    end   
end