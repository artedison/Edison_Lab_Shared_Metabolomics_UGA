function [ridges] = ridgeTracing_clusterPeaks_interactive_2(input,numberOfRidges,timeWeight,ppmWeight,intensityWeight)
%% Filtering by PeakPick1D (simple)
% Usage:
%
%     ppm = sampleData(1).ppmR_1h1d;
%     matrix = sampleData(1).Xcollapsed_1h1d;
%     timepoints = sampleData(1).timesCollapsed_1h1d;
%     ROI = [2.5,2.8];
%     peakPickThreshold = 0.9;
%     numberOfRidges = 8; %this is ~number of peaks picked*1.25
%     timeWeight = 4;
%     ppmWeight = 1;
%     intensityWeight = 1;
%
%
%     [ridges] = ridgeTracing_PeakPick1D(matrix,ppm,timepoints,ROI,peakPickThreshold,numberOfRidges,timeWeight,ppmWeight,intensityWeight);

% MJ and YW 20APR2018

%% Trace the ridges

        peakInds = input.peakInds;
        window = input.window;
        windowPPMs = input.windowPPMs;
        timepoints = input.times;
        ROI = input.ROI;
        ROIinds = input.ROIinds;
        peakPickThreshold = input.peakPickThreshold;
        ppm = input.ppm;
        matrix = input.matrix;
        inputLittleVars = input;
        inputLittleVars.ppm = [];
        inputLittleVars.matrix = [];
    % Agglomerative clustering of peaks based on position, where time dimension
    % is weighted higher and also considering peak height

        heights = window(peakInds);
        rows = input.rows; % in terms of the window
        cols = input.cols; % in terms of the window
            %figure,scatter3(cols,rows,heights)
        ridgePoints = [(rows/max(rows)/timeWeight)',(cols/max(cols)/ppmWeight)',((heights/max(heights))/intensityWeight)'];
        % Count Ridges
            %numberOfRidges = ceil(median(cellfun(@length,col))*1.25);

        clusters = clusterdata(ridgePoints,numberOfRidges);%'Distance','euclidean','Linkage','single');
        %figure,scatter3(cols,rows,heights,100,clusters,'filled')
% Transform clusters to ridges

%figure,plot3(cols(inds),rows(inds),heights(inds))



        clusterNumbers = unique(clusters);
    ridges = struct();
    for i = 1:length(clusterNumbers)
        % Get the indices in 'window'
            inds = find(clusters==clusterNumbers(i)); % Problem here?: linearIndices should not be in this order
            ridges(i).RowInds = rows(inds);         % in terms of the window
            %ridges(i).ColumnInds = cols(inds);      % in terms of the window
            ridges(i).ColumnInds = cols(inds) + matchPPMs(min(windowPPMs),ppm);
            ridges(i).WindowIndices = sub2ind(size(window),ridges(i).RowInds,cols(inds)); % in terms of the window
        % Get the ridge info (matrix context)
            ridges(i).ppms = windowPPMs(cols(inds));
            ridges(i).times = timepoints(ridges(i).RowInds);
            ridges(i).intensities = window(ridges(i).WindowIndices);
        % Record other relevant info
            ridges(i).inputs.outputs_ridgeTracing_PeakPick1D = inputLittleVars; % (did this in ridgeTracing_PeakPick1D)
            ridges(i).inputs.timeWeight = timeWeight;
            ridges(i).inputs.ppmWeight = ppmWeight;
            ridges(i).inputs.intensityWeight = intensityWeight;

        % Remove spurious ridges (not necessary)
             % Second derivative filter in ppm dimension

             % Filter out ridges for which there are two datapoints at a
                % single timepoint (THIS IS NECESSARY)

             %
    end

% Plot the resulting ridges, and select the desired ridges by clicking.
    % The highlighted ones get selected.
    % Click once to select. Click again to de-select
    %fig = plotRidges(window,ppm,ROIinds,timepoints,ridges,'Interactive Ridge Picking. Press any key to get menu.');
%         f2 = figure; hold on
%             surf(ppm(ROIinds(1:end)),timepoints,window,'FaceColor','Interp');%,hold on
%             shading interp
%             set(gca,'xdir','reverse')
%             xlabel('ppm')
%             zlabel('Scaled Intensity')
%             ylabel('Time (h)')
%             %title(titleStr)
%             for i = 1:length(ridges)
%                 h(i) = plot3(ridges(i).ppms,ridges(i).times,ridges(i).intensities,'LineWidth',4);
%             end

    %% Do interactive stuff
            answer = 0;
            while ~or(answer==2,answer==4)
                plotRidges(window,ppm,ROIinds,timepoints,ridges,'Interactive Ridge Picking');
                answer = menu('Interactive Ridge Picking Menu','Pick Clusters to Join','Pick Final Clusters','Delete Ridge(s)','Cancel');
                % Globals MUST be cleaned up after each iteration and
                % initialized ONCE in the code.
                    global clickedRidges;
                        clickedRidges = [];
                    global lineNumber;
                        lineNumber = 1;
                switch answer
                    case 1
                                                    
                            ridges = joinRidges(ridges);
                            %%
                                    % Get the ridges to join:
                                    %                             title('Select two ridges to join by clicking on them, then hit Return (follows maxima between the points):')
                                    %                                 clickedRidges = [];
                                    %                                 lineNumber = 1;
                                    %                             selectLine(gcf)
                                    %                             pause()
                                                            % Convert to useful indices
                                                            % variables set automatically?)


                                    % % %                             inds = clickedRidges;
                                    % % %                             inds = inds(inds~=1)-1;  % shift (surface plot = 1)
                                    % % %                             % Get the odd-frequency ridges, only take these
                                    % % %                             clear('clickedRidges','lineNumber')
                                    % % %                                 ridgesToConnect = unique(inds(find(mod(sum(inds==inds'),2))));
                                    % % %                                         newRidge = joinRidges(ridges,ridgesToConnect,input,window);
                                                                            % Concatenate with the ridges
                                                                            % structure:
                                                                                %ridges = [ridges newRidge];
                                                                            % Get rid of the old ridges
                                                                                %ridges(ridgesToConnect) = [];

                    case 2
                        % Pick the final ridges
                        
                            title('Select the final ridges by clicking on them, then hit Return')
                            
                            % Develop ridge clicking method based on
                            % general proximity (not exact clicks)
                            
                                [ridgeInds] = clickRidge(ridges);
                                close(gcf)                                                 
                            
                        % Save them
                            if ~isempty(ridgeInds)
                                ridges = ridges(ridgeInds);
                            end
                            
                            
                            
                    case 3
                        % Delete ridge(s):
                            title('Select the final ridges by clicking on them, then hit Return')
                                clickedRidges=[];
                                lineNumber=1;
                            selectLine(gcf)
                            pause()
                            inds = clickedRidges;
                            inds = inds(inds~=1)-1;  % shift (surface plot = 1)
                            % Count the number of odd ridges
                                oddRidges = inds(find(mod(sum(inds==inds'),2)));
                            ridges(oddRidges) = [];
                            ridgePoints(oddRidges,:) = [];
                            % Redo the clustering
                            clusters = clusterdata(ridgePoints,numberOfRidges);%'Distance','euclidean','Linkage','single');
                            %clusters = clusterdata(ridgePoints,numberOfRidges-length(oddRidges));%'Distance','euclidean','Linkage','single');
                            % Remake ridges structure?
%                                 ridges = struct();
%                                 for i = 1:length(clusterNumbers)
%                                     % Get the indices in 'window'
%                                         inds = find(clusters==clusterNumbers(i)); % Problem here?: linearIndices should not be in this order
%                                         ridges(i).RowInds = rows(inds);         % in terms of the window
%                                         %ridges(i).ColumnInds = cols(inds);      % in terms of the window
%                                         ridges(i).ColumnInds = cols(inds) + matchPPMs(min(windowPPMs),ppm);
%                                         ridges(i).WindowIndices = sub2ind(size(window),ridges(i).RowInds,cols(inds)); % in terms of the window
%                                     % Get the ridge info (matrix context)
%                                         ridges(i).ppms = windowPPMs(cols(inds));
%                                         ridges(i).times = timepoints(ridges(i).RowInds);
%                                         ridges(i).intensities = window(ridges(i).WindowIndices);
%                                     % Record other relevant info
%                                         ridges(i).inputs.outputs_ridgeTracing_PeakPick1D = inputLittleVars; % (did this in ridgeTracing_PeakPick1D)
%                                         ridges(i).inputs.timeWeight = 1;
%                                         ridges(i).inputs.ppmWeight = .05;
%                                         ridges(i).inputs.intensityWeight = 30;
%
%                                     % Remove spurious ridges (not necessary)
%                                          % Second derivative filter in ppm dimension
%
%                                          % Filter out ridges for which there are two datapoints at a
%                                             % single timepoint (THIS IS NECESSARY)
%
%                                          %
%                                 end


                    case 4
                        % Quit without plotting
                            close(gcf)
                        break
                end
                clear('clickedRidges','lineNumber')
                close(gcf)
            end
            if answer == 2
                        % Plot and save the selected ridges:
                            plotRidges(window,ppm,ROIinds,timepoints,ridges,'Final Ridges. Program will exit now. Run again to trace more stuff.');
            end

            end

    function [h] = plotRidges(window,ppm,ROIinds,timepoints,ridges,titleStr)
        h = figure; hold on
            surf(ppm(ROIinds(1:end)),timepoints,window,'FaceColor','Interp');%,hold on
            shading interp
            set(gca,'xdir','reverse')
            xlabel('ppm')
            zlabel('Smoothed Intensity')
            ylabel('Time (h)')
                set(gcf, 'InvertHardCopy', 'off');
                set(gca,'fontsize',20)
                set(gca,'box','off')
                title('Clustered peaks on Gaussian-smoothed spectra')
            title(titleStr)
            for i = 1:length(ridges)
                h(i) = plot3(ridges(i).ppms,ridges(i).times,ridges(i).intensities,'LineWidth',5);
            end
    end

    function [ridges] = joinRidges(ridges)
            
            title('Join Ridges')
            
            m(1) = msgbox([{'Join ridges into one by doing one of the following'};...
               {'1) Press ''b'' for box drawing mode (in development)'};...
               {'2) Press any other letter key to activate click'}],'WindowStyle','modal');
           
            waitforbuttonpress;
                key = get(gcf,'CurrentCharacter');
                if strcmp(key,'b')
                    waitfor(msgbox('Option ''b'' - Box drawing mode - is not allowed at the moment, still under development','WindowStyle','modal'));
                    % Box method
%                         % Get box 
%                             [xbds,ybds] = drawROI(); 
% 
%                         % Collect all ridge points in box
% 
%             %                 [~,lineObjs] = extractLineData();
%                             [~,...
%                                 lineObjData,...
%                                 inds_myLines,...
%                                 ~] = mapLines2Lines(findall(gca,'Type','Line'),{ridges.ppms},...
%                                                                                                {ridges.times});
% 
%                             inXbounds = xbds(1)<lineObjData(:,1) & xbds(2)>lineObjData(:,1);
%                             inYbounds = ybds(1)<lineObjData(:,2) & ybds(2)>lineObjData(:,2);
% 
%                             ridgesToConnect = unique(inds_myLines(inXbounds & inYbounds));
                            
                else % Use click method
                        % Click one or more ridges to join
                        
                            m(3) = msgbox([{'Join Ridges:'};...
                                   {'1) Press any letter key to activate click'};...
                                   {'2) Click a ridge to select/deselect'};...
                                   {'Repeat to select ridges to join into one.'};...
                                   {'Press ''Enter/Return'' to update data and return to menu'}],'WindowStyle','modal');
                            
                            % Develop ridge clicking method based on
                            % general proximity (not exact clicks)

                                [ridgeInds,~,inds_myLines,~,lineObjData] = clickRidge(ridges);
                                close(gcf)                                                 
                            
                        % Update ridges
                            if ~isempty(ridgeInds)
                               matches = inds_myLines == ridgeInds';
                               newRidgePts = any(matches,2);
                               newppms = lineObjData(newRidgePts,1);
                               newtimes = lineObjData(newRidgePts,2);
                                   [t,~,uniquetimeInds] = unique(newtimes); % times need to be unique, not ppms. Unique default is sorted
                                   compmat = zeros(length(uniquetimeInds),length(ridgeInds));
                                   
                                   % Find overlapping times
                                        % Inds to fill out matrix
                                            [~,ridgePtID] = find(matches);
                                            lininds = sub2ind(size(compmat),uniquetimeInds,ridgePtID);
                                            
                                            % Get the intensities
                                                ints = {ridges.intensities};
                                                ints = [ints{:}]';
                                                
                                            compmat(lininds) = ints(lininds);
                                            
                                            
                                        % Pick the highest ridge point at
                                        % each time
                                        
                                            [newInts,maxinds] = max(compmat,[],2);
                                            
                                        % Apply this choice to ppms and
                                            ridges(end+1).times = t; % has to be after [ints{:}] to avoid empty issues
                                            ridges(end).ppms = newppms(maxinds)';
                                            ridges(end).intensities = newInts';
                                            ridges(end).joined = ridges(ridgeInds);

                                            ridges(ridgeInds) = [];
                            else
                                waitfor(msgbox('No ridges were selected','WindowStyle','modal'))
                                
                            end                    
                    
                    
                end
            
        % clean up msgboxes
            
            close(m(ishandle(m)))
            
    end
    
    
    
    
    
    
    
    
    
    
%%    
%     function [ridgeNumber] = pickRidges
%         clickedRidges = [];
%         lineNumber = 1;
%         selectLine(gcf)
%         pause()
%         % Convert to useful indices
%         % variables set automatically?)
%         inds = clickedRidges;
%         inds = inds(inds~=1)-1;  % shift (surface plot = 1)
%     end

% function [newRidge] = joinRidges(ridges,ridgesToConnect,input,window)
% 
%     figure,hold on,surf(input.ppm(input.ROIinds),input.times,window),shading interp
%             scatter3(input.ppm(v1(2,:)),input.times(v1(1,:)),v1(3,:),'r')
%             scatter3(input.ppm(v2(2,:)),input.times(v2(1,:)),v2(3,:),'y')
% 
% %% This is behaving as expected
% %     v2 = v1(:,1:5);
% %     v1 = v1(:,6:15);
% %     figure,hold on,surf(input.ppm(input.ROIinds),input.times,window),shading interp
% %             scatter3(input.ppm(v1(2,:)),input.times(v1(1,:)),v1(3,:),'r')
% %             scatter3(input.ppm(v2(2,:)),input.times(v2(1,:)),v2(3,:),'y')
% %     figure,hold on,surf(input.ppm(input.ROIinds),input.times,window),shading interp
% %             scatter3(input.ppm(v1(2,1)),input.times(v1(1,1)),v1(3,1),'r')
% %             scatter3(input.ppm(v2(2,end)),input.times(v2(1,end)),v2(3,end),'y')
% 
%     % Find points where the ridges are closest. Get the indices of those points.
%             allRidgePointDistances = pdist2(v1',v2');
%             %[allRidgePointDistances,minLinearInd] = pdist2(v1',v2','euclidean','Smallest',1); % pdist2 will return the smallest distance and its index
%         [~,minLinearInd] = min(allRidgePointDistances(1:end));
%         [v1Ind,v2Ind] = ind2sub(size(allRidgePointDistances),minLinearInd); % these index the points on each respective ridge
% %%
%     % Define a window by the closest points of the current ridges
%             ridges(ridgesToConnect(1)).ppms(v1Ind)
%             ridges(ridgesToConnect(2)).ppms(v2Ind)
% 
%         connectWindowPPMBounds = sort([,]);
%             connectWindowPPMBoundInds = matchPPMs(connectWindowPPMBounds,input.ppm);
% 
%             ridges(ridgesToConnect(1)).RowInds(v1Ind)
%             ridges(ridgesToConnect(2)).RowInds(v2Ind)
%         [connectWindowTimeBounds,timeSort] = sort([,]);
% 
%             % *** timeSort is used to ensure that the ridges are ordered correctly when concatenated later on
%     figure,hold on,surf(input.ppm(input.ROIinds),input.times,window),shading interp
%             scatter3(input.ppm(v1(2,:)),input.times(v1(1,:)),v1(3,:),'r')
%             scatter3(input.ppm(v2(2,:)),input.times(v2(1,:)),v2(3,:),'y')
% 
%  %% Positive control for the correct region:
%     a = [29732,29714,29706,29718,29779,29709,29709,29718,29775,29805,29830,29861,29885,29904];
%     b = [2,3,4,5,6,7,8,9,10,11,12,13,14,15];
% %     figure,hold on,surf(input.ppm(a),input.times(b),input.matrix(b,a),'FaceColor','Interp');
% %     %scatter3(input.ppm(windowCols(newPointsWindowInds)),input.times(windowRows),newPointsIntensities,'*r');
% 
%     aa = matchPPMs(input.ppm(a),input.ppm(input.ROIinds));
%     figure,hold on,surf(input.ppm(a),input.times(b),window(b,aa),'FaceColor','Interp');
%     %scatter3(input.ppm(windowCols(newPointsWindowInds)),input.times(windowRows),newPointsIntensities,'*r');
%    %%
%     % If the ridges are adjacent in the ppm or time dimension, simply combine
%     % them. Otherwise, follow the highest ridge between them.
% 
%         if or(abs(connectWindowTimeBounds(1)-connectWindowTimeBounds(2))<=1    ,   abs(connectWindowPPMBoundInds(1)-connectWindowPPMBoundInds(2))<=1)
%             % Concatenate and re-sort the ridges. No Window necessary.
% 
%         else
%             % Chop off the first and last point
%                 windowRows = connectWindowTimeBounds(1):connectWindowTimeBounds(2);
%                 windowRows = windowRows(2:end-1);
%                 % (We can certainly have two points on the same ppm, but not the same
%                 % timepoint.)
%                 windowCols = connectWindowPPMBoundInds(1):connectWindowPPMBoundInds(2);
%              % Make the window
%                 %connectWindow = input.matrix(windowRows,windowCols); %this is not the right matrix
%                 connectWindow = window(windowRows,matchPPMs(input.ppm(windowCols),input.ppm(windowCols)));
%              % Find the maximum of each row in that window
%                 [newPointsIntensities,newPointsWindowInds] = max(connectWindow,[],2);
%     plotRidges(connectWindow,input.ppm,windowCols,input.times(windowRows),ridges,'test');
%     figure,hold on,surf(input.ppm(windowCols),input.times(windowRows),connectWindow,'FaceColor','Interp');
%     scatter3(input.ppm(windowCols(newPointsWindowInds)),input.times(windowRows),newPointsIntensities,'*r');
%                 newPointsPPMs = input.ppm(newPointsWindowInds + matchPPMs(connectWindowPPMBounds(1),input.ppm)); %%%% this could be incorrect offset, need to check
%                                                                                                                  % should be correct if first timepoint was removed
%                 %a = windowCols(newPointsWindowInds);
%                 %b = windowRows;
%             % Concatenate and re-sort ridges
%                 % Big issue: how do you arrange
%                 % the points such that they
%                 % follow a single line?
%                     % We collect all the points
%                     % into a vector (just the
%                     % rows), then sort them and
%                     %
%                         %newRidge.RowInds = unique(ridges(ridgesToConnect(1)).RowInds);
%                 newRidge = struct();
%                 newRidge.RowInds = [ridges(ridgesToConnect(timeSort(1))).RowInds,...
%                                     windowRows,...
%                                     ridges(ridgesToConnect(timeSort(2))).RowInds];
%                 newRidge.ColumnInds = [ridges(ridgesToConnect(timeSort(1))).ColumnInds,...
%                                        matchPPMs(newPointsPPMs,input.ppm),...
%                                        ridges(ridgesToConnect(timeSort(2))).ColumnInds];
%                 newRidge.WindowIndices = ridges(ridgesToConnect(1)).WindowIndices;
%                 newRidge.ppms = [ridges(ridgesToConnect(timeSort(1))).ppms,...
%                                  newPointsPPMs,...
%                                  ridges(ridgesToConnect(timeSort(2))).ppms];
%                 newRidge.times = input.times(newRidge.RowInds);
%                 newRidge.intensities = [ridges(ridgesToConnect(timeSort(1))).intensities,...
%                                         newPointsIntensities',...
%                                         ridges(ridgesToConnect(timeSort(2))).intensities];
% %                 [newRidge.RowInds,sortInds] = sort(newRidge.RowInds);
% %                 newRidge.ColumnInds = newRidge.ColumnInds(sortInds);
% %                 newRidge.ppms = newRidge.ppms(sortInds);
% %                 newRidge.times = newRidge.times(sortInds);
% %                 newRidge.intensities = newRidge.intensities(sortInds);
%         end
% end

% function [newRidge] = expandRidge(ridges,input)
%     answer = 0;
%     while answer~=3
%         % Regenerate the figure with the expanded ridge
%         plotRidges(window,ppm,ROIinds,timepoints,ridges,'Interactive Ridge Picking');
%         answer = menu('Interactive Ridge Picking Menu','Pick Clusters to Join','Pick Final Clusters','Delete Ridge(s)','Cancel');
%         % Draw the ROI
%
%         % Get the indices of that region in terms of whole matrix
%
%         % Calculate the maxima
%
%         % Figure out which ridge is being added to
%
%         % Add the maxima to the ridge
%
%         % Sort by row index
%
%         close(gcf)
%     end


%end
