function [newRidges,adjustedRidges] = ridgeCorrection(matrix,currentppm,newRidges,windowWidth,viewWidth,figs)

if ~strcmp(figs,'plotFigs')
     set(0,'DefaultFigureVisible','off');
end

    for ridgenumber = 1:length(newRidges)
        % Get the information about the ridge
%                 matrix = sampleData(samples(sample)).Xcollapsed_1h1d;
%                 currentppm = sampleData(samples(sample)).ppmR_1h1d;
                smoothedWindow = newRidges(ridgenumber).inputs.outputs_ridgeTracing_PeakPick1D.window;
                windowPPMs = newRidges(ridgenumber).inputs.outputs_ridgeTracing_PeakPick1D.windowPPMs;

                ppms = newRidges(ridgenumber).ppms;
                ridgeInds_row = newRidges(ridgenumber).RowInds;
                ridgeIntensities = newRidges(ridgenumber).intensities;   
                ridgeTimes = newRidges(ridgenumber).times;
    % 
                fullTimes = newRidges(ridgenumber).inputs.outputs_ridgeTracing_PeakPick1D.times;

            % Deal with little side chains
                    [newRidge_Intensities,newRidge_times,removedPointInds] = ridge_tracing_pruneRidge(ridgeIntensities,ridgeTimes);
                % Apply the same filter to the position data
                    ridgeInds_row(removedPointInds) = [];
                    ppms(removedPointInds) = [];

            % Handle discontinuous ridges using linear interpolation on smoothed data
                % Simplest method: linear interpolation of peak position
                    [newRidge_rowInds,newRidge_times,newRidge_colInds,newRidge_PPMs,newRidge_intensities,newRidge_gapInds] = ridge_tracing_fillGaps(matrix,currentppm,smoothedWindow,windowPPMs,ppms,ridgeInds_row,newRidge_times,fullTimes,newRidge_Intensities,'plotFigs'); % 'noFigs'   

            % Adjust Ridges to fit original (non-smoothed) matrix
                % Draw the window around the preliminary ridge
%                     windowWidth = 5;
%                     viewWidth = 0.1; % in ppms around ridge extrema
                    [~,~,windowInds,ppmtester] = getAdjacentTrajectories_simpler(matrix,currentppm,newRidge_rowInds,newRidge_colInds,windowWidth,viewWidth);                        

                % Look in the window and get the maximum of each row:
                    [newRidgeLinearInds,newRidge_RowInds,newRidge_ColInds] = adjustRidge(matrix,currentppm,newRidge_rowInds,newRidge_colInds,windowInds,ppmtester,viewWidth);
                    
            %        Fit the ridges iteratively (option)
            %        for i = 1:4             
            %         %% Adjust Ridges to fit original (non-smoothed) matrix
            %             % Draw the window around the preliminary ridge
            %                 windowWidth = 1;
            %                 viewWidth = 0.1; % in ppms around ridge extrema
            %                 set(0,'DefaultFigureVisible','off')
            %                 [ridgeInds_linear,contours,windowInds,ppmtester] = getAdjacentTrajectories_simpler(matrix,currentppm,ridgeInds_row,ridgeInds_col,windowWidth,viewWidth);                        
            %                 
            %             % Look in the window and get the maximum of each row:
            %                 [newRidgeLinearInds,newRidge_RowInds,newRidge_ColInds] = adjustRidge(matrix,currentppm,ridgeInds_row,ridgeInds_col,windowInds,ppmtester,viewWidth);
            %                 
            %        end     
            %        set(0,'DefaultFigureVisible','on')
                    
            % Sort the fields:   
                [adjustedRidges(ridgenumber).RowInds,sortOrder] = sort(newRidge_RowInds);  

                newRidgeLinearInds = newRidgeLinearInds(sortOrder);
                newRidge_ColInds = newRidge_ColInds(sortOrder);       
                newRidge_intensities = zeros(1,length(fullTimes));  

            % Add heads and tails to time series (use last value to impute
            % time). These MUST be sorted.
                intensities = matrix(newRidgeLinearInds);
                newRidge_intensities(1:newRidge_RowInds(1)) = intensities(1); % fill to the beginning with the first value
                newRidge_intensities(newRidge_RowInds(1):newRidge_RowInds(end)) = intensities;    % don't forget to concatenate the middle values, which are the actual measured values
                newRidge_intensities(newRidge_RowInds(end):length(fullTimes)) = intensities(end); % fill to the end with the last value
            % Assign to struct:   
                %[adjustedRidges(ridgenumber).RowInds,sortOrder] = sort(newRidge_RowInds); 
                adjustedRidges(ridgenumber).RidgeName = '';
                adjustedRidges(ridgenumber).Annotation = '';
                adjustedRidges(ridgenumber).windowWidth = windowWidth;
                adjustedRidges(ridgenumber).viewWidth = viewWidth;
                adjustedRidges(ridgenumber).LinearInds = newRidgeLinearInds;
                adjustedRidges(ridgenumber).ColumnInds = newRidge_ColInds;
                adjustedRidges(ridgenumber).RowInds = newRidge_RowInds;   
                adjustedRidges(ridgenumber).ppms = currentppm(newRidge_ColInds); 
                adjustedRidges(ridgenumber).ExtendedIntensities = newRidge_intensities;  
                adjustedRidges(ridgenumber).Intensities = intensities;  
                adjustedRidges(ridgenumber).Times = fullTimes(newRidge_RowInds)';  
                adjustedRidges(ridgenumber).FilledGapInds = newRidge_gapInds; % these correspond to the row number (not the point number), so no sorting necessary
                adjustedRidges(ridgenumber).EndPointsAdded = setdiff(1:length(fullTimes),newRidge_RowInds); % these should be correct
    end     

    set(0,'DefaultFigureVisible','on');
    
end