function [newRidge_linearInds,newRidge_RowInds,newRidge_ColInds] = adjustRidge(matrix,currentppm,ridgeInds_row,ridgeInds_col,windowInds,ppmtester,viewWidth)

% Author: Michael T. Judge
% Version: 0.1
% Date: 04/12/2018

% Description:
%       Given a window around a peak across spectra, this function maps the
%       current ridge obtained for smoothed data to unsmoothed data, or
%       some other spectral matrix. This is accomplished by identifying the
%       maximum within the window at each row. 
%       Accessory function
%       to ridgeCorrection.m, applied after
%       getAdjacentTrajectories_simpler.m. 
%
% Input: 
%       matrix: spectroscopic data, where spectra = rows, ppms = columns
%       currentppm: ppm vector for matrix
%       ridgeInds_row: row of matrix for each ridge point
%       ridgeInds_col: currentppm/matrix indices providing peak position in 
%           each row
%       halfWindowWidth: number of data points to include to left and right
%           of the ridge indicating peak position at different rows
%       viewWidth: spectral width considered for plotting (in ppms)
%
% Output: 
%       linInds: linear indices defining the peak position across spectra,
%           such that matrix(linInds) returns the ridge trajectory
%       contours: the result of matrix(linInds), giving matrix intensity at
%               each point in the window
%       windowInds: linear indices defining the window across spectra,
%           such that matrix(windowInds) returns the entire window
%       ppmtester: useful for extracting the ppm positions from windowInds,
%           accomplished by ppmtester(windowInds) 
%
% Log:
%
% Example run: 
%       See ridgeCorrection.m for context:
%       [linInds,contours,windowInds,ppmtester] = getAdjacentTrajectories_simpler(matrix,currentppm,ridgeInds_row,ridgeInds_col,halfWindowWidth,viewWidth);
%       [newRidge_linearInds,newRidge_RowInds,newRidge_ColInds] = adjustRidge(matrix,currentppm,ridgeInds_row,ridgeInds_col,windowInds,ppmtester,viewWidth);


    %% Autoset the ROI:
        ROI = [min(currentppm(ridgeInds_col))-viewWidth,max(currentppm(ridgeInds_col)+viewWidth)];
        ROIinds = matchPPMs(ROI(1),currentppm):matchPPMs(ROI(2),currentppm);

    %% Make the correction:        
        [~,newRidgeInts] = max(matrix(windowInds),[],2);
        newRidge_linearInds = windowInds(sub2ind(size(windowInds),1:size(windowInds,1),newRidgeInts'));

    %% Give the rows and cols in terms of matrix
        [newRidge_RowInds,newRidge_ColInds] = ind2sub(size(matrix),newRidge_linearInds);
        
    %% Plot the result:
        % Check to see that these bounds are right
            figure, hold on
                surf(currentppm(ROIinds),1:size(matrix,1),matrix(:,ROIinds),'FaceColor','Interp');%,hold on
                shading interp
                xlabel('trajectory')
                zlabel('Scaled Intensity')
                ylabel('Time (h)')
                title('Corrected Ridge on Unsmoothed Data')
                set(gca,'xdir','reverse')
                % Plot the ridge in blue:
                  %scatter3(ppmtester(newRidge_linearInds), ridgeInds_row, matrix(newRidge_linearInds))
                  plot3(ppmtester(newRidge_linearInds), ridgeInds_row, matrix(newRidge_linearInds),'b','linewidth',2)
                  plot3(ppmtester(windowInds(:,1)),ridgeInds_row,matrix(windowInds(:,1)),'r','linewidth',1)
                  plot3(ppmtester(windowInds(:,end)),ridgeInds_row,matrix(windowInds(:,end)),'r','linewidth',1)                    
            hold off                 
% 
end                