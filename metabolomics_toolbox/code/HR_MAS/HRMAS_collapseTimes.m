function [matrixCollapsed,timesCollapsed,totalTime,resolution] = HRMAS_collapseTimes(matrix,times,binsize)
%% HRMAS_collapseTimes

% Author: Michael T. Judge
% Version: 0.1
% Date: 04/12/2018

% Description:
%       Takes time-series spectra (i.e. HR-MAS) and sums every [binsize] 
%       spectra for improved signal-to noise. 
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
%       [linInds,contours,windowInds,ppmtester] = getAdjacentTrajectories_simpler(matrix,currentppm,ridgeInds_row,ridgeInds_col,halfWindowWidth,viewWidth)

% We can't actually give the desired resolution because the data are 
% discrete. We can only report the resulting resolution.

            matrixCollapsed = zeros(length(binsize:binsize:size(matrix,1)),size(matrix,2));
            timesCollapsed = zeros(size(matrixCollapsed,1),1);
            list = binsize:binsize:size(matrix,1);
            for i = 1:length(list)
                matrixCollapsed(i,:) = sum(matrix(list(i)-(binsize-1):list(i),:),1);
                timesCollapsed(i) = times(list(i));
            end   
    totalTime = num2str(round(times(end)-times(1),1));
    if size(timesCollapsed,1)>1
        resolution = num2str(round((timesCollapsed(2)-timesCollapsed(1))*60,1));
    else
        resolution = 'Resolution is not defined for a single timepoint';
        msgbox('Warning: resolution is not defined for a single timepoint');
    end
end