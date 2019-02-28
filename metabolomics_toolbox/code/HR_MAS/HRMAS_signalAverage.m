function [matrixCollapsed,timesCollapsed,totalTime,resolution] = HRMAS_signalAverage(matrix,times,binsize)
%% HRMAS_signalAverage

% Author: Michael T. Judge
% Version: 0.1
% Date: 13FEB2019

% Description:
%       Takes time-series spectra (i.e. HR-MAS) and does a moving average 
%       along the time dimension for improved signal-to noise. 
%
% Input: 
%       matrix: spectroscopic data, where spectra = rows, ppms = columns
%       times: timepoint vector (should be complete and linear)
%       binsize: number of spectra to consider in the averaging
%       model: type of smoothing to do
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

    matrixCollapsed = movmean(matrix,binsize,1);
    timesCollapsed = times;
    
    totalTime = num2str(round(times(end)-times(1),1));
    if size(timesCollapsed,1)>1
        resolution = num2str(round((timesCollapsed(2)-timesCollapsed(1))*60,3));
    else
        resolution = 'Resolution is not defined for a single timepoint';
        msgbox('Warning: resolution is not defined for a single timepoint');
    end
end