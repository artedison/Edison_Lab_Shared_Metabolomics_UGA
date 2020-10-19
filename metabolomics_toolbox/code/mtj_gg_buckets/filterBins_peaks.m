function [binsWithPeaks,matchingPeaks,binFilt,peakFilt] = filterBins_peaks(bins,peaks)
%% filterBuckets_Peaks_opt

% Author: MTJ and GG
% Version: 0.1
% Date: 2020
%
% Description:
%
%   Filters individual opt_bucket results, or individual rows from 
%   optimize_optBucket() based on whether or not each bucket has a peak as 
%   provided in peaks structure. Buckets with no peaks are excluded. 
%   Typically called/used within filterBuckets_Peaks_opt(). 
%
% Inputs:
%
%       bins        bucket boundaries (n x 2 matrix, ppm values)
%       peaks       peak list (ppm vals)
%
% Output:
%       
%       binsWithPeaks   bins with peaks
%       matchingPeaks   peaks that are within bins
%       binFilt         logical values indicating which bins have peaks
%       peakFilt        logical values indicating which peaks fall within
%                       bins
%
% Usage: 
%         
%      [binsWithPeaks,matchingPeaks,binFilt,peakFilt]  = filterBins_peaks(binBounds,peaks.shifts);
%                 
%   MTJ & GG 2020

    %% Filter for the bins with peaks
    
        peakFilt = any(peaks >= bins(:,1) & peaks <= bins(:,2),1); % faster (if no overlapping bins is no issue)
        binFilt = any(peaks >= bins(:,1) & peaks <= bins(:,2),2); % faster (if no overlapping bins is no issue)
        binsWithPeaks = bins(binFilt,:);
        matchingPeaks = peaks(peakFilt);
        
end