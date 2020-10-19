function [buckets] = filterBuckets_Peaks_opt(ppm,buckets,peaks)

%% filterBuckets_Peaks_opt

% Author: MTJ
% Version: 0.1
% Date: 2020
%
% Description:
%
%   Filters 'buckets' result from optimize_optBucket() based on whether or
%   not each bucket has a peak as provided in peaks structure. Buckets with
%   no peaks are excluded. Calls filterBins_peaks(). 
%
% Inputs:
%
%       ppm         ppm vector
%       buckets     buckets structure, e.g. from optimize_optBucket()
%       peaks       peaks structure, e.g. from Peakpick1D()
%
% Output:
%       
%       buckets     updated buckets structure, with added field in 
%                   buckets.optimization.results
%
% Usage: 
%         
%       [buckets] = filterBuckets_Peaks_opt(ppm,buckets,peaks);
%                 
% MTJ 2020   
%%
        p = reportParams('exclude',{'buckets','ppm'});
        
    % 
        optResults = buckets.optimization.results; % make a copy
    
    % Filter out the bins with no peaks
        for i = 1:length(optResults) % for each optimization result
            [optResults(i).binsWithPeaks,...
                optResults(i).matchingPeaks,...
                optResults(i).binFilt,...
                optResults(i).peakFilt] = filterBins_peaks(optResults(i).I_b,peaks.shifts);

            
            optResults(i).matchingPeakInts = peaks.ints(:,optResults(i).peakFilt);
            optResults(i).matchingPeakInds = matchPPMs(optResults(i).matchingPeaks,ppm);
            
            % find the max peak ind, shift, and intensity in each bin
            
                peakints = optResults(i).matchingPeakInts;
                peakinds = optResults(i).matchingPeakInds;
                optResults(i).matchedPeaks_binNumbers = zeros(size(peakinds));
                
                for j = 1:length(optResults(i).binsWithPeaks)
%                     optOB_out.results(i).binsWithPeaks;
%                         bin = optOB_out.results(i).binsWithPeaks(j,:);
%                         selectInds = peakinds>matchPPMs(bin(1),ppm) & peakinds<matchPPMs(bin(2),ppm);
                        selectInds = peakinds>=matchPPMs(optResults(i).binsWithPeaks(j,1),ppm) & peakinds<=matchPPMs(optResults(i).binsWithPeaks(j,2),ppm);

                        [maxInt,maxMatchInd] = max(peakints(selectInds)); % result indexed in peakinds
                       
%                     maxInds = optOB_out.results(i).matchingPeakInds;
                    optResults(i).maxMatchedPeakInds(j) = maxMatchInd;
                    optResults(i).maxMatchedPeakInts(j) = maxInt;
                    optResults(i).matchedPeaks_binNumbers(selectInds) = j;
                end
                    
        end
        
        % Assign the copy of the data (updated)
            buckets.optimization.results = optResults;
            
        % Store params
            buckets.optimization.filterParams = p;
        
end