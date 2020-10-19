function [ints,shifts,p,ax] = optimize_Peakpick1D(X,ppm,represent,peakthreshes,mode)
%% optimize_Peakpick1D

% Author: Michael T. Judge
% Version: 0.1
% Date: 2020

% Description:
%       Runs peakpick1D function (with results plotted and axes linked) so
%       a range of different thresholds can be tested and compared. 
%
% Input:
%       X: stack of 1D spectra
%       ppm: chemical shift vector
%       represent: representative spectrum for peakpick1D
%       peakthreshes: list of thresholds to explore (0,1)
%       mode: 'simple' or 'complex' according to peakpick1D
%
% Output:
%       A plot for each peakpick result, linked in ax array 'ax'
%           handles returned in 'ax' 
%       ints: peak intensities from peakpick1D
%       shifts: peak chemical shift (location) 
%       p: reportParams() structure for recording params from calling envt
%       ax: handles for generated plots
%
%
% Usage: 
%   
%              optimize_Peakpick1D(matrix,ppm,'max',0.05:0.05:0.5,'Complex');
%   
%
%%
        p = reportParams();
        
        shifts = cell(1,length(peakthreshes));
        ints = shifts;
        
        for i = 1:length(peakthreshes)
             [ints{i},shifts{i}]= Peakpick1D(X,ppm,represent,peakthreshes(i),mode);
                    xlabel('Chemical Shift (ppm)')
                    ylabel('Signal Intensity')
                    title(['thresh = ',num2str(peakthreshes(i))])
                    ax(i) = gca;
        end 
    
        linkaxes(ax);
end