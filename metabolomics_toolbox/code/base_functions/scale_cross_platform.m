function scaleX=scale_cross_platform(X_NMR,data)

% Author: Rahil Taujale
% Version: 0.1
% Date: 01/30/2017

% Description:
%       Finds the max peak in 2 spectra and scales the peaks in 2nd to the 
%       intensity of the 1st one.
%
% Input: 
%       X_NMR   : stack of 1D spectra
%       data    : 2nd dataset to be scaled to X_NMR
%
% Output: 
%       Data matrix with scaled intensities for the 2nd data.
%
% Example run: data_scaled=scale_cross_platform(X,data);

    mXNMR=max(max(X_NMR));
    mXdata=max(max(data));
    factor=mXdata/mXNMR;
    scaleX=(data/factor);
end
