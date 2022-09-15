function [X,ppm,XTitles]=Setup1D(spectra,shiftpoints,resolution)
% Author: Steven Robinette
% Version: 0.1
% Date: 30NOV2018
%
% Description:
%       Preprocess Bruker NMR data from array spectra. This function
%       assembles data matrix X from array spectra. Takes
%       the spectra structure (e.g. obtained from loadallft, Load1D, or
%       loadStudyFTdata) and produces a spectral matrix (X), a ppm
%       vector (ppm), and a cell array of strings for titles (Xtitle). 
%
% Input:
%       spectra:	Array of spectra produced by Load1D.m
%       shiftpoints:	(Optional)- ppm domain range ie [10.5 -0.5]
%                   % We rarely use this option
%       resolution:	(Optional)- Data resolution, ie 65536 for 64k data
%                   % We rarely use this option
%
%
% Output:
%       X:  Spectral matrix
%       ppm:    Chemical shift vector
%       Titles: Titles of NMR spectra
%
% Log:    
%       Ver 0.1: RT, MTJ, YW annotated and checked
%
% Example run:
%
%   [X,ppm,XTitles]=Setup1D(spectra,shiftpoints,resolution)
%
    % Get the first and last ppm values for each spectrum, as well as the
    % resolution.
    if nargin==1
        % For each spectrum
        for ind=1:size(spectra,2)           
            last(ind)=[spectra(ind).ppm(1)]; % get the max ppm value
            first(ind)=[spectra(ind).ppm(length(spectra(ind).ppm))]; % get min ppm value
            spectres(ind)=length(spectra(ind).real);    % get the number of data points
        end
        shiftpoints=[max(first), min(last)]; % get the intersection of ppm values across all spectra
        resolution=max(spectres);   % get the highest number of datapoints across all spectra
    end
    
    % Calculate a common ppm vector based on those parameters
    disp('Constructing spectral dataset')
    X=zeros(size(spectra,2),resolution);
    ppm=shiftpoints(2):(shiftpoints(1)-shiftpoints(2))/(resolution-1):shiftpoints(1);
    XTitles=cell(size(spectra,2),1);
    
    % Interpolate the data for each spectrum based on this new ppm vector:
    for k=1:size(spectra,2)
        clear h k1 k2 
        [h,k1]=min(abs(spectra(k).ppm-shiftpoints(2)));
        [h,k2]=min(abs(spectra(k).ppm-shiftpoints(1)));
        
        X(k,:)=interp1(spectra(k).ppm(k1:k2),spectra(k).real(k1:k2),ppm,'spline','extrap');
        XTitles(k)=cellstr(spectra(k).Title);
    end
    % Reverse the ppm order
    ppm=ppm(length(ppm):-1:1);
    X=X(:,length(ppm):-1:1);

end