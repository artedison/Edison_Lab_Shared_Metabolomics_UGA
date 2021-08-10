function [newmat,newtimes,totaltime,resolution,p] = interp2dmatrix(matrix,ppm,timepoints,resampFact,interpMethod)
%% interp2dmatrix

% Author: Michael T. Judge
% Version: 0.1
% Date: MAY2021

% Description:
%       Resample values based on interpolation (for unevenly spaced
%       timepoints in CIVM datasets, or HPLC fraction data). 
%
% Input: 
%       matrix: spectroscopic data, where spectra = rows, ppms = columns
%       ppm: ppm vector corresponding to columns of matrix
%       timepoints: time vector corresponding to rows of matrix
%
% Output: 
%       newmat:         resampled (smoothed) matrix
%       newtimes:       resampling timepoints
%       resolution:     
%
% Log:
%
% Example run: 

%%

if ~exist('interpMethod','var')
    interpMethod = 'linear';
end
if ~exist('resamplingFactor','var')
    resampFact = 1;
end

p = reportParams('exclude',{'matrix','ppm','timepoints'});

    % Interp strategy
      
        % Linearize ppm (dim1) and time (dim2) coordinates

        linppms = meshgrid(ppm,1:size(matrix,1));
        lintimes = meshgrid(timepoints,1:length(ppm))';
       
        % Define timepoints to resample from
            % Linearly Interpolate timepoints based on the minimum
            % resolution (smallest spacing between two timepoints)
            
            resolution = min(diff(timepoints))/resampFact;
            numpoints = ceil(max(timepoints)/resolution);
            
                alltimes = meshgrid(linspace(min(timepoints),...
                                            max(timepoints),...
                                            numpoints),...
                                    1:length(ppm))';
                                
                allppms = meshgrid(ppm,1:numpoints);

        newmat = interp2(linppms,lintimes,matrix,allppms,alltimes,interpMethod);
        newtimes = alltimes(:,1)';
        totaltime = alltimes(end)-alltimes(1);
        
            % figure,imshow(scale(matrix,'pareto'))
            % figure,imshow(scale(newmat,'pareto'))
% V = single(X(200:300,1:25));
% Display the image region.
% 
% imagesc(V);
% axis off
% title('Original Image')
    % Report the fit of the interpolation, as well as the %
    % interpolated for returned data. 
    
    % Optimal resampling scheme is likely of the same (minimal) resolution
    % as the original data. 
    
    % Report the hotspots
    
end