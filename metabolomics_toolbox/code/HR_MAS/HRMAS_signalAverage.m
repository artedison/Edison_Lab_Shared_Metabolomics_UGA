function [matrixCollapsed,timesCollapsed,totalTime,resolution] = HRMAS_signalAverage(matrix,timepoints,windowsize,varargin)
%% HRMAS_signalAverage

% Author: Michael T. Judge
% Version: 0.1
% Date: 13FEB2019

% Description:
%       Takes time-series spectra (i.e. HR-MAS) and does a moving average 
%       along the time dimension for improved signal-to noise. Higher
%       resolution than simply summing consecutive spectra as done in 
%       HRMAS_collapseTimes. Do not use this to average a matrix that spans
%       multiple treatments (each must be done separately).
%
% Input: 
%       matrix: spectroscopic data, where spectra = rows, ppms = columns
%       times: timepoint vector (should be complete and linear)
%       binsize: number of spectra to consider in the averaging
%
%       varargin: 
%
%            'interleaved', spacing (integer)
%
% Output: 
%       matrixCollapsed:    time-averaged matrix
%       timesCollapsed:     time-averaged time vect
%       totalTime:          experiment length
%       resolution:         collapsed resolution
%
% Log:
%
% Example run: 

%%

interleaved = false;

    if ~isempty(varargin)
        ind = strcmp('interleaved',varargin);
        if ~isempty(ind)
            interleaved = true;
            spacing = varargin{ind+1};
        end
    end

%%

% Simple
    
    if ~interleaved
        matrixCollapsed = movmean(matrix,windowsize,1,'Endpoints','shrink');
        timesCollapsed = timepoints;
    else
        
        
    end
        
    totalTime = num2str(round(timepoints(end)-timepoints(1),1));
    if numel(timesCollapsed)>1
        resolution = num2str(round((timesCollapsed(2)-timesCollapsed(1))*60,3));
    else
        resolution = 'Resolution is not defined for a single timepoint';
        msgbox('Warning: resolution is not defined for a single timepoint');
    end
    
% Interleaved

    
    
end