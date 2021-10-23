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

% Old method (summing)
%   [matrixCollapsed,timesCollapsed,totalTime,resolution] = HRMAS_collapseTimes(matrix,times,binsize)

% Simple
    
    if ~interleaved
        matrixCollapsed = movmean(matrix,windowsize,1,'Endpoints','shrink');
        timesCollapsed = timepoints;
    else
        
        
    end
    
% Interleaved

    if interleaved
        % Use local interpolation to fill in the gaps for timepoints
            % How many points to interp for?
                                        
            si = varargin{3};
            [~,b] = ismember('noesypr1d',{si.sample(3).expType.type});
                [~,sortInd] = sort(cellfun(@str2num,{si.sample(3).dataFiles.name}));
                si.sample(3).expTypeIndex,
                
                % Check to see if exptKey matches the timepoints
                    theseTps = si.sample(3).exptKey(sortInd) == b;
                    ttps = double(theseTps);
%                     sum(theseTps);  
                    ttps(theseTps) = timepoints;
                    
                % *We know where the gaps are
                % Diff the timepoints
                    dtp = diff(timepoints);
                    
                % Fill each gap with round((elapsedTime) / avgdiff)
                % timepoints
                    % Elapsed time during each gap
                        
                        gaps = find(~theseTps);
                            lgap = gaps-1;
                            rgap = gaps+1;
                                hasboth = lgap>0 & rgap<length(theseTps);
                                rgap = rgap(hasboth);
                                lgap = lgap(hasboth);
                        gapDiffs = ttps(rgap) - ttps(lgap);
%                             figure,plot(gapDiffs)
                            
                 % Avg diff for each normal tp
                    % Remove the diffs at the gaps
                        gaps = gaps - (1:length(gaps));
                        dtp(gaps(gaps<length(dtp))) = [];
                            
                % Generate new timepoints new data
                    meanRes = mean(dtp);
                    newTimes = timepoints(1):meanRes:timepoints(end);
                    
                    
%                     figure,plot(ttps)
%                     figure,plot(dtp)
%                     figure,hold on,plot(timepoints), plot(newTimes)
%                     
                % Interpolate the surface
                    x = 1:size(matrix,2);
                    y = timepoints;
                    v = matrix;
                    
                    [xq,yq] = meshgrid(x, newTimes);
                    newmat = griddata(x,y,v,xq,yq);
                    Vq = interp2(X,Y,V,Xq,Yq);
    end
            
            
        
    totalTime = num2str(round(timepoints(end)-timepoints(1),1));
    if numel(timesCollapsed)>1
        resolution = num2str(round((timesCollapsed(2)-timesCollapsed(1))*60,3));
    else
        resolution = 'Resolution is not defined for a single timepoint';
        msgbox('Warning: resolution is not defined for a single timepoint');
    end
    
    
end