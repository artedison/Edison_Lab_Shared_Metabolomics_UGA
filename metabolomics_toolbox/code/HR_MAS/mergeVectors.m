function [unifiedVect, inds_unified_inEach,indsCells,repeatedVals] = mergeVectors(vectList)

% Merge similar, overlapping vectors (e.g. ppm vectors, or time vectors)
% into a mean vector in the overlapping region. Inds of each list are
% provided in a matrix, in which the rows correspond to the cells of
% vectList. Still not sure if this is the best way to do this, yet. 
% 
% MTJ 2021
% a = 1:11;
% b = 2:13;
% c = 0:10;
% vectList = {a,b,c};

% Borrowed from Setup1D (interpolate data):

%     for ind=1:size(length(vectList)   
%         last(ind)=[spectra(ind).ppm(1)]; % get the max ppm value
%         first(ind)=[spectra(ind).ppm(length(spectra(ind).ppm))]; % get min ppm value
%         spectres(ind)=length(spectra(ind).real);    % get the number of data points
%     end
% 
%     shiftpoints=[max(first), min(last)]; % get the intersection of ppm values across all spectra
%     resolution=max(spectres);   % get the highest number of datapoints across all spectra
% 
%
% Calculate a common ppm vector based on those parameters
% 
%     disp('Constructing spectral dataset')
%     X=zeros(size(spectra,2),resolution);
%     ppm=shiftpoints(2):(shiftpoints(1)-shiftpoints(2))/(resolution-1):shiftpoints(1);
%     
%     
%     
%     % Interpolate the data for each spectrum based on this new ppm vector:
%     
%     for k=1:size(spectra,2)
%         clear h k1 k2 
%         [h,k1]= matchPPMs(shiftpoints(2),spectra(k).ppm);
%         [h,k2]= matchPPMs(shiftpoints(1),spectra(k).ppm);
%         
%         X(k,:)=interp1(spectra(k).ppm(k1:k2),spectra(k).real(k1:k2),ppm,'spline','extrap');
%     end


% Find overlapping range across all vects

    bounds = [cellfun(@min,vectList);cellfun(@max,vectList)];
    
% Using average resolution, fill out those times - this is the unified vector

    meanRes = mean(cellfun(@range, vectList)./cellfun(@length,vectList));
    unifiedVect = linspace(max(bounds(1,:)),min(bounds(2,:)),  round(  ( min(bounds(2,:))-max(bounds(1,:)) ) / meanRes)  );
    
% get the inds in timevects for unifiedTimes, store the inds

    inds_unified_inEach = zeros(length(vectList),length(unifiedVect));
    
    for i = 1:length(vectList)
        inds_unified_inEach(i,:) = matchPPMs(unifiedVect,vectList{i});
    end
    
    % Check for any repeated values (happens sometimes)
    
        [x,y] = find(~diff(inds_unified_inEach,1,2));
        repeatedVals.vect = x;
        repeatedVals.ind = y;
        
%     for i = 1:length(vectList)
%         vectList{i}(inds_unified_inEach(i,:))
%     end

    indsCells = num2cell(inds_unified_inEach,2);
        
end
