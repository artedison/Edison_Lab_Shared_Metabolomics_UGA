function [ppmMatches] = matchPPMs(ppmVector)
%% Find closest ppm values to a set
    ppmMatches = nan(size(ppmVector));
    for i = 1:length(ppmMatches)    
        [~,ppmMatches(i)] = min(abs(ppm-ppmVector(i)));
    end
