function [ppmMatches] = matchPPMs(selectPPMs,allPPMs)
%% Find closest ppm values to a set
% MJ 6MAR2017 update: handles multidimensional arrays of ppms
% MJ 22MAY2017 update: use numel() instead of length()
% MJ 26OCT2017 update: make sure allPPMs is a 1D matrix (vector). This does
%                       NOT mean that allPPMs is sorted, which would be the
%                       standard.

    allPPMs = allPPMs(1:numel(allPPMs)); % make sure that allPPMs is 1D
    ppmMatches = nan(size(selectPPMs));
    if ~or(isempty(selectPPMs),isempty(allPPMs))
        for i = 1:numel(ppmMatches) % index
                [~,ppmMatches(i)] = min(abs(allPPMs-selectPPMs(i)));
        end
    else
        ppmMatches = [];
    end
