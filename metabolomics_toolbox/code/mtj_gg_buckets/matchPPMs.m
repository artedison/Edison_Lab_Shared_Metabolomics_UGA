function [ppmMatches] = matchPPMs(selectPPMs,allPPMs)
%% 
    % Author: Michael T. Judge
    % Version: 0.4
    % Tested on Matlab Versions 2017b-2020a
    % Date: 24AUG2020
    %
    % Description:
    %       Converts array of ppm values to the indices of their closest 
    %       matches in a complete ppm vector. 
    %       ppm(ppmMatches) converts back to ppms. 
    %       Uses min(abs(difference)).
    %
    % Input:
    %       selectPPMs  array of ppm values to convert
    %       allPPMs     whole ppm vector
    %
    % Output:
    %       ppmMatches  indices of the matches in the ppm vector. 
    %
    % Log:
    %       MJ 6MAR2017 update: handles multidimensional arrays of ppms
    %       MJ 22MAY2017 update: use numel() instead of length()
    %       MJ 26OCT2017 update: make sure allPPMs is a 1D matrix (vector). 
    %       MJ 24AUG2020 update: inline allppms vector unlisting
    %
    % Example run:
    %       
    %       [ppmMatches] = matchPPMs([1.23,7.5,9.4],ppms);
    %       [ppmMatches] = matchPPMs([1.23,7.5,9.4;1.4,7.6,9.55],ppms);
    %
%%


    allPPMs = allPPMs(:); % make sure that allPPMs is 1D
    ppmMatches = nan(size(selectPPMs));
    if ~or(isempty(selectPPMs),isempty(allPPMs))
        for i = 1:numel(ppmMatches) % index
                [~,ppmMatches(i)] = min(abs(allPPMs-selectPPMs(i)));
        end
    else
        ppmMatches = [];
    end
