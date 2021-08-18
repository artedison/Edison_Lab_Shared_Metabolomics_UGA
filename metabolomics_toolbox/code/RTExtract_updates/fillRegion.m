function [inds] = fillRegion(bounds,ppm)
% e.g. ROI (ppm) -> vector of indices
% MTJ 2020

    if bounds(2)<bounds(1)
        bounds = flip(bounds);
    end

    inds = matchPPMs(bounds(1),ppm):matchPPMs(bounds(2),ppm);
    
end