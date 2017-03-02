function spectra=Loadcomp2D(thresh)

% Uses visspec(spectra) to load plots of spectra with linked axes.
% 
% Arguments:
% thresh (optional)

if exist('thresh','var')==0
    thresh=200;
end

loadallft()
assignin('caller', 'spectra', spectra)
visspec(spectra,thresh)   