function Loadcomp1D()

% Uses visspec(spectra) to load plots of spectra with linked axes.

loadallft()
assignin('caller', 'spectra', spectra)
overlayNMR(spectra)


    