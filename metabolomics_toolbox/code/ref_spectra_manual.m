function spectra = ref_spectra_manual(spectra, ppmrange, ppm, thresh)
% ************
% Reference spectra by clicking the reference peak in each spectrum
% Required arguments
%     spectra - Structure array containing fields: real, ppm
% Optional arguments:
%     ppmrange - [-0.2, 0.2] (default)
%     thresh - 2 (default)
%     ppm - 0 (default)
% *******************
% Author: Greg Stupp May 2014
% modified by Chaevien: added ppm in case you don't want to reference to
% zero
% MJ edit 2JAN2018: Happy New Year! Print progress to Command Window (Line
%                   29)


if ~exist('ppmrange','var')
    ppmrange = [-0.2, 0.2];
end

if ~exist('thresh','var')
    thresh = 2;
end

if ~exist('ppm','var')
    ppm = 0;
end

for i=1:length(spectra);
    fprintf(['\n\t\tReferencing Spectrum ',num2str(i),'/',num2str(length(spectra))])
    threshI = thresh*abs(median(spectra(i).real));
    peaks=peakpick(spectra(i).real,30,threshI);
    plot(spectra(i).ppm,spectra(i).real,'k'), hold,
    plot(spectra(i).ppm(peaks),spectra(i).real(peaks),'ro')
    set(gca,'xlim',ppmrange)
    set(gca,'xdir','rev')
    [x,y]=ginput(1);
    close
    [~,idx]=min(abs(x-spectra(i).ppm(peaks)));
    ref_ppm(i)=spectra(i).ppm(peaks(idx));
    corrected_ref(i)=ppm-ref_ppm(i);
    spectra(i).ppm=spectra(i).ppm+corrected_ref(i);
end
