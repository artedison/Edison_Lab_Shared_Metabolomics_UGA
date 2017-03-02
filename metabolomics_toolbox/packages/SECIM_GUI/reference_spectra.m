function reference_spectra(sample_number, axes1, spectra, ppmrange, thresh)
% ************
% Reference spectra by clicking the TSP peak in each spectrum
% Required arguments
%     spectra - Structure array containing fields: real, ppm
% Optional arguments:
%     ppmrange - [-0.2, 0.2] (default)
% *******************
if ~exist('ppmrange','var')
    ppmrange = [-0.2, 0.2];
end

if ~exist('thresh','var')
    thresh = 2;
end

%for i=1:length(spectra);
%     threshI = thresh*abs(median(spectra(sample_number).real));
%     peaks=peakpick(spectra(sample_number).real,30,threshI);
     plot(axes1,spectra(sample_number).ppm,spectra(sample_number).real,'k')%, hold,
%     %h = get(gcf,'CurrentAxes')
%     plot(axes1, spectra(sample_number).ppm(peaks),spectra(sample_number).real(peaks),'ro')
     set(axes1,'xlim', ppmrange)
     set(axes1,'xdir','rev')
%     %[x,y]=ginput(1);
%     hold off;
%     %[~,idx]=min(abs(x-spectra(i).ppm(peaks)));
%     %TSP_ppm(i)=spectra(i).ppm(peaks(idx));
%     %spectra(i).ppm=spectra(i).ppm-TSP_ppm(i);
end
%sample_number = sample_number + 1;
%disp(sample_number);