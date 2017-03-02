function [X,ppm,XTitles]=Setup1D(spectra,shiftpoints,resolution)

% [X,ppm,XTitles]=Setup1D(spectra,shiftpoints,resolution)
%
% Preprocess Bruker NMR data from array spectra.  This function
% assembles data matrix X from array spectra.
%
% Arguments:
% spectra              Array of spectra produced by Load1D.m
% shiftpoints          (Optional)- ppm domain range ie [10.5 -0.5]
% resolution           Data resolution, ie 65536 for 64k data
%
% Return Values:
% X                    Spectral matrix, aligned and normalized if specified
% ppm                  Chemical shift vector
% Titles               Titles of NMR spectra

% Steven Robinette

if nargin==1
    
    for ind=1:size(spectra,2)
        last(ind)=[spectra(ind).ppm(1)];
        first(ind)=[spectra(ind).ppm(length(spectra(ind).ppm))];
        spectres(ind)=length(spectra(ind).real);
    end
    
    shiftpoints=[max(first), min(last)];
    resolution=max(spectres);
end

disp('Constructing spectral dataset')
X=zeros(size(spectra,2),resolution);
ppm=shiftpoints(2):(shiftpoints(1)-shiftpoints(2))/(resolution-1):shiftpoints(1);
XTitles=cell(size(spectra,2),1);
for k=1:size(spectra,2)
    clear h k1 k2
    [h,k1]=min(abs(spectra(k).ppm-shiftpoints(2)));
    [h,k2]=min(abs(spectra(k).ppm-shiftpoints(1)));
    
    X(k,:)=interp1(spectra(k).ppm(k1:k2),spectra(k).real(k1:k2),ppm,'spline','extrap');
    XTitles(k)=cellstr(spectra(k).Title);
end
ppm=ppm(length(ppm):-1:1);
X=X(:,length(ppm):-1:1);



end