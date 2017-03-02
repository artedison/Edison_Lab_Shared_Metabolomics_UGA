function Setup2D(spectra,shiftpoints1,shiftpoints2)

% Setup2D(spectra,shiftpoints)
%
% Create 3D matrix of Bruker NMR data from array of spectra.  This function
% assembles data matrix X from array spectra.  
%
% Arguments:
% spectra              Array of spectra produced by Stack2D.m
% shiftpoints1         Specify F2 chemical shift range eg; [-1 10] for 1H
% shiftpoints2         Specify F1 chemical shift range eg; [-5 170] for 13C
%
% Return Values:
% X                    Spectral matrix
% XTitles              Bruker format spectrum titles
% XNoise               Local Noise matrix
% ppm1                 Chemical shift vector along F2
% ppm2                 Chemical shift vector along F1

% Steven Robinette

%% Assemble data matrix

disp('Creating spectral matrix and chemical shift vectors')
for ind=1:size(spectra,2)
    last1(ind)=[spectra(ind).ppm1(1)];
    last2(ind)=spectra(ind).ppm2(1);
    first1(ind)=[spectra(ind).ppm1(length(spectra(ind).ppm1))];
    first2(ind)=spectra(ind).ppm2(length(spectra(ind).ppm2));
    spectres(ind,:)=size(spectra(ind).real);
end

if exist('shiftpoints1')==0
    shiftpoints1=[max(first1), min(last1)];
    shiftpoints2=[max(first2), min(last2)];
end

if length(spectra)==1;
    resolution=spectres;
else
resolution=max(spectres);
end
ppm2=shiftpoints2(2):(shiftpoints2(1)-shiftpoints2(2))/(resolution(1)-1):shiftpoints2(1);
ppm1=shiftpoints1(2):(shiftpoints1(1)-shiftpoints1(2))/(resolution(2)-1):shiftpoints1(1);
[XI,YI]=meshgrid(1:resolution(2),1:resolution(1));


for ind1=1:size(spectra,2)
    clear h k1 k2 k3 k4
    [h,k1]=min(abs(spectra(ind1).ppm1-shiftpoints1(2)));
    [h,k2]=min(abs(spectra(ind1).ppm1-shiftpoints1(1)));
    [h,k3]=min(abs(spectra(ind1).ppm2-shiftpoints2(2)));
    [h,k4]=min(abs(spectra(ind1).ppm2-shiftpoints2(1)));
    [X2,Y2]=meshgrid(1:(resolution(2)-1)/(k2-k1):resolution(2),1:(resolution(1)-1)/(k4-k3):resolution(1));
    
    X(:,:,ind1)=interp2(X2,Y2,spectra(ind1).real(k3:k4,k1:k2),XI,YI,'spline');
    XNoise(:,:,ind1)=Koradi(X(:,:,ind1));
    XTitles{ind1}=spectra(ind1).Title;
    assignin('caller','X',X);
    assignin('caller','XNoise',XNoise);
    assignin('caller','XTitles',XTitles);
    assignin('caller','ppm1',ppm1);
    assignin('caller','ppm2',ppm2);
end


