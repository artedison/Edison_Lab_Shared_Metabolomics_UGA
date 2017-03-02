function [XR,ppmR]=remove_ends(X,ppm,shift1,shift2)

% [XR,ppmR]=remove_ends(X,ppm,shift1,shift2)
% 
% Removes region between -Inf and ppm value of "shift1" ,and 
% region between ppm value of "shift2" and Inf.
%
% Arguments:
% X                    Data matrix of spectra
% ppm                  Chemical shift vector corresponding to X
% shift1               upfield ppm of region to be removed
% shift2               downfield ppm of region to be removed
%
% Return Values:
% XR                   Spectral matrix, with region removed
% ppmR                 Chemical shift vector with region removed
%
% Bing Wang

if shift1>shift2 %flip them
    d=shift2;
    shift2=shift1;
    shift1=d;
end

disp(['Removing region from ',num2str(ppm(1)),' ppm to ', num2str(shift1),' ppm']);
disp(['Removing region from ',num2str(shift2),' ppm to ', num2str(ppm(end)),' ppm']);
[~,k1]=min(abs(ppm-shift1));
XR1=X(:,k1:end);
ppmR=ppm(k1:end);
[~,k2]=min(abs(ppm-shift2));
XR=X(:,1:k2);
ppmR=ppmR(1:k2);
