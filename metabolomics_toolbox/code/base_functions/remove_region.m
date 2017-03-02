function XR=remove_region(X,ppm,shift1,shift2)

% XR=remove_region(X,ppm,shift1,shift2)
% 
% Removes region between ppm value of "shift1" and ppm value of "shift2"
% such as water region. You can put Inf or -Inf to specify "to the end".
%
% Arguments:
% X                    Data matrix of spectra
% ppm                  Chemical shift vector corresponding to X
% shift1               upfield ppm of region to be removed
% shift2               downfield ppm of region to be removed
%
% Return Values:
% XR                   Spectral matrix, with region removed

if shift1>shift2 %flip them
    d=shift2;
    shift2=shift1;
    shift1=d;
end

if isinf(shift1)
    shift1=ppm(1);
end
if isinf(shift2)
    shift2=ppm(end);
end

XR=X;
disp(['Removing region from ',num2str(shift1),' ppm to ', num2str(shift2),' ppm']);
[~,k1]=min(abs(ppm-shift1));
[~,k2]=min(abs(ppm-shift2));
XR(:,k1:k2)=zeros(size(X,1),k2-k1+1);
