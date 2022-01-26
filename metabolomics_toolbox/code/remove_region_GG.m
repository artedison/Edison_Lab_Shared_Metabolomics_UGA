function XR=remove_region_GG(X,ppm,shift1,shift2)

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
% MJ    ** Future: simulate spectral noise along the added line**
% YW documeted 10/10/2018
% modified by GG to create a flat line between the removed regions 


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
for i=1:size(X(:,1))
XR(i,k1:k2)=zeros(size(X(i),1),k2-k1+1) + linspace ( X(i,k1), X(i,k2),(k2-k1+1) );
end
