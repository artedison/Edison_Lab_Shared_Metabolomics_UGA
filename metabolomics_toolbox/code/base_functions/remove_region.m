function XR=remove_region(X,ppm,shift1,shift2)

    % Author: Edison Lab
    % Ver 0.1
    % Tested on Matlab Version R2017b
    % Date: 25FEB2019
    %
    %
    % Description:
    %   Removes region between ppm value of "shift1" and ppm value of "shift2"
    %   such as water region. You can put Inf or -Inf to specify "to the end".
    %
    % Input:
    %   X: Data matrix of spectra
    %   ppm: Chemical shift vector corresponding to X
    %   shift1: upfield ppm of region to be removed
    %   shift2: downfield ppm of region to be removed
    %
    % Output:
    %   XR: Spectral matrix, with region removed (ppm vector not modified)
    %
    % Log:
    %   Edited by : MTJ,LM,YW,SZ
    %   Date      : 25FEB2019
    %   Ver       : 0.1
    %   YW documented 10/10/2018
    % Example run:
    %


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
