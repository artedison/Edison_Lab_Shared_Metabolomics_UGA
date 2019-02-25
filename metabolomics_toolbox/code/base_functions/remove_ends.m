function [XR,ppmR]=remove_ends(X,ppm,shift1,shift2)

    % Author: Bing Wang
    % Ver 0.1
    % Tested on Matlab Version R2017b
    % Date: 25FEB2019
    %
    %
    % Description:
    %   Removes region between -Inf and ppm value of "shift1",and
    %   region between ppm value of "shift2" and Inf. "shift1" < "shift2".
    %
    % Input:
    %   X: Data matrix of spectra
    %   ppm: Chemical shift vector corresponding to X
    %   shift1: upfield ppm of region to be removed
    %   shift2: downfield ppm of region to be removed
    %
    % Output:
    %   XR: Spectral matrix, with region removed
    %   ppmR:Chemical shift vector with region removed
    %
    % Log:
    %   Edited by : MTJ,LM,YW,SZ
    %   Date      : 25FEB2019
    %   Ver       : 0.1
    % Example run:
    %

    disp(['Removing region from ',num2str(ppm(1)),' ppm to ', num2str(shift1),' ppm']);
    disp(['Removing region from ',num2str(shift2),' ppm to ', num2str(ppm(end)),' ppm']);
    [~,k1]=min(abs(ppm-shift1));
    XR1=X(:,k1:end);
    ppmR=ppm(k1:end);
    [~,k2]=min(abs(ppmR-shift2));
    XR=XR1(:,1:k2);
    ppmR=ppmR(1:k2);
end