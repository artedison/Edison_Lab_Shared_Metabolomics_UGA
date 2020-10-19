function [X, axis] = select_region(XOriginal, AxisOriginal,shift1,shift2)
%% select_region

    % Author: MTJ/RMB?
    % Version: 0.1
    % Tested on Matlab Version ?
    % Date: ?
    %
    % Description:
    %       
    %       Takes a pair of ppm values and returns the data from XOriginal
    %       in that region. 
    %
    % Input:
    %
    %       XOriginal: spectra you're pulling data from
    %       AxisOriginal: ppm vector going with XOriginal
    %       shift1: ppm shift boundary 1
    %       shift2: ppm shift boundary 2
    %
    % Output:
    %       
    %       X:      selected region of X
    %       axis:   ppm vector for selected region
    %
    % Log:
    %
    % Example run:
    %
    %       [X, axis] = select_region(XOriginal,AxisOriginal,shift1,shift2);
    %
    %

%% 
%
%   
% 
%   Hope this helps,
%   MJ 6MAR2017

if shift1>shift2 %flip them
    d=shift2;
    shift2=shift1;
    shift1=d;
end

if isinf(shift1)
    shift1=AxisOriginal(1);
end
if isinf(shift2)
    shift2=AxisOriginal(end);
end

disp(['Selecting region from ',num2str(shift1),' to ', num2str(shift2),' from the original Axis']);
[~,k1]=min(abs(AxisOriginal-shift1));
[~,k2]=min(abs(AxisOriginal-shift2));

axis = AxisOriginal(1,k1:k2);

X = XOriginal(:,k1:k2);




