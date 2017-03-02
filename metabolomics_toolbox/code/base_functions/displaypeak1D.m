function displaypeak1D(X,ppm,shift,Y)

% Author: Edison Lab
% Ver 0.2

%displaypeak1D(X,ppm,shift,Y)
%
%Input: X: stack of 1D spectra
%       ppm: chemical shift vector
%       shift: one or more chemical shifts, takes scalar or vector 
%       Y: Optional- if response vector Y is known, use it to color the scores plot
%
% displaypeak plots a set of chemical shift values as circles on a set of
% 1D spectra.  Use this script to assign correlated peaks.

% Edited by : Rahil Taujale
% Date      : 03/02/2017
% Ver       : 0.2
%           Changed colormap to generate it using distinguishable_colors
%           function.
if exist('Y')==0
    Y=zeros(size(X,1),1);
    Ycolor=ones(size(X,1),1);
% else
%     % convert Y values to colors 
%     %         mean-center      range of Y values              why? scale
% %Ycolor=ceil((  [(Y-mean(Y)) / (2.01*max(  abs(Y-mean(Y))  ))]   +.5)*100);
% Ycolor = ceil( [ Y ./ max(Y)-min(Y) + 0.01 ] .* 95);
end
% cmap=jet(100);
cmap=distinguishable_colors(100);

figure, hold on;
for k=1:size(X,1)
%     plot(ppm,X(k,:),'Color',cmap(Ycolor(k),:))
    plot(ppm,X(k,:),'Color',cmap(Y(k),:))
end
set(gca,'XDir','reverse')
for k=1:length(shift)
[a,b(k)]=min(abs(ppm-shift(k)));
end
scatter(shift,max(X(:,b)));