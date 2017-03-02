function colors = YtoRGB(Y)
Y=remGaps(Y);

rgbspec = [1 0 0;0 0 1;0 1 0;0 1 1;1 0 1;1 1 0;0 0 0];
cspec = 'rbgcmyk';
colors=[];
for i=1:length(Y)
    colors = [colors;rgbspec(Y(i),:)];
end