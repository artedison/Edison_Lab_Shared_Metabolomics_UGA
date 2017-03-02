function c=randColor(option)
if ~exist('option','var')==1
    option = 1;
end
if option==1
    colors='krgbcmy';
    c=colors(randi([1, length(colors)]));
else
    colors = jet(32);
    c=colors(randi([1, length(colors)]),:);
end