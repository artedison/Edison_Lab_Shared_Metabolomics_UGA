function h=display1D(X,ppm,Y)

%display1D(X,ppm,Y)
%
%Input: X: stack of 1D spectra
%       ppm: chemical shift vector
%       Y: Optional- if response vector Y is known, use it to color the scores plot
%
% display 1D NMR spectra

if exist('Y')==0
    Y=zeros(size(X,1),1);
    Ycolor=ones(size(X,1),1);
else
Ycolor=ceil(([(Y-mean(Y))/(2.01*max(abs(Y-mean(Y))))]+.5)*100);
end
cmap=jet(100);

h=figure; hold on;
for k=1:size(X,1)
    p(k)=plot(ppm,X(k,:),'Color',cmap(Ycolor(k),:));
end
set(gca,'XDir','reverse')
hold off;
Yu=unique(Y);
for i=1:length(Yu)
    g(i)=hggroup;
    set(p(Y==Yu(i)),'Parent',g(i));
end
grps=('A':char(64+length(Yu)));
legend(g,grps(:));
