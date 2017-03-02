function h=compare1D(X1,X2,ppm,Y)

%
%Input: X1: the first stack of 1D spectra
%       X2: the second stack of 1D spectra
%       ppm: chemical shift vector
%       Y: Optional- if response vector Y is known, use it to color the scores plot
%
%show1D plots two sets of 1D spectra to compare the effects of baseline correction or alignment.

if exist('Y')==0
    Y=zeros(size(X1,1),1);
    Ycolor=ones(size(X1,1),1);
else
Ycolor=ceil(([(Y-mean(Y))/(2.01*max(abs(Y-mean(Y))))]+.5)*100);
end
cmap=lines(100);

h=figure;
h1=subplot(2,1,1);
hold on;
for k=1:size(X1,1)
    p(k)=plot(ppm,X1(k,:),'Color',cmap(Ycolor(k),:));
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
h2=subplot(2,1,2);
hold on;
for k=1:size(X2,1)
    plot(ppm,X2(k,:),'Color',cmap(Ycolor(k),:));
    
end
set(gca,'XDir','reverse')
linkaxes([h1,h2],'xy')
