%This is another change

ax(1)=subplot(2,1,1);
subplot(2,1,1); plot(ppmMixture,Xmixture(6,:))
 set(gca,'XDir','rev')
xlim([-0.5,200])

for i=1:20
    ax(2)=subplot(2,1,2);
    subplot(2,1,2); plot(ppm,X(i,:))
    title(XTitles(i))
    set(gca,'XDir','rev')
xlim([-0.5,200])

     linkaxes(ax,'x');
   pause
end

%%
Xmedian=median(X,1);
    peaks=peakpick(Xmedian,30,0.03*max(Xmedian)); %0.001 is the threshold
plot(ppm,Xmedian,'k')
    oplot(ppm(peaks),Xmedian(peaks),'ro')
    
    for i=1:length(spectra)
low=find(XALlow(i,:)<2e6);
XALlow(i,low)=0.01;
end
figure, plot(ppm,XALlow)
set(gca,'XDir','reverse')
title('star alignment')

XALSN==0;
XALSN(ans)=0.001;