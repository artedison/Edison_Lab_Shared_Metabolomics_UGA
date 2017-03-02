function overlayNMR1D(spectra)
%accepts a structure array containing real, ppm, Title

colors=('ymcrgbkrbgymk');
numSpectra=length(spectra);
%if length(spectra)>7
    %disp('Only showing first 7')
    %numSpectra=7;
%end
figure('Position',[50 50 1800 900])
subplot(2,numSpectra,1)

for i=1:numSpectra
    h(i)=subplot(2,numSpectra,i);
    plot(spectra(i).ppm, spectra(i).real, colors(i))
    set(gca,'XDir','rev')
end

[n,xout]=hist(spectra(1).real,100);
offset=xout(1)*10;
middle=round(numSpectra+(numSpectra/2));
h(numSpectra+1)=subplot(2,numSpectra,numSpectra+1:middle);
hold on
for i=1:numSpectra
    plot(spectra(i).ppm,(offset*i-offset)+spectra(i).real,colors(i))
    set(gca,'XDir','rev')
end
hold off

h(numSpectra+2)=subplot(2,numSpectra,middle+1:numSpectra*2);
hold on
for i=1:numSpectra
    plot(spectra(i).ppm,spectra(i).real,colors(i))
    set(gca,'XDir','rev')
end
pause(1)
linkaxes([h(:)],'x')
set([h(:)],'ylimmode','manual');







