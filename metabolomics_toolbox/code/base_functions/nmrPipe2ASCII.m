function nmrPipe2ASCII()
loadallft
for i=1:length(spectra)
    s=spectra(i);
    x=s.ppm;
    x(:,2)=s.real';
    dlmwrite([s.Title, '.txt'],x,'delimiter','\t','precision','%.6f');
    x=[];
end