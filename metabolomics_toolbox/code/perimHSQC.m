function out=perimHSQC(file,showPlot)
if ~(exist('showPlot','var')==1)
    showPlot=0;
end
assert(strcmp(file(end-2:end),'.ft'),'File must be a ft file')

spectra=pipe2matlab(file);
ppm1=spectra.ppm1;
ppm2=spectra.ppm2;
spectra.real=abs(spectra.real);
%spectra.real(spectra.real<1e4)=0;
if showPlot,figure,imagesc(spectra.ppm1,spectra.ppm2,spectra.real>5e5),colormap('blue'),title('Orig'),end
img=mat2gray(spectra.real);
img=img>=.1;
img=bwmorph(img,'clean');

maxSize=max(size(img));
ppm1 = linspace(min(ppm1),max(ppm1),maxSize);
ppm2 = linspace(min(ppm2),max(ppm2),maxSize);
img=imresize(img,[maxSize,maxSize]);

img=imclose(img,strel('square',2));
perim=bwperim(img);

out.perim=flipud(fliplr(double(perim)));
out.ppmH=ppm1;
out.ppmX=ppm2;

if showPlot,figure,imagesc(out.ppmH,out.ppmX,out.perim),colormap('blue'),title('perim'),end
end