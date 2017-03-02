function vis2D(spectra,disp,thresh)
% 
% vis2D(spectra,disp,thresh)
% 
% Display plot of spectrum after loading spectra.  
%
% Arguments:
% spectra            Index of spectra to see: ex. spectra(1)
% disp		         'f' for fast, 'c' for contour
% thresh		    Contour/color threshold- try 5 to start

if exist('thresh')==0
    thresh=5;
end

if disp(1)=='f'
cmap =[0         0    0.5625
         0         0    0.6250
         0         0    0.6875
         0         0    0.7500
         0         0    0.8125
         0         0    0.8750
         0         0    0.9375
         0         0    1.0000
         0    0.0625    1.0000
         0    0.1250    1.0000
         0    0.1875    1.0000
         0    0.2500    1.0000
         0    0.3125    1.0000
         0    0.3750    1.0000
         0    0.4375    1.0000
         0    0.5000    1.0000
         0    0.5625    1.0000
         0    0.6250    1.0000
         0    0.6875    1.0000
         0    0.7500    1.0000
         0    0.8125    1.0000
         0    0.8750    1.0000
         0    0.9375    1.0000
         0    1.0000    1.0000
    0.1250    1.0000    1.0000
    0.2500    1.0000    1.0000
    0.3750    1.0000    1.0000
    0.5000    1.0000    1.0000
    0.6250    1.0000    1.0000
    0.7500    1.0000    1.0000
    0.8750    1.0000    1.0000
    1.0000    1.0000    1.0000
    1.0000    1.0000    1.0000
    1.0000    1.0000    1.0000
    1.0000    1.0000    0.8333
    1.0000    1.0000    0.6667
    1.0000    1.0000    0.5000
    1.0000    1.0000    0.3333
    1.0000    1.0000    0.1667
    1.0000    1.0000         0
    1.0000    0.9375         0
    1.0000    0.8750         0
    1.0000    0.8125         0
    1.0000    0.7500         0
    1.0000    0.6875         0
    1.0000    0.6250         0
    1.0000    0.5625         0
    1.0000    0.5000         0
    1.0000    0.4375         0
    1.0000    0.3750         0
    1.0000    0.3125         0
    1.0000    0.2500         0
    1.0000    0.1875         0
    1.0000    0.1250         0
    1.0000    0.0625         0
    1.0000         0         0
    0.9375         0         0
    0.8750         0         0
    0.8125         0         0
    0.7500         0         0
    0.6875         0         0
    0.6250         0         0
    0.5625         0         0
    0.5000         0         0];

%figure,
%imagesc(spectra.ppm2,spectra.ppm1,spectra.real);
imagesc(spectra.ppm1,spectra.ppm2,spectra.real);
caxis([-thresh*std(std(spectra.real)),thresh*std(std(spectra.real))]);
set(gcf,'Colormap',cmap)
set(gca,'XDir','rev')

elseif disp(1)=='c'

range=3;
levels=5;

figure
vector=(2.^[-1*range:(range-(-1*range))/(levels-1):range])*(thresh*std(std(spectra.real)));
h=contour(spectra.ppm1,spectra.ppm2,spectra.real,vector,'EdgeColor','k');
set(gca,'XDir','rev')
set(gca,'YDir','rev')
end