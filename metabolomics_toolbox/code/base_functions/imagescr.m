function imagescr(x,y,z,thresh)
if exist('thresh')==0
    thresh=5;
end
imagesc(x,y,z)
caxis([-thresh*std(std(z)),thresh*std(std(z))]);
set(gcf,'Colormap',redblue)
set(gca,'XDir','reverse')
set(gca,'YDir','reverse')
end