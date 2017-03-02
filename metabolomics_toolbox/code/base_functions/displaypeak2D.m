function displaypeak2D(X,XNoise,label,indices,ppm1,ppm2)

%displaypeak2D(X,XNoise,label,targetindices,ppm1,ppm2)
%
%Inputs: 
% X            Data matrix of spectra
% XNoise       Local Noise matrix
% label        label matrix from segment2D.m 
% indices      index in label matrix of peak to display - eg. [1] or
%              find(sigs==1)
% ppm1         Chemical shift vector of F2
% ppm2         Chemical shift vector of F1
% 
%
% displaypeak2D plots a set of bin regions on the mean pseudo-spectrum of X.
% Use this script to assign correlated bin regions.

spect=mean(X,3);
noise=mean(XNoise,3);

labelcopy=zeros(size(label));
for k=1:length(indices)
    labelcopy(label==indices(k))=indices(k);
end

thresh=1;
levels=10;
range=3;
[FXlab,FYlab]=gradient(labelcopy);
figure, imagesc(ppm1,ppm2,FXlab+FYlab)
cmapregions=[0.847058832645416,0.160784319043159,0;1,1,1;0.847058832645416,0.160784319043159,0];
colormap(cmapregions);
vector=(2.^[-1*range:(range-(-1*range))/(levels-1):range])*(thresh*std(std(spect./noise)));
caxis([-1 1])
hold on
h=contour(ppm1,ppm2,spect./noise,vector,'EdgeColor','k');
set(gca,'XDir','rev')
set(gca,'YDir','rev')

