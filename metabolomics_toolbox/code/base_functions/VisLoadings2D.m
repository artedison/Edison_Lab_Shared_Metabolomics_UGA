function [corr,covarX]=VisLoadings2D(X,loadings,ppm1,ppm2,disp,label)

% [corr,covarX]=VisLoadings2D(X,loadings,ppm1,ppm2,disp,label)
% 
% Plots correlation/covariance projection of 2D NMR spectra using a vector of
% loadings from PCA, PLS, or a correlation vector from univariate
% regression.  
%
% Arguments:
% 
% X            Data matrix of spectra
% loadings     Vector of loadings or correlations from binmat or full spectra- eg PCA.loadings
% ppm1         Chemical shift vector of F2
% ppm2         Chemical shift vector of F1
% disp         string specifying 'contour' plot or 'fast' plot 
% label        Optional spectral segmentation matrix if using binning.
%              Length of loadings vector must equal number of segments.
%
% Return Values:
% corr         Full spectral size correlation matrix
% covarX       Full size back-projected covariance matrix.  


cmap=jet(20);

if exist('label')==0
    for k=1:size(X,3)
        linesX(k,:)=reshape(X(:,:,k),1,[]);
    end
    precorr=loadings./max(abs(loadings));
    precovarX=loadings.*std(linesX);
    corr=reshape(precorr,size(X,1),[]);
    covarX=reshape(precovarX,size(X,1),[]);
else  
    corr=zeros(size(X,1),size(X,2));
    for k=1:max(max(label))
        corr(find(label==k))=loadings(k)/max(abs(loadings));
    end
    covarX=std(X,[],3).*corr;
end

if disp(1)=='c'
    
    ind=1;
    for z=-1:1.99/(size(cmap,1)-1):.99;
        proj=NaN(size(X,1),size(X,2));
        proj(find(corr>z))=abs(covarX(find(corr>z)));
        project(:,:,ind)=proj;
        ind=ind+1;
    end
    
    figure
    range=5;
    thresh=10;
    levels=10;
    vector=(2.^[-1*range:(range-(-1*range))/(levels-1):range])*(thresh*std(std(covarX)));
    contour(ppm1,ppm2,project(:,:,1),vector,'EdgeColor',cmap(1,:))
    hold on
    for z=2:size(project,3)
        contour(ppm1,ppm2,project(:,:,z),vector,'EdgeColor',cmap(z,:));
    end
    set(gca,'YDir','rev')
    set(gca,'XDir','rev')
    xlabel('F2 (ppm)')
    ylabel('F1 (ppm)');
    colormap(jet);
    t=colorbar;
set(get(t,'ylabel'),'String', 'Correlations/Loadings');
caxis([min(loadings) max(loadings)])
set(gca,'XDir','rev');
    
elseif disp(1)=='f'
    figure, imagesc(ppm1,ppm2,covarX)
    set(gca,'XDir','rev')
    set(gca,'YDir','rev')
    xlabel('F2 (ppm)')
    ylabel('F1 (ppm)');
else
    error('disp must be take string "contour" or "fast"')
end
