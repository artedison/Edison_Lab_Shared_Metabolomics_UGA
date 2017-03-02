function [corr,covarX]=VisLoadings2DST(X1,X2,ppm1,ppm2,disp,label)

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

% if exist('label')==0
%     for k=1:size(X,3)
%         linesX(k,:)=reshape(X(:,:,k),1,[]);
%     end
%     precorr=loadings./max(abs(loadings));
%     precovarX=loadings.*std(linesX);
%     corr=reshape(precorr,size(X,1),[]);
%     covarX=reshape(precovarX,size(X,1),[]);
% else  
%     corr=zeros(size(X,1),size(X,2));
%     for k=1:max(max(label))
%         corr(find(label==k))=loadings(k)/max(abs(loadings));
%     end
%     covarX=std(X,[],3).*corr;

target_vect=X2;
  
if size(target_vect,1)~=size(X1,1) && size(target_vect,2)==size(X1,1)
    target_vect=target_vect';
end

corr=zscore(target_vect)'*zscore(X1);
corr=corr';
covarX=(X1-repmat(mean(X1),size(X1,1),1))'*(X2-repmat(mean(X2),size(X2,1),1));
rem=abs(covarX)<1;
rem1=sum(rem,1)==8192;
rem2=sum(rem,2)==4096;
ppm12r=ppm12;
ppm12r(rem2)=[];
ppm22r=ppm22;
ppm22r(rem1)=[];

covarX(:,rem1)=[];
covarX(rem2,:)=[];
corr(:,rem1)=[];
corr(rem2,:)=[];

% end
project=[];
% if disp(1)=='c'
    
    ind=1;
    for z=-1:1.99/(size(cmap,1)-1):.99;
    %for z=-1:.5:1 
        proj=NaN(size(corr));
        proj(corr>z)=abs(covarX(corr>z));
        project(:,:,ind)=proj;
        ind=ind+1;
    end
    
    
    figure
    range=5;
    thresh=10;
    levels=10;
    vector=(2.^[-1*range:(range-(-1*range))/(levels-1):range])*(thresh*std(std(covarX)));
    contour(ppm22r',ppm12r,project(:,:,1),vector,'EdgeColor',cmap(1,:))
    hold on
    for z=2:size(project,3)
        contour(ppm22r,ppm12r,project(:,:,z),vector,'EdgeColor',cmap(z,:));
    end
    set(gca,'YDir','rev')
    set(gca,'XDir','rev')
    xlabel('F2 (ppm1)')
    ylabel('F1 (ppm2)');
    
elseif disp(1)=='f'
    figure, imagesc(ppm1,pm2,corr)
    set(gca,'XDir','rev')
    set(gca,'YDir','rev')
    xlabel('F2 (ppm1)')
    ylabel('F1 (ppm2)');
else
    error('disp must be take string "contour" or "fast"')
end
