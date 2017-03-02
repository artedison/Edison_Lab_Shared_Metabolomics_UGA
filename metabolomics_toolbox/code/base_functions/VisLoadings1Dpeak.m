function VisLoadings1Dpeak(XALSN,loadings,features,range)

% VisLoadings1D(X,loadings,features,range)
% 
% Correlation/covariance plot of spectra modeled by loadings or
% correlations.  
%
% Arguments:
% 
% X                (N x P) Data matrix of spectra 
% loadings         (1 x P) Vector of coefficients from PLS betas, PCA loadings, or
%                  correlations
% features         (1 x P) Vector of features (chemical shifts or m/z)
% range            Optional: range of coefficients , default [min max]

[feature_sort,order]=sort(features);

corr=loadings(order);
covar=loadings(order).*std(XALSN(:,order));

if exist('range')~=1
    range=[-1*max(abs(corr)),max(abs(corr))];
end

cmap=jet;

lines=NaN(size(cmap,1),size(corr,2));
ind=1;
for k=range(1):(range(2)-range(1))/size(cmap,1):range(2)
    lines(ind,(corr>k))=covar((corr>k));
    ind=ind+1;
end


hmult=bar(feature_sort,lines(1,:),'FaceColor',cmap(1,:));
set(hmult,'LineWidth',1);
hold on
for k=2:size(cmap,1)
    hmult(k)=bar(feature_sort,lines(k,:),'FaceColor',cmap(k,:),'EdgeColor',cmap(k,:)); 
end
%xlim([-0.5,9]);
set(hmult,'LineWidth',2);
set(gca,'XDir','rev')

xlabel('ppm')
ylabel(['Back projected UV scaled loading intensities'])


t=colorbar;
set(get(t,'ylabel'),'String', 'UV scaled loading coefficients');
caxis([-1 1])
end


