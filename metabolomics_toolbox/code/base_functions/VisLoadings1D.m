function VisLoadings1D(XALSN,loadings,features,range)

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

cmap=jet(100);

lines=NaN(size(cmap,1),size(corr,2));
ind=1;
for k=range(1):(range(2)-range(1))/size(cmap,1):range(2)
    lines(ind,find(corr>k))=covar(find(corr>k));
    ind=ind+1;
end

figure
%subplot(3,2,6); 
plot(feature_sort,lines(1,:),'Color',cmap(1,:))
hold on
for k=2:size(cmap,1)
    %subplot(3,2,6); 
    plot(feature_sort,lines(k,:),'Color',cmap(k,:));
end
xlabel('Chemical Shift')
ylabel(['Loadings coefficients * std(data)'])

colormap(cmap);
t=colorbar;
set(get(t,'ylabel'),'String', 'Loadings Coefficients');
caxis([-1 1])
set(gca,'XDir','rev');
hold off

end


