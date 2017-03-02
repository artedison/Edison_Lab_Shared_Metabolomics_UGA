function [covar]=getcovar(X,loadings,features)

[feature_sort,order]=sort(features);

corr=loadings(order);
covar=loadings(order).*std(X(:,order));
