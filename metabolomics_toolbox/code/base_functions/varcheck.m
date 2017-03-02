function varcheck(X,features)

% varcheck(X,features)
%
% Displays histogram and box plots of log-fold change vs. median for all
% features in each spectrum.  Dilution / normalization effects are often
% visible as distributions not centered at 0.  
%
% Arguments:
% 
% X            N x P matrix of spectral features for each sample
% features     1 x N vector of spectral features - chemical shifts or m/z

if exist('features')~=1
    features=1:size(X,2);
end

variance_X=var(X);
[h,k]=sort(mean(X));

figure, plot(1:length(features),variance_X(k))
title('Variance by feature rank')
xlabel('Feature Rank')
ylabel('Variance')

[h2,k2]=sort(features);
figure, plot(h2,variance_X(k2))
title('Variance by feature')
xlabel('Feature - chemical shift or m/z')
ylabel('Variance')

logvar=log10(variance_X);
x=0:0.1:max(logvar);
n=histc(logvar',x);
figure, plot(x,n)
title('Histogram of variances')
xlabel('Log 10 Variance')
ylabel('Features (count)')