function [corrcoefs,p]=statsy(X1,X2,shifts,ppm1,ppm2,method)

%function [corrcoefs,p]=statsy(X1,X2,shifts,ppm1,ppm2,method)
% 
% Statistical Spectroscopy - Calculates correlation coefficients between
% features in spectral matrix 'X' and feature 'target'.  'Target' can be
% either an index of a feature in X, or a vector or features of size (X,1).
% When applied to NMR data, this is STOCSY.  When applied between NMR and
% MS, this is SHY.
%
% Arguments:
% 
% target       Chemical shift of driver peak or response vector with length
%              equal to size(X,1) 
% X            Data matrix of spectra  1=norm 2= driver
% method       Correlation metric, 'Pearson', 'Spearman', or 'Kendall'.
%              Default 'Pearson'
%
for i=1:length(shifts)
[~,index(i)]=min(abs(ppm2-shifts(i)));
end


% if length(target)==1
%     target_vect=X(:,target);
% else
%     target_vect=target;
% end
    target_vect=X2(:,index);
if size(target_vect,1)~=size(X1,1) && size(target_vect,2)==size(X1,1)
    target_vect=target_vect';
end
if exist('method')~=1
    method='Pearson';
end

for i=1:size(target_vect,2)
[corrcoefs(:,i),p(:,i)] = corr(X1,target_vect(:,i),'type',method);
end

VisLoadings1D(X1,corrcoefs',ppm1);