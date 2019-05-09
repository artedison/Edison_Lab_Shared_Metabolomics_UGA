function [corr,covar,lines]=StackSTOCSYall_v3(targetlist,X,ppm,rb,lb)

% Author: Sicong Zhang, modified from function STOCSY.m
% Date: 05/07/2019
%
% StackSTOCSYall_v3(targetlist,X,ppm)
% Description:
%       Stacked STOCSY for all 78 glycans as drivers, version 3: normalize
%       on absolute covariance, only in a short chemical shift region
%
% Input: 
% targetlist   A list of driver peak
% X            Data matrix of spectra
% ppm          Chemical shift vector corresponding to X
% rb           Right boundry of NMR spectra for calculation
% lb           Left boundry of NMR spectra for calculation
%
% Output:
%       corr: Correlation matrix to target driver peak 
%       covar: Covariance matrix to target driver peak 
%       lines: Matrix of covariance, orgnized by correlation for coloring


corr=NaN(size(X));
covar=NaN(size(X));
cmap=jet(100);

for ind1=1:length(targetlist)
    [~,k]=min(abs(ppm-targetlist(ind1)));
    target_vect=X(:,k);
    
    corr(ind1,:)=(zscore(target_vect')*zscore(X))./(size(X,1)-1);
    covar(ind1,:)=(target_vect-mean(target_vect))'*(X-repmat(mean(X),size(X,1),1))./(size(X,1)-1);
end
covar=normalize(abs(covar(:,rb:lb)),ppm(rb:lb),'pqn');
corr=corr(:,rb:lb);
lines=NaN(size(cmap,1),length(rb:lb),size(targetlist,2));
for ind1=1:length(targetlist)
    ind=1;
    for k=-1:2/size(cmap,1):.99
        lines(ind,find(corr(ind1,:)>k),ind1)=covar(ind1,find(corr(ind1,:)>k));
        ind=ind+1;
    end
end
end




