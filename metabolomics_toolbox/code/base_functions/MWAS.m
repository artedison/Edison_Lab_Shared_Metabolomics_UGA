function [p,sigs]=MWAS(X,Y,correction)

% [p,sigs]=MWAS(X,Y,correction)
% 
% Metabolome-wide association study design; t-test or logistic regression
% of spectral features against classes Y.  Correction for multiple
% hypothesis testing using either Sidak or Bonferroni correction.  
%
% Arguments:
% 
% X            Data matrix of spectra 
% Y            Vector of classes (0 and 1)
% correction   Either Bonferroni, 'bonferroni', or Sidak, 'sidak'
%
% Ouputs:
% p            P values for t-test or logistic regression model
% sigs         logical vector of significant features

if exist('correction')~=1
    beta=0.05/size(X,2);
elseif strcmp(correction,'bonferroni')==1
    beta=0.05/size(X,2);
elseif strcmp(correction,'sidak')==1
    beta=1-((1-0.05)^(1/size(X,2)));
else
    Error('Correction must be either bonferroni or sidak')
end

pos=find(Y==1);
neg=find(Y==0);
for k=1:size(X,2)
    [sigs(k),p(k)] = ttest2(X(neg,k),X(pos,k),beta);
end
