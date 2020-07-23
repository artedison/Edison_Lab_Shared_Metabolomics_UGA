function [XN,factors,params]=normalize(X,ppmH,method,features)

    % Author: Edison Lab
    % Version: 0.1
    % Tested on Matlab Version R2017b
    % Date: 25MAR2019
    % 
    % Description:
    %       Normalize spectral matrix X to total area, single feature, or integral of set of features.        
    % 
    % Input:
    %       X        : stack of 1D spectra
    %       ppmH     : chemical shift vector
    %       method   : 'total' for Total Area, 
    %                  'PQN' for Probablistic Quotient,
    %                  'quantile' for Quantile Normalization,
    %                  'intensity' for normalization to single feature,
    %                  'integral' for normalization to sum of set of features
    %       features : only required for 'intensity' or 'integral'.  
    %                  For 'intensity', the ppm of the feature in X to normalize to - e.g. [10] for X(:,10).  
    %                  For 'integral', the range of features in X that span the peak to normalize to - eg [-.05,0.05]
    %
    %                  If method is set to 'quantile',
    %                  optionally set feature to 'median' to take median of the ranked values instead of the mean.
    %
    % Output:
    %       XN      : N x P matrix of normalized data
    %       factors : N calculated normalization factors: XN = X / factors
    %
    % Log:
%             'factors' output for PQN now includes the initial 'total' normalization factors
%             'factors' output for quantile returns nan to avoid error
%             
    %
    % Example run:
    % [XN,factors]=normalize(X,ppmR,'PQN');
params = reportParams('exclude',{'X','ppmH'});
XN=zeros(size(X));

if nargin<2
    error('incorrect number of input parameters');
end

[obs dim]=size(X);

switch lower(method)
    case 'total'
        factors=repmat(NaN,[1 obs]);
        for i=1:obs
            factors(i)=sum(X(i,:));
            XN(i,:)=X(i,:)./factors(i);
        end
    case 'pqn'
        [X,factorsT] = normalize(X,ppmH,'total');      % MTJ edit 20APR2020
        X = X.*100;                                    % MTJ edit 20APR2020
        X(0==X)=0.00000001;
        normRef=repmat(median(X),size(X,1),1);
        F=X./(normRef);
        for i=1:obs
            factors(i)=median(F(i,:));
            XN(i,:)=X(i,:)./factors(i);
        end
        factors = factors .* factorsT;                 % MTJ edit 20APR2020
    case 'intensity'
        [~,feature]=min(abs(ppmH-features));
        factors=X(:,feature);
        for i=1:obs
            XN(i,:)=X(i,:)./factors(i);
        end
    case 'integral'
        [~,feature(1)]=min(abs(ppmH-features(1)));
        [~,feature(2)]=min(abs(ppmH-features(2)));
        for i=1:obs
            factors(i)=trapz(X(i,feature(1):feature(2)));
            XN(i,:)=X(i,:)./factors(i);
        end
    case 'quantile'
        if exist('features','var');
            XN = transpose(quantilenorm(X',features,'true'));
        else
            XN = transpose(quantilenorm(X'));
        end           
        factors = nan;                                  % MTJ edit 20APR2020
end
end