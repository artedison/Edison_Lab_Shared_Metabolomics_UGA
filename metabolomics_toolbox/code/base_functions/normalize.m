function [XN,factors]=normalize(X,ppmH,method,features)

% [XN,factors]=normalize(X,method,features)
%
% Normalize spectral matrix X to total area, single feature, or integral of set of features. 
%
% Arguments:
% 
% X               N x P matrix of spectral features for each sample
% method          'total' for Total Area, 'PQN' for Probablistic Quotient,
%                 'intensity' for normalization to single feature,
%                 'integral' for normalization to sum of set of features
% features        only required for 'intensity' or 'integral'.  For
%                 'intensity', the ppm of the feature in X to normalize
%                 to - e.g. [10] for X(:,10).  For 'integral', the range of
%                 features in X that span the peak to normalize to - eg
%                 [-.05,0.05]
%
% Outputs:
% XN              N x P matrix of normalized data
% factors         N calculated normalization factors: XN = X / factors


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
        X=normalize(X,ppmH,'total').*100;
        X(0==X)=0.00000001;
        normRef=repmat(median(X),size(X,1),1);
        F=X./(normRef);
        for i=1:obs
            factors(i)=median(F(i,:));
            XN(i,:)=X(i,:)./factors(i);
        end
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
        
end
end