function PCA=nipalsPCA(X,components)

% PCA=nipalsPCA(X,components)
%
% Calculates first 'components' principal components of X by NIPALS algorithm and
% outputs an array with scores, loadings, and variance for each principal
% component
%
% Arguments:
% 
% X               N x P matrix of spectral features for each sample
% components      scalar, number of components to calculate
%
% Outputs:
% PCA.scores      Scores for each spectrum in X along each principal component
% PCA.loadings    Coefficients of each variable contributing to generation of the scores in each principal component.
% PCA.variance    Total variance accounted for by each principal component

% Xi=zscore(X);
Xi=X;
T=zeros(size(X,1),components); %T will be PCA scores
P=zeros(size(X,2),components); %P will be PCA loadings
Var=zeros(1,components);

thresh=1e-15;

for i=1:components
    residual=1;
    p_initial=zeros(size(X,2),1); %initalize loadings as zero vector
    t=Xi(:,1); %initialize scores as a column of Xi
    
    while residual > thresh
        p=Xi'*t/(t'*t); %Project Xi onto t to find the corresponding loading p
        p=p/sqrt(p'*p); %Normalise loading vector p to length 1
        t=Xi*p/(p'*p);  %Project X onto p to find corresponding score vector t
        E=(p_initial-p);
        p_initial=p;
        residual=E'*E;
    end
    
    T(:,i)=t;
    P(:,i)=p;
    Xi=Xi-t*p'; %Remove the estimated PC component from Xi
    Var(i)=(sum(sum((t*p').*(t*p')))/sum(sum(X.*X)));
end

PCA.scores=T;
PCA.loadings=P';
PCA.variance=Var;

