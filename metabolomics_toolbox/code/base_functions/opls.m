function PLS=opls(X,Y,nOSC)
%+++ Function: orthogonal to latent sturctures(O-PLS)
%+++ Input:  X: m x n  (Sample matrix)
%            Y: m x 1  (measured property)
%         nOSC: The number of OSC components to remove
%+++ Reference: Johan Trygg and Svante Wold, Orthogonal projections to latent
%    structures (O-PLS),J. Chemometrics 2002; 16: 119-128
%+++ Author: by H.D. Li, Fool's day, 2012.
%+++ Contact:lhdcsu@gmail.com
%    Modified May 2014 Greg Stupp

if nargin<3;nOSC=1;end

[n,p]=size(X);
Xorig=X;
ssqXcal=sum(sum((X.^2))); 
w=X'*Y;         % calculate weight vector
w=w/norm(w);    % normalization
for i=1:nOSC
    prevssq = sum(sum(X.^2));
    t=X*w;                     % calculate scores vector
    p=X'*t/(t'*t);             % calculate loadings of X
    wosc=p-(w'*p)/(w'*w)*w;    % orthogonal weight
    wosc=wosc/norm(wosc);      % normalization
    tosc=X*wosc;               % orthogonal components
    posc=X'*tosc/(tosc'*tosc); % loadings 
    
    X=X-tosc*posc';            % remove orthogonal components
    W_orth(:,i)=wosc;          % record results
    scores(:,i)=tosc;          % record results
    loadings(i,:)=posc';          % record results
    variance(i)=1-sum(sum(X.^2))/prevssq; % ratio of explained variance of OSC components
end

Zcal=X;   %+++ data with orthogonal signal components removed

%+++ Output
PLS.W=W_orth;
PLS.scores=scores;
PLS.loadings=loadings;
PLS.Zcal=Zcal;
PLS.variance=variance;

%
% beta = PLS.W*PLS.loadings';
% beta = [meanY - meanX*beta; beta];
