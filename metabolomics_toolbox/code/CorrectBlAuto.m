function XNew=CorrectBlAuto(X)

% Correct the basline for the spectra.
% Call [A,S] = showBaseline(ppm, X) to get A and S
% Arguements:
% X  N x P matrix of spectral features for each sample
% Outputs:
% XNew  N x P baseline-corrected spectra
%
% Bing Wang

B = 1.25;
XNew=zeros(size(X));
for i = 1:size(X,1)
    [A,S] = showBaseline(X(i,:));
    bd=Baseline(X(i,:)',A,B,S);
    XNew(i,:)=X(i,:)-bd';
end
