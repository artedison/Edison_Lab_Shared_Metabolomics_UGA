function XNew=CorrectBl(X, A, S, R)

% Correct the basline for the spectra.
% Call [A,S] = showBaseline(ppm, X) to get A and S
% Arguements:
% X  N x P matrix of spectral features for each sample
% R = range of rows, only to which Baseline correction is applied
% Outputs:
% XNew  N x P baseline-corrected spectra
%
% Bing Wang

% instantiation
B = 1.25;
XNew = X;

% only one argument, use default parameters and throw warning
if nargin == 1
    warning('Using default fit paramters.')
    [A,S] = showBaseline(X);
    %set range to all rows
    R = (1:size(X,1));
end

% no R supplied
if nargin == 3
    %set range to all rows
    R = (1:size(X,1));
end

% convert R to column list
R = R';

% iterate through R and alter the rows specified
for i = 1:size(R,1)
    % select row specified in R
    Row = R(i);
    % call baseline on row specified
    bd = (Baseline(X(Row,:)',A,B,S))';
    % change respective row in XNew
    XNew(Row,:) = X(Row,:)-bd;    
end

    