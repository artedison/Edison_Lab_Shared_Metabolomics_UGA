function [A,S] = showBaseline(X, ppm, A, S)
% Function for playing with and vizualizing baseline corrections
% To get the defaults for A and S, call with only 1 argument: [A,S] = showBaseline(X)
% Then, you can adjust A and S until you get a nice fit. showBaseline(X, ppm, A, S);
% Use these values to call XNew=CorrectBl(X, A, S)

% S: Noise Level. Raise to rise above the noise
% A: Smoothing Factor. lower = fit more closely, higher = fit less (smoother)

% GSS, May 5, 2014

if nargin == 1 % if only 1 input, use defaults
    A=5e-9*size(X,2)^4;
    S=17000;
    if S==0
        S=1;
    end
    return
end
B = 1.25;
if size(X,1)>1
    warning('Only showing 1st spectrum')
end
bd=Baseline(X(2,:)',A,B,S);

figure,h1=subplot(2,1,1);
plot(ppm,X(2,:)), hold
plot(ppm,bd,'r')
set(gca,'XDir','reverse')
h2=subplot(2,1,2);
plot(ppm,X(2,:)-bd')
set(gca,'XDir','reverse')
linkaxes([h1, h2]);