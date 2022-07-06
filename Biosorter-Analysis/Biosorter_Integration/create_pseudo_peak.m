function [X] = create_pseudo_peak(mean, sigma, intensities)

% [X] = create_pseudo_peak(0, 1,intensities)
% Description:
%       Create gaussian curve for each input intensity
% Input:
%       mean:         1-by-1 double, mean for normal distribution, affect peak center
%       sigma:        1-by-1 double, sigma for normal distribution, affect peak width
%       intensities:  1-by-1 double, intensities related to any given activity 
x=-100:100;
    X = 1/sqrt(2*pi)/sigma*exp(-(x-mean).^2/2/sigma/sigma);
   X = X * intensities; 
end