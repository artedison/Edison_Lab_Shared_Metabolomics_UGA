function [X] = create_pseudo_peak(mean, sigma, intensities)

% [X] = create_gaussian_curve_Activity(x,0, 1,intensities)
% Description:
%       Create gaussian curve for each input intensity
% Input:
%       mean:         Mean for normal distribution, affect peak center
%       sigma:        Sigma for normal distribution, affect peak width
%       intensities:  Intensities related to any given activity 
x=-100:100;
    X = 1/sqrt(2*pi)/sigma*exp(-(x-mean).^2/2/sigma/sigma);
   X = X * intensities; 
end