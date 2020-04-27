function [A,B] = PeakBounds(cursor_info) 
%
% Determine the peak boundries that will be used within the 'integrate'
% function, based upon the cursor_info data selected by the user.
% 
% Arguments:
%       cursor_info    
%               The cursor info table generated after exporting selected points between which you wish to integrate
% Output:
%       A- Array of all the left (or right) peak boundaries selected for
%       each peak of interest
%       B- Array of all the right (or left) peak boundaries selected for
%       each peak of interest    
% 
[x,y]=size(cursor_info);
A = zeros(1,y/2);
B= zeros(1,y/2);
counter = 1;
for i = y:-2:2
    X = cursor_info(1,i).Position(1);
    A(counter) = X;
    counter = counter + 1;
end
counter = 1;
for i = y-1:-2:1
    X = cursor_info(1,i).Position(1);
    B(counter) = X;
    counter = counter + 1;
end