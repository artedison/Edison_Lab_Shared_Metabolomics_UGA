function [FWHM] = findFWHM(X,ppm,shift1,shift2)
%findFWHM: Finds the full width half max (FWHM) of a function

 %Optional arguments:
%     ppmrange - [-0.2, 0.2] (default)
%
%	FWHM = findFWHM(x, fx);

% if shift1 || isempty
     shift1= -0.195
     shift2= 0.2 
%i=X(1:5,:)
for i=1
    [XT,ppmR]=remove_ends2(XAL,ppm,shift1,shift2);
    [m, n] = max(XT(i,:));		%	Find maximum value and index
   %FWHM = interp1(fx(end:-1:n), x(end:-1:n), m/2, 'spline') - interp1(fx(1:n), x(1:n), m/2, 'spline');
    ind = find(XT(i,:)>=m(i)/2);	%	Find indicies where I>=max(I)/2
    nl = min(ind(i,:));			%	Leftmost index
    nr = max(ind(i,:));			%	Rightmost index

%	Linear interpolate x positions
    xl = ((X(i,:))*(nl(i))-((X(i,:))*(nl(i)-1))*(m(i)/2-(XR(i,:))*(nl(i)-1))/((XR(i,:))*(nl(i))-(XR(i,:))*(nl(i)-1)) + (X(i,:))*(nl(i)-1));
    xr = ((X(i,:))*(nr)-(X(i,:))*(nr-1))*(m/2-(XR(i,:))(nr-1))/((XR(i,:))(nr)-(XR(i,:))(nr-1)) + (X(i,:))*(nr-1);
end 
%	Get FWHM

FWHM = abs(xr-xl);






x = data(:,1);
y= data(:,2);
maxy = max(y); 
f = find(y==maxy); 
cp = x(f);% ignore Matlabs suggestion to fix!!!
y1= y./maxy;
ydatawr(:,1) = y1;
ydatawr(:,2) = x;
newFit1=find(x>= cp);
newFit2=find(x < cp);
ydatawr2 = ydatawr(min(newFit1):max(newFit1),:);
ydatawr3 = ydatawr(min(newFit2):max(newFit2),:);
sp1 = spline(ydatawr2(:,1),ydatawr2(:,2),0.5);
sp2 = spline(ydatawr3(:,1),ydatawr3(:,2),0.5);
Fullw = sp1-sp2;
