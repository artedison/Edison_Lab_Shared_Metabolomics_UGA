function FWHM = findFWHM(x, fx)
%findFWHM: Finds the full width half max (FWHM) of a function
%
%	FWHM = findFWHM(x, fx);


[m, n] = max(fx);		%	Find maximum value and index
%FWHM = interp1(fx(end:-1:n), x(end:-1:n), m/2, 'spline') - interp1(fx(1:n), x(1:n), m/2, 'spline');
ind = find(fx>=m/2);	%	Find indicies where I>=max(I)/2
nl = min(ind);			%	Leftmost index
nr = max(ind);			%	Rightmost index

%	Linear interpolate x positions
xl = (x(nl)-x(nl-1))*(m/2-fx(nl-1))/(fx(nl)-fx(nl-1)) + x(nl-1);
xr = (x(nr)-x(nr-1))*(m/2-fx(nr-1))/(fx(nr)-fx(nr-1)) + x(nr-1);

%	Get FWHM
FWHM = abs(xr-xl);