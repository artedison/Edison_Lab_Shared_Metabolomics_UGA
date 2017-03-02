function whichLine()
%GSS
dcm = datacursormode(gcf);
datacursormode on;
set(dcm,'UpdateFcn',@whichLineCallback);
end