function highlightLine(names)
%Example:
%figure,plot(rand(10,4))
%highlightLine({'a','b','c','d'})
%Copyright 2014, Gregory Stupp Productions
dcm = datacursormode(gcf);
datacursormode on;
set(dcm,'UpdateFcn',{@highlightLineCallback,names});
end