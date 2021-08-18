function selectLine(names,par1,par2)
%Example:
%figure,plot(rand(10,4))
%highlightLine({'a','b','c','d'})
%Copyright 2014, Gregory Stupp Productions
% add return function by yue wu
global clickedRidges;
global lineNumber;
%     eval(sprintf('global %s;'),par1)
%     eval(sprintf('global %s;'),par2)
dcm = datacursormode(gcf);
datacursormode on;
set(dcm,'UpdateFcn',{@highlightLineCallback2,names});
%set(dcm,'UpdateFcn',{@highlightLineCallback2,names,par1,par2});
end
