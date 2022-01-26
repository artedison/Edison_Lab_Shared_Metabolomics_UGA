function [hax,hlines] = plotr(varargin)
%Plot with x-axis reversed by default. Any arguments get passed into 'plot'
h = plot(varargin{:});
set(get(h(1),'Parent'),'xdir','rev')
set(gcf,'color','w');
box off
if nargout > 0
    hax=get(h(1),'Parent');
end
if nargout ==2
    hlines=h;
end