function output_txt = highlightLineCallback(obj,event_obj,names)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

%In a plot with a bunch of lines, figure out which line I'm clicking
%GSS

lineHandle=get(event_obj, 'Target');
parent=get(lineHandle,'parent');
children=get(parent,'children');
idx=((children)==lineHandle);
set(children(idx),'Linewidth',3)
set(children(~idx),'Linewidth',0.5)
output_txt=['Sample: ', names{find(flipud(idx))}];
end