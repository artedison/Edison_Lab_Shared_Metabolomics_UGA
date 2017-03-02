function output_txt = whichLineCallback(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

%In a plot with a bunch of lines, figure out which line I'm clicking
%GSS

lineHandle=get(event_obj, 'Target');
parent=get(lineHandle,'parent');
children=get(parent,'children');
idx=find(flipud(children)==lineHandle);
output_txt=['Line: ', num2str(idx)];
end