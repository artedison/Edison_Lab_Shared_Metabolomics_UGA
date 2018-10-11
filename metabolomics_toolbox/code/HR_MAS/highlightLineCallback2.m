function output_txt = highlightLineCallback2(obj,event_obj,names,hObject)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

%In a plot with a bunch of lines, figure out which line I'm clicking
%GSS
% eval(sprintf('global %s;'),par1)
% eval(sprintf('global %s;'),par2)
global clickedRidges;
global lineNumber;
%%%
lineHandle=get(event_obj, 'Target');
parent=get(lineHandle,'parent');
children=get(parent,'children');
idx=((children)==lineHandle);
if get(children(idx),'Linewidth') < 9
    set(children(idx),'Linewidth',9)
else
    set(children(idx),'Linewidth',7)
end
%%return index
% ind=struct();
% ind.x=lineHandle.XData;
% ind.y=lineHandle.YData;
% ind.z=lineHandle.ZData;
%line_i_ch=strcat('line',num2str(lineNumber));
    %eval(sprintf('line_i_ch=strcat(''line'',num2str(%s));'),par1)
clickedRidges(lineNumber)=find(flipud(idx));
    %eval(sprintf('%s.(line_i_ch)=idx;'),par2)
lineNumber=lineNumber+1;

%%%%%%%%%%%%%%%%%%
output_txt=['Sample: ', names{find(flipud(idx))}]

end
