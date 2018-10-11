function reorderFigures(order)

%% 
% Reorder figures open in Matlab so that they are all brought to the front
%   (even if minimized) and lined up back-to-back. Use Command-'~' to switch 
%   between them in order front-to back. The order is determined by the
%   order in which they were opened:
%   
%   order:  'highToLow' - list figures with latest figures up front
%           'lowToHigh' - list figures with latest figures at the back
%   
%   Hope this helps,
%   MJ 3MAR2017


    figures = findall(0,'type','figure');
    figNums = [];
    figNums = sort(vertcat(figNums,figures.Number));
    if strcmp('highToLow',order)
        for i = 1:length(figures)
           figure(figNums(i));
        end        
    else
        if strcmp('lowToHigh',order)
            for i = 1:length(figures)
               figure(figNums((length(figures)-i)+1));
            end
        end
    end
end