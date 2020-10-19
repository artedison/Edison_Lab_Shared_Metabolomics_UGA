function reorderFigures(order)
%% reorderFigures

    % Author: MTJ
    % Version: 0.1
    % Tested on Matlab Version R2020a
    % Date: JUN2020
    %
    % Description:
    %       
    %       Reorder figures open in Matlab so that they are all brought to the front
    %       (even if minimized) and lined up back-to-back. Use Command-'~' to switch 
    %       between them in order front-to back. The order is determined by the
    %       order in which they were opened:
    %
    %       order:  'highToLow' - list figures with latest figures up front
    %               'lowToHigh' - list figures with latest figures at the back
    %
    % Input:
    %
    %       order:  'highToLow' - list figures with latest figures up front
    %               'lowToHigh' - list figures with latest figures at the back
    %
    % Output:
    %       
    %       None; figure windows are simply re-ordered. 
    %
    % Log:
    %
    %       MTJ SEPT2020 - may not work with hidden or minimized figures 
    %
    % Example run:
    %
    %       reorderFigures('lowToHigh')
    %       reorderFigures('highToLow')
    %   

%%

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