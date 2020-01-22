function stackSpectra_highlightLines(timepoints,timeVector,units,extraText)
    % Provided a stackSpectra figure, a time vector corresponding to rows 
    % of the matrix, and a list of timepoints of interest, make the
    % lines closest to those timepoints bold. 
 
    % Params (for dvpt)
%         timepoints = [3,6];
%         timeVector = HSQCETGPSISP.timesCollapsed(1:10:end);
%         units = 'h';

    % Assign Current Figure handle to gca
        fig = gca;

    % Find the indices of the timepoints of interest
        timeInds = matchPPMs(timepoints,timeVector);
        
    % Find the corresponding lines
        %props = get(fig);
        children = get(fig,'Children');
        
    % Make those lines bold
        set(children(timeInds),'LineWidth',2)
        
    % Add text labels to those times
        %figure,plot(currentppm,children(timeInds(1)).YData);
        for n = 1:length(timeInds)
            xpos = children(timeInds(n)).XData(1)-0.25;
            ypos = children(timeInds(n)).YData(1);
            textLabel = num2str(timeVector(timeInds(n)));
            text(xpos,ypos,[textLabel,' ',units,'; ',extraText{n}])
        end

end