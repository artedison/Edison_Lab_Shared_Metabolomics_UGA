function [ROIs,figureFileName] = manual_feature_selection_1()
%{
    %% Get the data from the figure
    h = gcf; %current figure handle
                axesObjs = get(h, 'Children');  %axes handles
                dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
                %objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
                xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
                    ppm = xdata{1,:};
                    clear('xdata')
                matrix = nan(size(get(dataObjs, 'YData'),1),length(ppm));
                    mat = get(dataObjs, 'YData');
                    for i = 1:size(get(dataObjs, 'YData'),1)
                        matrix(i,:) = mat{i};
                    end
%}

%% Manually select features using selectROIsFromFigure()
    timestamp = num2str(now());
     
%while ~strcmp(response, 'exit')
    % Save the current figure and regenerate
    %{
        % Get rid of heavy patch objects except last one
            D=get(gcf,'Children'); %get the handle of the line object
            oldRecs = findall(D,'Type','Patch');  % make a list of all patch boxes in the figure
            oldRecs = oldRecs(2:end);
            delete(oldRecs)
    %}
 
    % Pick each region manually

        ROIs = selectROIsFromFigure();
        %y = horzcat(y,selectROIsFromFigure()); % picked from 10-0ppm in sections. 
                                    % Gets really heavy after ~50 regions. 
                                    % I like to pick 50 or so, choose
                                    % 'S' to at prompt to save and re-open,
                                    % in order to continue picking. If you
                                    % save the figures, you can always
                                    % get the ROIs from the patch boxes
                                    % that save with it (the pink
                                    % boxes) using the code below.
%{
   % Save and continue?
        prompt = '''S'' to save and continue, or ''exit'' to save and exit (Default: save and continue): ';
        response = inputdlg(prompt);
        %i = i + 1; % keep track of which region we're on
        fprintf([num2str(size(ROIs,2)),' regions were picked so far\n\n'])
     %}
     %Save again
        figureFileName = [timestamp,'_FeatureSelection','.fig'];          
            %savefig(gcf,figureFileName,'compact')
            savefig(gcf,figureFileName)
            fprintf([figureFileName,' was saved in "',cd(),'" ...\n']);
            fprintf(['The last feature selected was at ',num2str(ROIs(1,end)),'-',num2str(ROIs(2,end)),'ppm.\n\n'])
            
        %open(filename)
%end
    %features = gatherROIsfromFigures(listOfFigures,matrix,ppm);

end
