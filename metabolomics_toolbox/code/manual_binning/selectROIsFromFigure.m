function [x] = selectROIsFromFigure()   

%%
%{
    This function allows you to select multiple regions in your
    spectrum/spectra currently displayed and store their ppm boundaries in a 
    vector for easy manipulation. Run the output directly through matchPPMs() 
    to get ppm vector indices, if desired. 

To Run:
    Open a figure with a single spectrum or multiple spectra:
    call the function: 
    [x,y] = selectROIsFromFigure();

Return Values:
    x is the range of ppms selected
    y is usually not useful, but is the range of y values selected. 
        * better to calculate this as a max for a given spectrum, such as
        in inspectANDremovePeaks.m

    The active figure is used, so open or click inside that before running this. 
    Follow the pop-up boxes. Zoom into region of interest before clicking the 
    'Click to draw rectangle around ROI.' button. You can't just hit Return 
    to submit a response unless you hit Tab first. It is easy to get used to 
    doing that quickly. 
    
    Hope this helps,
    MJ 3MAR2017

    MJ 1MAY2017 Update: removed y values from output. 
    MJ 17MAY2017 Update: Added "Save" and "Continue" buttons.
%}
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
%%
    x = [];
    y = [];
    %listOfFigures = {};
            response = 'Y';
            i = 1;
            fprintf('\n\n');
                while ~strcmp(response, 'N')%or(strcmp(response, 'Y'),or(strcmp(response, 'y'),or(strcmp(response, 'yes'),strcmp(response, 'yea boi')))) % keep going until response in dlgbox is not one of these
            %{
                if strcmp(response,'S')            
                %% Save the current figure and regenerate
                    % Get current figure zoom level

                     % Save the figure

                    filename = ['FeatureSelection_stoppedAtFeature_',num2str(i),'.fig'];          
                        savefig(gcf,filename,'compact')
                        fprintf([filename,' was saved in "',cd(),'" ...\n']);
                    listOfFigures = vertcat(listOfFigures,filename);
                    close(gcf)
                    open(filename)
                    % Get rid of heavy patch objects except last one
                        D=get(gcf,'Children'); %get the handle of the line object
                        oldRecs = findall(D,'Type','Patch');  % make a list of all patch boxes in the figure
                        oldRecs = oldRecs(2:end);
                        delete(oldRecs)
                 end
                %}
                    
                    % Ask if ROI bound selection is wanted at the moment
                    uicontrol('Position',[2 4 300 30],'String','Click to draw rectangle around ROI.',...
                              'Callback','uiresume(gcbf)');
 
                    uiwait(gcf);% Wait for input, continue code when button is clicked
 
                    % Get a pair of points from the user
                    % [xpos,ypos] = ginput();
                    % [xmin ymin width height]
                    try % in case no region is selected
                        rect = getrect; % must drag left->right

                        xpos = [rect(1),rect(1)+rect(3)];
                        ypos = [rect(2),rect(4)];
                            % Sort ppms (with intensities) of chosen point pair so
                            % that they are low->high in value
                                xtemp = [];
                                ytemp = [];
                                xtemp = xpos';
                                ytemp = ypos';
                                    if xtemp(2) < xtemp(1)
                                        xtemp = flip(xtemp);
                                        ytemp = flip(ytemp);
                                    end

                            x = [x,xtemp]; %% then store in variable outside loop for later
                            y = [y,ytemp];

                        %% Draw a box around each region
                                llowerpt = [x(1,i),y(1,i)];
                                width = rect(3);
                                height = rect(4); %max height on spectrum
                                llx = x(1,i); % 
                                lux = llx;
                                rlx = llx + width; % = rux
                                rux = rlx;
                                lly = y(1,i);
                                luy = height;
                                rly = lly;
                                ruy = luy;
                                xcoords = [llx rlx rux lux];
                                ycoords = [lly rly ruy luy];
                                %plot(llowerpt(1),llowerpt(2),'r*')
                                %plot(llowerpt(1)+width,height,'r*')
                                %rectangle('Position',[llowerpt,width,height]);
                                p=patch(xcoords,ycoords,'r'); %light red % pink [1    0.5  1]
                                set(p,'FaceAlpha',0.1); % transparency
                    catch   % if no region is recorded
                        fprintf('\nWarning: No feature was selected. Follow prompts to continue or save and exit.\n\n');
                    end
                %% Do this again?
                    prompt = 'Enter Another ROI? (tab and return to select another ROI; ''N'' to save and exit: ';
                    response = inputdlg(prompt);
                    i = i + 1; % keep track of which region we're on

                end
            %x = x'; % this is easier to look at
            %y = y'; % ditto
            mTextBox = uicontrol('style','text','Position',[400 4 100 30]); % let us know the regions were selected
            set(mTextBox,'String',[num2str(size(x,2)),' region(s) selected']);
    hold off
    %Save again
    %{
                    filename = ['FeatureSelection_stoppedAtFeature_',num2str(i),'.fig'];          
                        savefig(gcf,filename,'compact')
                        fprintf([filename,' was saved in "',cd(),'" ...\n']);
                    listOfFigures = vertcat(listOfFigures,filename);
    %}
    %close(f);
    % Clean up (Matlab should do this anyways if it is a function)
        clearvars('matrix','maxHeight','i','xpos','ypos','xtemp','ytemp','prompt','mTextBox','llowerpt','height','width','response','llx','lux','rlx','rux','lly','luy','rly','ruy','xcoords','ycoords','p')
        
        %% Pick peaks in the data using data from figure
        %{
        h = gcf; %current figure handle
            axesObjs = get(h, 'Children');  %axes handles
            dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
            objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
            xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
                ppm = xdata;
                clear('xdata')
            X = get(dataObjs, 'YData');
        [peaks] = basicPeakPickingROIs(X,ppm,features)
        %}
        
        end