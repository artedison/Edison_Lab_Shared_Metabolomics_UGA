function [multiAlignMat,multiAlign_out] = multiAlign_interactive_draft_2(X,ppm,outdir,preAligns,lastFig,patches)

%% multiAlign_interactive_draft_2 (replaces multiAlign_interactive_draft_1)
% NOTE: you should always use preAligns as input, as data extraction from
% figures is not currently working.

% This draft saved on 15NOV2017, but previous version comments kept.
% MJ 1NOV2017
% MJ 2NOV2017 - edit to handle pre-calculated alignments (dynamically)
% MJ 15NOV2017 - edit to allow for import of previous boxes from a figure
%                   supplied in lastFig parameter. If program crashes, the
%                   current figure handle for the interactive figure is
%                   returned, allowing the user to simply re-call the
%                   function to seamlessly continue picking. In the
%                   combined plot, the regions are now drawn by the
%                   regions determined by the extracted patch objects from
%                   the current combined figure rather than from ROIs,
%                   which is a temporary field that happened to work.
% MJ 16NOV2017 - bug fix with axes not loading for previous figure. Also
%                   added a 'remove last region' button so that regions can
%                   more easily be erased. Comparison figure with patches
%                   is now timestamped to eliminate overwriting. 
% MJ 4DEC2017 - bug fix for lines not being extracted from axes in the
%                   correct order. Rely only on preAligns for replacements
%                   and default alignment data instead of data in the
%                   figure. See lines 218,229. Did not change in option 4
%                   yet.
%{
% Usage:
% %% Pre-align spectra and build prealigns input structure (gives user more flexibility with alignment params)
% currentppm = ppmR;
% matrix = XR;
%     % Suppress figures:
%         set(0,'DefaultFigureVisible','off')
%     % Plot the unaligned data for comparison:
%             prealigns(1).Data = matrix;
%             prealigns(1).Titles = 'Unaligned';
%     
%     % Try the Pearson PAFFT
%         XRALg_guide_PAFFT_correl=guide_align1D(matrix,currentppm,'correlation','PAFFT');
%             % HCA on raw spectra shows expected groupings of blanks and samples!!
%         % Plot it:
%             prealigns(2).Data = XRALg_guide_PAFFT_correl;
%             prealigns(2).Titles = 'XRALg - Guide Align - PAFFT (correlation)';
%  
%    % Try the Spearman PAFFT     
%         XRALg_guide_PAFFT_spear=guide_align1D(matrix,currentppm,'spearman','PAFFT');
%             % HCA on raw spectra shows expected groupings of blanks and samples!!
%         % Plot it:
%             prealigns(3).Data = XRALg_guide_PAFFT_spear;
%             prealigns(3).Titles = 'XRALg - Guide Align - PAFFT (spearman)';
% 
%    % Try the Star Align PAFFT
%         XRALg_PAFFT_mean=star_align1D(matrix,currentppm,'mean','PAFFT');
%         % Plot it:
%             prealigns(4).Data = XRALg_PAFFT_mean;
%             prealigns(4).Titles = 'XRALg - Star Alignment - PAFFT (mean)';
% 
%     % Try the Star Align CCOW
%         XRALg_CCOW=star_align1D(matrix,currentppm,'mean','CCOW');
%         % Plot it:
%             prealigns(5).Data = XRALg_CCOW;
%             prealigns(5).Titles = 'XRALg - Star Alignment - CCOW';
%                 
%         set(0,'DefaultFigureVisible','on') % turn figures back on
% %% Compare alignments, choose region-by-region, and gather consensus XAL
% 
%     matrix = XR;
%     currentppm = ppmR;
%     outdir = 'output_multiAlign_interactive';
%     
%     [XALg,lastFig] = multiAlign_interactive_draft_2(matrix,currentppm,outdir,preAligns,lastFig,patches);
%   
%}
%% Initialize return params and directory in the event the program is terminated prematurely. Also make the multipleplot:
    multiAlignMat = [];
   
%% Check if figure is provided:

    if ~isempty(lastFig)    % if lastFig is supplied
        if ischar(lastFig)  % if lastFig is supplied as name of a saved file in path, open it:
            fprintf('\n\tOpening the supplied figure to use boxes as a starting point ...\n')
            lastFig = hgload(lastFig);
        else                % if lastFig is supplied as a figure handle, bring it to front:
            figure(lastFig)
        end  
        % Check with user to make sure this is the right figure before
        % continuing!
        inpt = menu('Are you sure the figure provided is the multiple plot figure you want patch objects from?','Yes, continue','No, restart');
        switch inpt
            case 1
                
            case 2
                fprintf('\n\n\tQuitting multiAlign_interactive because incorrect figure was used.\n')
                return
        end
        % In either event, get the axes from the figure
            axes = flipud(findall(gcf,'Type','axes')); % reads them in reverse order
            titles = {};
            for i = 1:length(axes)
                titles{i} = axes(i).Title.String;
            end
    else                % if lastFig is NOT supplied, it needs to be generated:
        % Create a directory to store all this stuff:
            fprintf('\n\tCreating output directory ...\n')

            mkdir(outdir)
            addpath(genpath(outdir))
            cd(outdir)
        % If the alignments are supplied in preAligns structure, skip the calculations. If not, calculate them:
            if isempty(preAligns)
                % Calculate some alignments to compare (see function for
                % defaults)
                    [titles] = trySomeAlignments(X,ppm,'save');            
            else
                % Generate plots for the provided spectra using their titles:
                    fprintf(['\n\tCreating plots for ',num2str(length(preAligns)),' different alignments.\n\tDo not touch the figures!...\n'])
                    titles = cell(1,1);
                    for i = 1:length(preAligns)
                            figure, plot(ppm,preAligns(i).Data)
                            set(gca,'XDir','reverse')
                            titles{i} = [preAligns(i).Titles,' (',num2str(i),')'];
                            title(titles{i})
                            xlabel('Chemical Shift (ppm)')
                            ylabel('Signal Intensity')
                            if ~isempty(patches)
                                highlightROIs(patches{i},max(preAligns(i).Data(1:numel(preAligns(i).Data)))); %*****************************
                            end
                                fig = gca;
                                saveas(fig,fig.Title.String)    
                                close(gcf);
                    end 
            end
        % In either event (alignments are supplied or not), put everything on one plot for comparison
            addpath(genpath(outdir))
                fprintf('\n\tCreating interactive figure with multiple alignments...\n')
                axes = MultiplePlot_mod(2,titles); % replacing axes field here
                                    fig = gcf;
                                    saveas(fig,'Comparison of different alignment methods.fig');          
    end
    
    % Allow the user to use the datatip to select a line, making it bolded
        highlightLine({}) % names = cellfun(@num2str,{1:size(X,1)},'UniformOutput',0)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INTERACTIVE PART:

    % Scroll through the alignments to compare each of the methods side-by-side, including the unaligned data:
    
        ROIs = [];
        realigned = struct();
        i = 1;
        try                                             % this prevents dumb errors
                nextMove = 1;
                while nextMove ~= 0                     % while you are still browsing

                    % Select new region, or close?
                        zoom on; pause(); zoom off;     % let us zoom until button press   
                        nextMove = roiMenu();           % once a key is pressed, prompt user
                        switch nextMove
                            case 1 % if select new region:
                                % Zoom, then click on desired axis:
                                    ROIs(:,end+1) = drawROI();
                        
                            case 2 % if re-draw last region:
                                % delete patch object for last box    
                                    D=get(gca,'Children'); %get the handle of the line object
                                    patchobjs = findall(D,'Type','Patch');  % make a list of all patch boxes in the current axis only
                                    delete(patchobjs(1)); 
                                % Record new ROI  
                                    ROIs(:,end) = drawROI();
                            case 3 % if delete:
                                % delete patch object for last box    
                                    D=get(gca,'Children'); %get the handle of the line object
                                    patchobjs = findall(D,'Type','Patch');  % make a list of all patch boxes in the current axis only
                                    delete(patchobjs(1)); 
                            case 4 % if Re-align using one of the default methods:
                                % Redo the alignment locally within ROI
                                    ROIs(:,end+1) = drawROI(); % record new ROI
                                    [matrix,thisppm] = getDataFromAxes(gca);
                                    [matrix,thisppm] = getDataFromAxes(gca);
                                    [regionRealigned,realigned(i).ROI,realigned(i).Method] = localAlign(matrix,thisppm,ROIs(:,end)); % redo alignment
                                        i = i + 1;
                                    %highlightROIs(ROI,max(regionRealigned(1:numel(regionRealigned))));                                % Show user the altered alignment in context
                                    
                                % Prompt to accept new alignment
                                
                                % How do we keep track of every change that was made?
                            case 5 % if done
                                nextMove = 0;
                            case 6 % do nothing
                                
                            %break    
                        end
                end
        catch
            if or(nextMove == 2,nextMove == 3)
                % delete patch object for last box (this takes care of failed
                % redraw attempt:
                    D=get(gca,'Children'); %get the handle of the line object
                    patchobjs = findall(D,'Type','Patch');  % make a list of all patch boxes in the current axis only
                    delete(patchobjs(1));
            end
            lastFig = gcf;
            saveas(lastFig,['Comparison of different alignment methods_patches_',num2str(now),'.fig']);
            msgbox(['The program crashed (likely, a key was pressed out of sequence). Your progress '...
                    'was saved in the current figure; simply re-run the function call to continue.'])
            return % this will close the program 
        end
        lastFig = gcf;
        saveas(lastFig,['Comparison of different alignment methods_patches_',num2str(now),'.fig']);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% For each axis, collect the corresponding patch objects
    
            fprintf('\n\tCollecting patch objects from muliple alignment plots...\n')

        fig = gcf;
        patches = gatherROIs_multiPlot(fig,axes); % 'axes' comes from the multiplePlot_mod function
        
    % Build the initial spectral matrix, replace patch regions
        % Assign/build default matrix
            % Prompt user to select default matrix
                defaultAxes = inputdlg('Which alignment scheme should be used as the default? (enter plot number:)');
            % Get the data from the figure *** This could be done safer by
                % passing the matrices corresponding to each axes ***
                fprintf('\n\tBuilding default matrix...\n')
                %[multiAlignMat,thisppm] = getDataFromAxes(axes(str2num(defaultAxes{:})));
                multiAlignMat = preAligns(str2num(defaultAxes{:})).Data;
            % Make note of what the default matrix is                 
                default = axes(str2num(defaultAxes{:}));
                default = default.Title.String;
    % For each axes (alignment method), replace all patch regions with the 
        % corresponding region from the data in that axes (alignment
        % method)
        for i = 1:length(axes) % same as length(patches)
            % Get the data from the figure *** This could be done safer by
                % passing the matrices corresponding to each axes ***
                %[matrix,~] = getDataFromAxes(axes(i));
                matrix = preAligns(i).Data;
            % For every patch in that axes, get the region in the axes and
            % store it over top of the same region in multiAlignMat
                        fprintf('\n\tCombining alignments...\n')

            currentPatches = patches{i};
            currentBoundsInds = matchPPMs(currentPatches,ppm);
                for j = 1:size(currentPatches,2)
                    replaceInds = currentBoundsInds(2,j):currentBoundsInds(1,j);
                    multiAlignMat(:,replaceInds) = matrix(:,replaceInds);
                    
                end
        end
                        % Plot to check
                            matrix = multiAlignMat;
                            figure, plot(ppm,matrix)
                                set(gca,'XDir','reverse')
                                default = regexprep(regexprep(default,'XRALg - ',''),' - All Samples \(\w\)','');
                                title(['XRALg - Combined, "',default,'" as the default'])
                                xlabel('Chemical Shift (ppm)')
                                ylabel('Signal Intensity')
                                %highlightROIs(ROIs,max(matrix(1:numel(matrix))));  % ROIs was recorded during picking, not directly dependent on patches
                                    for i = 1:numel(patches)
                                        if ~isempty(patches{i})
                                            highlightROIs(patches{i},max(matrix(1:numel(matrix))));  % ROIs directly dependent on patches
                                        end
                                    end
                            fig = gca;
                            saveas(fig,fig.Title.String);
           %% Make output structure:                 
            multiAlign_out = struct();
            multiAlign_out.lastFig = lastFig;
            multiAlign_out.Patches = patches;
            multiAlign_out.Realigned = realigned;
                            
            fprintf(['\n\tFigures were saved in \n\t"',cd(),'"\n\tand combined alignment was returned!...\n'])
            
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Accessory functions:

    function xtemp = drawROI() 
                        %zoom on; pause(); zoom off;  % let us zoom until button press   
                        rect = getrect(gca); % must drag left->right. gca tells getrect which plot to consider for the rectangle (default is entire figure). it would be nice to let us know which axis is selected. 
                            xpos = [rect(1),rect(1)+rect(3)];
                            %ypos = [rect(2),rect(4)];
                                % Sort ppms (with intensities) of chosen point pair so
                                % that they are low->high in value
                                    %xtemp = [];
                                    %ytemp = [];
                                    xtemp = xpos';
                                    %ytemp = ypos';
                                        if xtemp(2) < xtemp(1)
                                            xtemp = flip(xtemp);
                                            %ytemp = flip(ytemp);
                                        end

                    %% Draw the box:
    %                 height = max(data(1:numel(data))); % maximum data point
                      height = abs(rect(2)-rect(4));
                            highlightROIs(xtemp,height)
    end

    function [matrix,ppm] = getDataFromAxes(axes)
        dataObjs = get(axes, 'Children'); %handles to low-level graphics objects in axes
            objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
            xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
                xdata = cell2mat(xdata(strcmp(objTypes,'line')));
                ppm = xdata(1,:);               % x-axis values
            ydata = get(dataObjs, 'YData');
                matrix = cell2mat(ydata(strcmp(objTypes,'line'))); % these are all the lines, in the order plotted 
    end

    function [choice] = roiMenu()
        % Ask 
%             choice = listdlg('PromptString','Menu:', ...
%                                 'SelectionMode','single',...
%                                 'ListString',...
%                                 {'Select another region',...
%                                  'Recalculate Alignment',...
%                                  'Delete last region',...
%                                  'I''m done picking regions',...
%                                  'Accidental button press'}); 
        % Alternatively: 
            choice = menu('What do you want to do?',...
                     'Select another region',...
                     'Redraw last region',...
                     'Delete last region',...
                     'Recalculate Alignment',...
                     'I''m done picking regions',...
                     'Accidental button press'); 
%          % Display button
%             % Construct a questdlg with four options
%                 choice = questdlg('Select another region?', ...
%                     'Region Selection Options', ...
%                     'Yes','Redraw last region','No, I''m done picking regions','No, I''m done picking regions');
%             % Handle response
%                 switch choice
%                     case 'Yes'
%                         button = 1;
%                     case 'Redraw last region'
%                         button = 2;
%                     case 'No, I''m done picking regions'
%                         button = 0;
%                 end
    end

    function [regionRealigned,ROI,method,figname] = localAlign(X,ppm,ROI)
        % Trim X to ROI (assumes ROI is column indices small-large)
            thismatrix = X(:,matchPPMs(ROI(1),ppm):matchPPMs(ROI(2),ppm));
        % Trim ppm to ROI
            currentppm = ppm(matchPPMs(ROI(1),ppm):matchPPMs(ROI(2),ppm));
            
        %     % Suppress figures:
                 set(0,'DefaultFigureVisible','off')
            % Unaligned data:
                    realigns(1).Data = thismatrix;
                    realigns(1).Title = 'Unaligned';

            % Try the Pearson PAFFT
                    realigns(2).Data = guide_align1D(thismatrix,currentppm,'correlation','PAFFT');
                    realigns(2).Title = 'Guide Align - PAFFT (correlation)';

           % Try the Spearman PAFFT    
                    realigns(3).Data = guide_align1D(thismatrix,currentppm,'spearman','PAFFT');
                    realigns(3).Title = 'Guide Align - PAFFT (spearman)';

           % Try the Star Align PAFFT
                    realigns(4).Data = star_align1D(thismatrix,currentppm,'mean','PAFFT');
                    realigns(4).Title = 'Star Alignment - PAFFT (mean)';

            % Try the Star Align CCOW
% Seg_ppm                Length of segment in ppm (dafault 0.08)
% MaxShift_ppm           Maximum distance a segment can be shifted in ppm
%                        (default is 0.05)
% slack_ppm              Slack parameter for CCOW in ppm- default is 0.005
                    %windowSize = peakwidth/10;
                        numPeaks = length(findpeaks(movmean(mean(thismatrix,1),3)));
                        numPeaks(numPeaks<2) = 2; 
                    numSegs = numPeaks;   % must be > 2. this appears to work best if ~# of peaks; or related to avg. peak width
                    Seg_ppm = (ROI(2)-ROI(1))/numSegs; % the width of each segment
                    MaxShift_ppm = 0.01;
                    slack_ppm = 0.0001; % must be smaller than the avg peak width by at least two data points
%                     realigns(5).Data = star_align1D(thismatrix,currentppm,'mean','CCOW',numSegs,MaxShift_ppm,slack_ppm);
%                     realigns(5).Title = 'Star Alignment - CCOW';
        % Plot the results:
              set(0,'DefaultFigureVisible','on') % turn figures back on
              ax = [];
              figure,hold on,
                    for i = 1:length(realigns)
                        ax(i) = subplot(length(realigns),2,[i*2-1,i*2]);
                            %figure, plot(currentppm,realigns(i).Data(10,:))
                            plot(currentppm,realigns(i).Data)
                            set(gca,'XDir','reverse')
                            titles{i} = ['Local Alignment for ',num2str(ROI(1)),' - ',num2str(ROI(2)),'ppm using ',realigns(i).Title];
                            title(titles{i})
                            xlabel('Chemical Shift (ppm)')
                            ylabel('Signal Intensity')
%                                 fig = gca;
%                                 saveas(fig,fig.Title.String)    
%                                 close(gcf);
                            highlightLine({})
                            %whichLine()                         
                    end 
                    linkaxes(ax,'xy')
        % Ask which method to use:
            choice = listdlg('PromptString','Which algorithm?', ...
                                'SelectionMode','single',...
                                'ListString',...
                                {'Return to picking without replacing region',...
                                 'Guide Align - PAFFT (correlation)',...
                                 'Guide Align - PAFFT (spearman)',...
                                 'Star Alignment - PAFFT (mean)'...
                                 });
                                 %'Star Alignment - CCOW'});
                    regionRealigned = realigns(choice).Data;
                    method = realigns(choice).Title;          
                    
                    fig = gcf;
                    set(fig, 'currentaxes', ax(choice));
                    highlightROIs(ROI,max(regionRealigned(1:numel(regionRealigned))));
                    figname = ['Local Alignment for ',num2str(ROI(1)),' - ',num2str(ROI(2)),'ppm.fig'];
                    saveas(fig,figname);
                    close(fig)
    end


    function patches = gatherROIs_multiPlot(fig,axes)
        patches = cell(length(axes),1);  

        for j = 1:length(axes) % loop through the axes in the figure
            set(fig,'CurrentAxes',axes(j)); % set the current axes to the one of interest

            % Get the current figure's data
                D=get(gca,'Children'); %get the handle of the line object
                patchobjs = findall(D,'Type','Patch');  % make a list of all patch boxes in the current axis only
                recBounds = nan(2,length(patchobjs));
                    for i = 1:length(patchobjs)       % for each patch box, get the vertices and store
                        verts = patchobjs(i).Vertices;
                        recBounds(:,i) = [verts(3);verts(1)];
                    %{
                        verts = [1 2; % ll % Only need one l and one r; 1&3
                                 3 4; % lr
                                 5 6; % ur
                                 7 8];% ul
                    %}
                    end
                patches{j} = recBounds;

        end 
    end

    function [titles] = trySomeAlignments(X,ppm,save)
            i = 1;
            titles = cell(1,1);
                    fprintf('\n\tRunning multiple alignments and creating plots (this will take about a minute for ~50 spectra).\n\tDo not touch the figures!...\n')
                % Plot the unaligned data for comparison:
                    titles{i} = ['Unaligned (', num2str(i),')'];
                        matrix = X;
                        figure, plot(ppm,matrix)
                            set(gca,'XDir','reverse')
                            title(titles{i})
                            xlabel('Chemical Shift (ppm)')
                            ylabel('Signal Intensity')
                            if strcmp(save,'save')
                                fig = gca;
                                saveas(fig,fig.Title.String)    
                                close(gcf);close(gcf)    
                                i = i+1;
                            end
                % Try the Pearson PAFFT
                    titles{i} = ['Guide Align - PAFFT (correlation) (', num2str(i),')'];
                    XRALg=guide_align1D(X,ppm,'correlation','PAFFT');
                        % HCA on raw spectra shows expected groupings of blanks and samples!!
                    % Plot it:
                        matrix = XRALg;
                        figure, plot(ppm,matrix)
                            set(gca,'XDir','reverse')
                            title(titles{i})
                            xlabel('Chemical Shift (ppm)')
                            ylabel('Signal Intensity')
                            if strcmp(save,'save')
                                fig = gca;
                                saveas(fig,fig.Title.String)    
                                close(gcf);close(gcf)
                                i = i+1;
                            end
               % Try the Spearman PAFFT     
                    titles{i} = ['Guide Align - PAFFT (spearman) (', num2str(i),')'];
                    XRALg=guide_align1D(X,ppm,'spearman','PAFFT');
                        % HCA on raw spectra shows expected groupings of blanks and samples!!
                    % Plot it:
                        matrix = XRALg;
                        figure, plot(ppm,matrix)
                            set(gca,'XDir','reverse')
                            title(titles{i})
                            xlabel('Chemical Shift (ppm)')
                            ylabel('Signal Intensity')
                            if strcmp(save,'save')
                                fig = gca;
                                saveas(fig,fig.Title.String)    
                                close(gcf);close(gcf)
                                i = i+1;
                            end
               % Try the Star Align PAFFT
                    titles{i} = ['Star Alignment - PAFFT (mean) (', num2str(i),')'];
                    XRALg=star_align1D(X,ppm,'mean','PAFFT');
                    % Plot it:
                        matrix = XRALg;
                        figure, plot(ppm,matrix)
                            set(gca,'XDir','reverse')
                            title(titles{i})
                            xlabel('Chemical Shift (ppm)')
                            ylabel('Signal Intensity')
                            if strcmp(save,'save')
                                fig = gca;
                                saveas(fig,fig.Title.String)    
                                close(gcf); 
                                i = i+1;
                            end
                % Try the Star Align CCOW
                    titles{i} = ['Star Alignment - CCOW (', num2str(i),')'];
                    XRALg=star_align1D(X,ppm,'mean','CCOW');
                    % Plot it:
                        matrix = XRALg;
                        figure, plot(ppm,matrix)
                            set(gca,'XDir','reverse')
                            title(titles{i})
                            xlabel('Chemical Shift (ppm)')
                            ylabel('Signal Intensity')
                            if strcmp(save,'save')
                                fig = gca;
                                saveas(fig,fig.Title.String)    
                                close(gcf);
                                i = i+1;
                            end
                % Try the ICOSHIFT Algorithm
              %      [xCS,ints,ind,target] = icoshift('average',xP,inter,n,options,Scal)
    end
end
