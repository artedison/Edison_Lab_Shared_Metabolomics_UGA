function [bucketStruct] = refineBuckets_GG_dec(matrix,ppm,bucketStruct,per_spec_peak_struct,,varargin)   
%% refineBuckets() 
%
%   Author: MTJ
%   Version: 0.2
%   Tested on Matlab Version R2020a
%   Date: JUL2020
%
%   Description:  
%
%           refineBuckets()  is a tool to assist users with refining buckets drawn by opt_bucket
% 
%           (Reference:  	SAA Sousa, A Magalhães, MMC Ferreira, Optimized Bucketing for NMR spectra: Three
%           				case studies. Chemometrics and Intelligent Laboratory Systems, 122, 93-102 (2013).)
% 					 
%           and modified by functions in the Edison Lab Metabolomics Toolbox.
% 
%           You will usually run this function after optimizing opt_bucket parameters, generating the buckets, and 
%           filtering to include only those buckets which include peaks (per Peakpick1d on each spectrum).
%           Once running, you can zoom around the spectrum and check out your buckets, which are highlighted.
%           When you find a region you want to modify, use key commands to
%           switch modes.
%
%           Key commands include: 
% 
%               'a'- add new bucket 
%               'r'- select buckets to remove
%               'c'- select buckets to combine
%               'e'- select buckets to expand/redraw
%               'h'- help
%               'q'- quit/finish refinement
% 
% Inputs: 
%
%         - matrix		'X' matrix, spectral matrix
%         - ppm			ppm value vector corresponding to matrix
%         - optOB_out 	output structure from opt_bucketing pipeline: (optimize_optBucket -> filterBuckets_Peaks_opt)
%
%         Optional (name-value pairs):
%             'expandedBuckets'   (usual case) using the valley-expanded buckets rather than traditional buckets 
%             'previousFigure'    use buckets refined from previous figure
%                                 as a starting point. Buckets from the figure
%                                 will be used, but data come from
%                                 optOB_out input
%             'figID'             If figID is provided, extractROIs will determine if it's
%                                 a figure filename or figure handle and act accordingly.
%                                 If not provided, it will default to empty, in which case gcf
%                                 will be used.
%
% 
% Outputs:
%
%         - bucketStruct		(adds refineBuckets to the optOB_out struct)
%             - refineBuckets contains:
%                 - originalBuckets: 	original bins (ppm bounds) provided (derived from optOB_out.results.binsWithPeaks)
%                 - refinedBuckets: 	updated bucket list (ppm bounds) after refining
%                 - removedBuckets: 	bucket (ppm bounds) which were removed during refinement
%                 - addedBuckets:       buckets (ppm bounds) which were added during refinement
%                 - figure:             name of the output figure (saved automatically in cd())
%             - bins_refined_ figure:   the figure is saved. You can use a function
%                                       like gatherROIsFromFigures() to get features from the figure
% 
% 
% MTJ 2020

%% Parse varargin

    % Set defaults for optional params to false:
    
        expandedBuckets = false;    % default is optOB_out.results.binsWithPeaks
        previousFigure = false;     % default is optOB_out.results will be used as bucket source
        %figID = [];                 % default is gcf will be used
    
    % Set optional params to true if flag is provided
    
        if ~isempty(varargin)
            if ismember('expandedBuckets',varargin)                          
                expandedBuckets = true;
            end

            if ismember('previousFigure',varargin)                          
                previousFigure = true;
            end
            
            if ismember('figID',varargin)                          
                % If figID is provided, extractROIs will determine if it's
                % a figure filename or figure handle and act accordingly.
                % If not provided, it will default to empty, for which gcf
                % will be used.
                [~,ind] = ismember('figID',varargin);
                figID = varargin{ind + 1};
            end
        end    
    
    % Record the params
    
        bucketStruct.refinedBuckets.params = reportParams('exclude',{'optOB_out','matrix','ppm','varargin'});
        
    %% Get the current Buckets

        % After this section, currentBuckets and bins will be set.
        
        if previousFigure
            
            % Extract from a figure and open it
            if exist('figID', 'var')
                [currentBuckets,~,fig] = extractROIs(figID);
                figure(fig);
            else
                 [currentBuckets,~,fig] = extractROIs(gcf);
                figure(fig);
            end
        else
            
            % Decide where to get the bins within optOB_out.results
            
            if expandedBuckets
                currentBuckets = bucketStruct.optimized.expandedBuckets;
            else
                currentBuckets = bucketStruct.results(ind).binsWithPeaks;
            end
            
            % Make the figure
            figure,
            
                        plotBuckets(matrix,ppm,currentBuckets,'')
                
                hold on
                for i= 1:size(matrix,1)
                    plotr( ceper_spec_peak_struct.ppmlist(i),  per_spec_peak_struct.iteration(i).ints,'c.','MarkerSize',9)
                    hold on
                end
               % plotr( ppeak.shifts, ppeak.ints, '.b', 'MarkerSize',12)
      
        end

        %selectLines('noPause'); % enable highlighting
        
        
        bins = currentBuckets;
            
    %% Run Interactive Loop
        
        key = '.';                      % dummy value for initialization
        memcycle = 1;
        mem = {currentBuckets};
        
        while ~strcmp(key,'q')          % run until 'quit' command

            % Find the existing buckets and apply updated bucket info to
            % figure
            
                [~,patches,~] = extractROIs();
                updateFigure(matrix,currentBuckets,patches,expandedBuckets)
                
            % Collect user input by keypress

                pause()
                key = get(gcf,'currentcharacter');
                key = lower(key);

            % Take an action depending on the keypress
            
                switch key
                    case 'a'

                        % Add a bucket
                            title('Click-drag a rectangle to draw new bucket boundaries (only x coords matter)')
                            currentBuckets = drawBucket(currentBuckets);

                    case 'd'

                        % Select one or more buckets by drawing box
                            title('Click and drag a rectangle to select buckets to delete') % provide instructions
                            remBuckets  = selectBuckets(currentBuckets);
                            currentBuckets(remBuckets,:) = [];

                    case 'c'

                        % Combine selected buckets
                            title('Click and drag a rectangle to select buckets to combine') % provide instructions
                            currentBuckets = combineBuckets(currentBuckets);

                    case 'e'

                        % Replace selected buckets with one new bucket
                            title('Redraw the boundaries by click-dragging a rectangle. Bounds included in the rectangle will be redrawn. ') % provide instructions
%                             currentBuckets = redrawBucket(currentBuckets);
                            [currentBuckets] = redrawBuckets_exp(currentBuckets);
                            
                    case 's'
                    
                        % Split a bucket into two
                            title('Click where you want to split the bucket') % provide instructions
                            currentBuckets = splitBucket(currentBuckets);

%                     case 'u'
%                     
%                         % Undo the previous action
%                             memcycle = memcycle-1;
%                             currentBuckets = mem{memcycle};                            
                            
%                     case 'r'
%                     
%                         % Redo the previous action
%                             currentBuckets = mem(memcycle);
%                             memcycle = memcycle-1;
                                                        
                    case 'h'
                        
                        % Present a box with the help menu
                        
                            msgbox(helpText())

                    case 'q'
                        
                        % Quit/Exit the interactive part of the program
                            bucketStruct.refinedBuckets.figure = ['Buckets_refined_',num2str(now),'.fig'];

                            title('Refined Buckets - Complete. Please wait; saving and figure in current directory...')

                    otherwise
                        
                        % Handle non-valid keypresses
                        
                            msgbox(['Sorry,   ','''',key,'''', '   isn''t a valid key command. Press ''h'' to see allowed commands.'])
                            
                end
                
                memcycle = memcycle+1;
                mem(memcycle) = {currentBuckets};
                
        end


        
    %% Clean up results, save the figure, 
           
        % Remove any single-point buckets
            bucketStruct.refinedBuckets.singlePoints = ~(cellfun(@numel,fillRegions(matchPPMs(currentBuckets,ppm)))>1);
            currentBuckets = currentBuckets(~bucketStruct.refinedBuckets.singlePoints,:);
            
        % Sort the buckets and return them:
            bucketStruct.refinedBuckets.originalBuckets = bins;
            bucketStruct.refinedBuckets.refinedBuckets = currentBuckets;
            bucketStruct.refinedBuckets.removedBuckets = bins(~ismember(bins,currentBuckets,'rows'),:);
            bucketStruct.refinedBuckets.addedBuckets = currentBuckets(~ismember(currentBuckets,bins,'rows'),:);
            
            
        % Save 
         
            fig = gca;
            saveas(fig,bucketStruct.refinedBuckets.figure)
           % close(fig)
            msgbox({['Refinement Completed Successfully! Figure saved as ''',bucketStruct.refinedBuckets.figure,''''];...
                    [];...
                    [' in'];...
                    [];...
                    [' ''',cd(),'''']})
            fprintf(['\n\tRefinement Completed Successfully! Figure saved as \n\t\t''',bucketStruct.refinedBuckets.figure,'''\n','\tin\n\t\t','''',cd(),'''.\n'])
            fprintf('\n\tResults saved in optOB_out.results.refinedBuckets.\n\n')

end

function txt = helpText()

   txt = {['HELP BOX:'];...
        [];...
        ['refineBuckets() is a tool to assist users with refining buckets drawn by opt_bucket (',...
        'Reference:  	SAA Sousa, A Magalhães, MMC Ferreira, Optimized Bucketing for NMR spectra: Three',...
        'case studies. Chemometrics and Intelligent Laboratory Systems, 122, 93-102 (2013).) ',...
        'and modified by opt_bucket pipeline functions in the Edison Lab Metabolomics Toolbox.',...
        ' You will usually run this function after optimizing opt_bucket parameters, generating the buckets, ',...
        ' filtering to include only those buckets which include peaks (per Peakpick1d on each spectrum), and ',...
        ' expanding bucket boundaries to valleys between buckets.',...
        ' Once running, you can zoom around the spectrum and check out your buckets, which are highlighted.',...
        ' When you find a region you want to modify, use key commands to switch modes.'];...
        [];...
        ['Key commands include: '];...
        [];...
        ['''a''- add new bucket by drawing it'];...
        ['''d''- remove by clicking or selecting bucket(s)'];...
        ['''c''- combine buckets by selecting them'];...
        ['''e''- expand/contract/redraw (combined remove/add)'];...
        ['''s''- split a bucket at the desired location'];...        
        ['''u''- undo previous action'];...        
        ['''r''- redo previous action'];...        
        ['''h''- help box'];...
        ['''q''- quit/finish refinement'];...
        [];...
%         ['At any time when zooming/panning (i.e. not during a key command action), '],...
%         ['use the data cursor in the MATLAB figure window to highlight/un-highlight individual '],...
%         ['spectra to better examine them (NOTE: the spectra will de-highlight after an action '],...
%         ['is taken.'],...
        [];...
        ['Hope this helps'];...
        [];...
        [];...
        ['MTJ 2020']};
end

function updateFigure(matrix,currentBuckets,patches,expandedBuckets)

    delete(patches)
    
    if expandedBuckets
        highlightROIs(currentBuckets',max(matrix(:)) ,'color','r','edgeColor','k')
    else
        highlightROIs(currentBuckets',max(matrix(:)) ,'color','r')
    end
    
    title('Press one of the allowed keys for a refining action. Press ''h'' for help')
    
end