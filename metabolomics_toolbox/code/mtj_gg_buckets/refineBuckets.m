function [buckets,refinedBounds] = refineBuckets(matrix,ppm,buckets,varargin)   
%% refineBuckets() 
%
%   Author: MTJ
%   Version: 0.3
%   Tested on Matlab Version R2020a
%   Date: JUL2020 (edited MAR2022 to add simpler run options)
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
%         The following are optional as of 2022 (i.e. refineBuckets() is
%         possible, no arguments if a figure with line data is provided):
%         
%           - matrix		'X' matrix, spectral matrix
%           - ppm			ppm value vector corresponding to matrix
%           - buckets       output structure from opt_bucketing pipeline: (optimize_optBucket -> filterBuckets_Peaks_opt)
%           
%           *** If the above are not provided, no optional params can be
%               used. 
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
%             'tempFigure'        using this flag sets figure saving to
%                                 'false'
%             'optset_ind'        if plotOptBucket_optResult() step AND  
%                                 expandBucketBounds() was skipped, an 
%                                 optimization (parameter set) from 
%                                 optimize_optBucket() still needs to be set using a 
%                                 Name-value pair: use the flag, then provide the index of the
%                                 desired parameter set as an integer in
%                                 string format. e.g:
%                                 (refineBuckets(matrix,ppm,buckets,...'optset_ind','4'...)
%
% 
% Outputs:
%
%         - buckets         struct() or 2 x n matrix of bounds (adds refineBuckets to the optOB_out struct)
%             - refineBuckets contains:
%                 - originalBuckets: 	original bins (ppm bounds) provided (derived from optOB_out.results.binsWithPeaks)
%                 - refinedBuckets: 	updated bucket list (ppm bounds) after refining
%                 - removedBuckets: 	bucket (ppm bounds) which were removed during refinement
%                 - addedBuckets:       buckets (ppm bounds) which were added during refinement
%                 - figure:             name of the output figure (saved automatically in cd())
%             - bins_refined_ figure:   the figure is saved. You can use a function
%                                       like gatherROIsFromFigures() to get features from the figure
%             ** if buckets is empty, then refineBuckets() will still allow
%             bucket drawing/refining, but buckets will simply contain the
%             bucket boundaries. Modified to handle this case 31MAR2022
%             (MTJ). This case will be equivalent to refinedBounds.
%
%         - refinedBounds               new bucket bounds (no bucketStruct)
% 
% MTJ 2020-2022 
% contact   @ judgemt@uga.edu 
%           or
%           @ mjudge@imperial.ac.uk (2022 onwards)

%% Parse varargin

    % Allow input vars to be skipped
    
        if ~(exist('matrix','var') && exist('ppm','var'))
            
            % Check for valid plot
            
                if isempty(gca)
                    error('Line data must be present in figure, or provided in ''matrix'' and ''ppm''.')
                end
            
            % Get the data from the current plot
            
                [matrix,ppm] = dataFromFig();
            
            % Now, check to make sure we got data
            
                if isempty(matrix) || isempty(ppm)
                    error('Line data must be present in figure, or provided in ''matrix'' and ''ppm''.')
                end
                
        end
 
        if ~exist('buckets','var')
            buckets = [];
            % We'll get these later
        end

        
    % Set defaults for optional params to false:
    
        expandedBuckets = false;    % default is optOB_out.results.binsWithPeaks
        previousFigure = false;     % default is optOB_out.results will be used as bucket source
        saveFig = true;
        ignoreBucketStruct = false;
        %figID = [];                 % default is gcf will be used
    % 
       
    % Set optional params to true if flag is provided
    
        if ~isempty(varargin)
            if ismember('expandedBuckets',varargin)                          
                expandedBuckets = true;
            end
            
            if ismember('optset_ind',varargin) 
                ind = str2double(varargin{find(ismember('optset_ind',varargin))+1}); % must be an integer as string
            end

            if ismember('previousFigure',varargin)                          
                previousFigure = true;
            end
            
            if ismember('tempFigure',varargin)                          
                saveFig = false;
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
    
        p = reportParams('exclude',{'optOB_out','matrix','ppm','varargin'});
        
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
                        
            if isstruct(buckets)
                if expandedBuckets
                    currentBuckets = buckets.optimized.expandedBuckets;
                else
                    currentBuckets = buckets.results(ind).binsWithPeaks;
                end
                
            else
                if isempty(buckets)
                    ignoreBucketStruct = true;
                    currentBuckets = [];
                else
                    if isa(buckets,'double')
                        % Buckets is provided as a list of regions
                            currentBuckets = buckets;
                            buckets = struct(); % set as empty struct to avoid issues later (some data will be stored in there as a struct)
                            saveFig = false;
                    else
                        error('Input ''buckets'' must be a struct (in opt_bucket pipeline) or a n x 2 array of doubles where buckets(n,1)<buckets(n,2) (if simply using bucket bounds), or []')
                    end
                end
            end
            
            % Make the figure
            
                plotBuckets(matrix,ppm,currentBuckets,'')
                
        end

        %selectLines('noPause'); % enable highlighting
        
        bins = currentBuckets;
            
    %% Run Interactive Loop
        % This is the core functionality of the tool. Written so that the
        % opt_bucket pipeline object "buckets" is not required, and it only
        % operates on currentBuckets and patch objects in the active
        % figure.
        
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
                            figname = ['Buckets_refined_',num2str(now),'.fig'];
%                             buckets.refinedBuckets.figure = ['Buckets_refined_',num2str(now),'.fig'];

                    otherwise
                        
                        % Handle non-valid keypresses
                        
                            msgbox(['Sorry,   ','''',key,'''', '   isn''t a valid key command. Press ''h'' to see allowed commands.'])
                            
                end
                
                memcycle = memcycle+1;
                mem(memcycle) = {currentBuckets};
                
        end


        
    %% Clean up results, save the figure, handle different output cases
    
        % Saving figure

                if saveFig
                    fig = gca;
                    title('Refined Buckets - Complete. Please wait; saving and figure in current directory...')
                    fprintf('\n\tSaving figure, please wait...\n')
                    saveas(fig,figname)
                   % close(fig)
                    msgbox({['Refinement Completed Successfully! Figure saved as ''',figname,''''];...
                            [];...
                            [' in'];...
                            [];...
                            [' ''',cd(),'''']})
                    fprintf(['\n\tRefinement Completed Successfully! Figure saved as \n\t\t''',figname,'''\n','\tin\n\t\t','''',cd(),'''.\n'])
                    fprintf('\n\tResults saved in optOB_out.results.refinedBuckets.\n\n')
                end
                
        % Don't do any of this if there are no currentBuckets, and not if
        % there buckets was empty
        
        if isempty(currentBuckets) % allow escape with no buckets
            refinedBounds = [];
            warning('No buckets were reported')
            return
        else
            if ignoreBucketStruct
                % Just put out the bucket bounds
                    buckets = currentBuckets;
                    refinedBounds = currentBuckets;
            else
                % Build the output struct refinedBuckets within buckets
                % struct:
                
                % Record figure name 
                    buckets.refinedBuckets.figure = figname;
                    
                % Report params
                    buckets.refinedBuckets.params = p; % pass the params now

                % Remove any single-point buckets
                    buckets.refinedBuckets.singlePoints = ~(cellfun(@numel,fillRegions(matchPPMs(currentBuckets,ppm)))>1);
                        currentBuckets = currentBuckets(~buckets.refinedBuckets.singlePoints,:);
                        
                    % If all of them were removed, then error out
                        if isempty(currentBuckets)
                        error(['Only single-point buckets detected. Try:  ',...
                                '(1) Check to ensure ppm axis is appropriate for bucket bounds (this can happen when all bucket bounds ',...
                                    'map most closely to one extreme of the axis); ',...
                                '(2) Try re-running refineBuckets() using the current figure as input with the ''previousFigure'' flag'])
                    end

                % Sort the buckets and return them:
                    [~,inds] = sort(mean(currentBuckets,2));
                    currentBuckets = currentBuckets(inds,:);

                    buckets.refinedBuckets.originalBuckets = bins;
                    buckets.refinedBuckets.refinedBuckets = currentBuckets;
                    buckets.refinedBuckets.removedBuckets = bins(~ismember(bins,currentBuckets,'rows'),:);
                    buckets.refinedBuckets.addedBuckets = currentBuckets(~ismember(currentBuckets,bins,'rows'),:);
                    refinedBounds = currentBuckets;
            
%                 % Save (already done)
% 
%                 if saveFig
%                     fig = gca;
%                     title('Refined Buckets - Complete. Please wait; saving and figure in current directory...')
%                     saveas(fig,buckets.refinedBuckets.figure)
%                    % close(fig)
%                     msgbox({['Refinement Completed Successfully! Figure saved as ''',buckets.refinedBuckets.figure,''''];...
%                             [];...
%                             [' in'];...
%                             [];...
%                             [' ''',cd(),'''']})
%                     fprintf(['\n\tRefinement Completed Successfully! Figure saved as \n\t\t''',buckets.refinedBuckets.figure,'''\n','\tin\n\t\t','''',cd(),'''.\n'])
%                     fprintf('\n\tResults saved in optOB_out.results.refinedBuckets.\n\n')
%                 end
            end
        end
        
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
    %set(gcf,'WindowState','fullscreen')
    set(gcf,'WindowState','maximized')
end

function [matrix,ppm] = dataFromFig()
    % Get data from current figure
        D = get(gca,'Children'); %get the handle of the line object
        ppm = get(D,'XData'); %get the x data
            ppm = ppm{1}; % just use first ppm axis
        matrix = get(D,'YData'); %get the y data
            matrix = cell2mat(matrix);
end