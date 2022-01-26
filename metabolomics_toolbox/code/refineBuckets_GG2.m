function [bucketStruct] = refineBuckets_GG2(matrix,ppm,bucketStruct,per_spec_peak_struct,ppeak,varargin)   
% refineBuckets() is a tool to assist users with refining buckets drawn by opt_bucket
% 
% 	(Reference:  	SAA Sousa, A Magalhães, MMC Ferreira, Optimized Bucketing for NMR spectra: Three
% 					case studies. Chemometrics and Intelligent Laboratory Systems, 122, 93-102 (2013).)
% 					 
% and modified by functions in the Edison Lab Metabolomics Toolbox.
% 
% You will usually run this function after optimizing opt_bucket parameters, generating the buckets, and 
% filtering to include only those buckets which include peaks (per Peakpick1d on each spectrum).
% Once running, you can zoom around the spectrum and check out your buckets, which are highlighted.
% When you find a region you want to modify, use key commands to switch modes.'
% Modified version of refineBuckets soo it displays individual peak heights

%% see at the end for usage and sequential worflow

% Inputs: 
%         - matrix		'X' matrix, spectral matrix
%         - ppm			ppm value vector corresponding to matrix
%         - optOB_out 	output structure from opt_bucketing pipeline: (optimize_optBucket -> filterBuckets_Peaks_opt)
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
% Outputs:
%         - optOB_out		(adds refineBuckets to the optOB_out struct)
%             - refineBuckets contains:
%                 - originalBuckets: 	original bins (ppm bounds) provided (derived from optOB_out.results.binsWithPeaks)
%                 - refinedBuckets: 	updated bucket list (ppm bounds) after refining
%                 - removedBuckets: 	bucket (ppm bounds) which were removed during refinement
%                 - addedBuckets:       buckets (ppm bounds) which were added during refinement
%                 - figure:             name of the output figure (saved automatically in cd())
%             - bins_refined_ figure:   the figure is saved. You can use a function
%                                       like gatherROIsFromFigures() to get features from the figure
%         - 

% Key commands include: 
% 
%         'a'- add new bucket 
%         'r'- select buckets to remove
%         'c'- select buckets to combine
%         'e'- select buckets to expand/redraw
%         'h'- help
%         'q'- quit/finish refinement
% 
% 
% Hope this helps
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
            
            
                plotBuckets(matrix,ppm,currentBuckets,'')
                
                hold on
                for i= 1:size(matrix,1)
                    plotr( per_spec_peak_struct.iteration(i).shifts,  per_spec_peak_struct.iteration(i).ints,'c.','MarkerSize',9)
                    hold on
                end
                plotr( ppeak.shifts, ppeak.ints, '.b', 'MarkerSize',12)
        end

        %selectLines('noPause'); % enable highlighting
        
        
        bins = currentBuckets;
            %%

                              
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

   %% Optimize Peak Picking (threshold; for representative spectrum)
%         
%                                 matrix = wrk_data.XRA;
%                                 ppm = wrk_data.ppmR;
%                                      optimize_Peakpick1D(matrix,ppm,'var',0.05:0.05:0.25,'Complex');   
%                             %% Do the actual Peak Picking (for representative spectrum)
%                                      peaks = struct();
%                                      [peaks.ints, peaks.shifts]= Peakpick1D(matrix ,ppm,'max',0.2,'Complex');
%                        
%                             %% Automatic Binning (Bucketing)
%  
%                             %% Generate buckets using a range of both params
% 
%                                 sb = 0.002:0.002:0.008;
%                                 sl = 0.3:0.1:0.6;
%                                 [optOB_out] = optimize_optBucket(matrix,ppm,sb,sl);
%                             %% Filter out the bins with no peaks
%                                [optOB_out] = filterBuckets_Peaks_opt(ppm,optOB_out, peaks);        
%                             %% Plot the results of optOB
%                                  [optOB_out] =  plotOptBucket_optResult(matrix,ppm,optOB_out,[3.6584    4.0], [min(matrix(:)) max( matrix(:,(ppm>3.6584  & ppm<4.0)), [],'all' ) ]);
%                          %% Peak pick every spectra independently, so that each spectrum contributes to the peak shape between the boundaries
%                                 peakthresh = 0.15 % adjust for moore or fewer peaks
%                                 mode = 'Complex'      
%                                 [itpeaks] = Peakpick1D_per_spectra(matrix,ppm,peakthresh,mode)
%                          %%  Expand buckets to each of the bins boundaries
%                                  [optOB_out] = expandBucketBounds (optOB_out, matrix, ppm, 'plotResult');
%                         %%  Manual Refinement of the boundaries
%                             %% [optOB_out] = refineBuckets(matrix,ppm,optOB_out,5);    
%                                  [optOB_out] = refineBuckets_GG2(matrix,ppm,optOB_out,itpeaks,peaks, 'expandedBuckets');    
%                              %% saving variables of the workspace 
%                                save('post_buckets_refined_xx_xx_2022.mat')
% %% This calculates the number of peaks in each bin in the original Peakpick1D_per_spectra data
%         % Calculates the max peak within each bin and its chemical shift 
%         % gap fills both ppm and intensities of peaks that were no present
%         % in the bin or not picked (not detected/below baseline)
%  matrix = wrk_data.XRA;
%  ppm = wrk_data.ppmR; 
%  buckets = optOB_out.refinedBuckets.refinedBuckets;
%  ppeak_struct = itpeaks;
% [itpeaks] = peaks_per_bin (matrix, ppm,buckets, ppeak_struct);
% %% Align metric
%  buckets = itpeaks.max_iteration.sorted_bucket_list;
%  ppeak_struct = itpeaks;
%  metadata = Td_pd_ugt;
% itpeaks = align_metric ( matrix, ppm, ppeak_struct, buckets, metadata );
% clearvars matrix ppm 
% %% getting blank matrix peak picked
% blank_matrix = wrk_data.X_blks;
% ppm = wrk_data.ppm; 
% [peaks.blank_ints, peaks.blank_shifts]= Peakpick1D(blank_matrix ,ppm,'max',0.9,'Complex');
% %% scoring blank peaks 
% spectra_struct = wrk_data.X_pd_ugt;
% pick_output = peaks;
% per_spectra_pick_output = itpeaks;
% sorted_bucket_list = itpeaks.max_iteration.sorted_bucket_list;
%  [itpeaks] = blank_feature_score (spectra_struct, blank_matrix, ppm, pick_output,  ...
%                                                         per_spectra_pick_output, sorted_bucket_list)
% clearvars spectra_struct ppick_output per_spectra_pick_output sorted_bucket_list
% clearvars blank_matrix ppm

