function [buckets] = plotOptBucket_optResult(matrix,ppm,buckets,xlims,ylims,varargin)

%% plotOptBucket_optResult

    % Author: MTJ
    % Version: 0.1
    % Tested on Matlab Version R2020a
    % Date: JUL2020
    %
    % Description:
    %
    %       Generates a grid (postage-stamp-style) figure, where each axes 
    %       (plot) contains the result for a different opt_bucket sb-sl 
    %       parameter combination. A pop-up box appears to allow the choice 
    %       of selecting a plot from the figure with the desired 
    %       combination to be used in downstream analyses.
    %
    % Input:
    %       
    %       matrix:     spectral matrix
    %       ppm:        ppm vector
    %       buckets:    result structure from filterBuckets_Peaks_opt 
    %       xlims:      xlims for zooming (zooming manually doesn't work)
    %       ylims:      ylims for zooming (zooming manually doesn't work)
    %       varargin:
    %           'expandedBuckets'   - if using expanded bucket bounds from 
    %                                 expandBucketBounds(), add this flag
    %           'saveFig'           - if saving the resulting figure, add
    %                                 this flag
    %
    % Output:
    %
    %       A postage-stamp-style figure that is interactive and allows
    %       a optimize_optBucket() param combination to be selected
    %
    % Log:
    %       
    %       SEP2020:    Reports that manual zooming doesn't work. Earlier
    %       versions of MATLAB may have issues with the interactive bit,
    %       but the function will still work for visualization. Make sure
    %       to manually set the following if that is the case:
    %           buckets.optimization.optInd (integer index of the plot)
    %           buckets.optimized = buckets.optimization
    %       
    % Example run:
    %
    %       buckets = plotOptBucket_optResult(matrix,ppm,buckets,[0.95 1.205],[-0.0021 10.7028]);
    %
        
%% Parse varargin
    
    expandedBuckets = false;
    saveFig = false;
    
    if ~isempty(varargin)
        if ismember('expandedBuckets',varargin)                          
            expandedBuckets = true;
        end
        if ismember('saveFig',varargin)
            saveFig = true;
        end
    end    

%% Report Params

    p = reportParams('exclude',{'matrix','ppm','buckets'});

%% Calculated Vars    

    bs = buckets.optimization.optParams_OB.bucketSizes;
    sn = buckets.optimization.optParams_OB.slacknessLevels;
    optResults = buckets.optimization.results;
    
    pinds = fillRegions(matchPPMs(xlims,ppm));
    pinds = pinds{:};
    matrixss = matrix(:,pinds);
    ppmss = ppm(pinds);

    
%%

figure;

clear('ax')
% Visualize
    for i = 1:length(sn)  
        for j = 1:length(bs)
            n =  i + (j-1) * length(sn); 
            % Make the ax optOB_out.resultsject
                ax(j,i) = subplot(length(bs),length(sn),n);
                    hold on
                % 
                    plot(ppmss,matrixss)
                    
                % Plot the peak points for the max of the spectra
                    % Only plot the max peak in each bin (exploring this)          
                        %pks = optOB_out.results(n).matchingPeaks(optOB_out.results(n).maxMatchedPeakInds);
                        %scatter(pks,optOB_out.results(n).maxMatchedPeakInts,'o','k')
                        
                    scatter(optResults(n).matchingPeaks,max(optResults(n).matchingPeakInts ,[],1),'o','k')
                        set(gca,'XDir','reverse')
                        %xlabel('Chemical Shift (ppm)')
                        %ylabel('Signal Intensity')
                        set(gca,'xlim',xlims)
                        set(gca,'ylim',ylims)
                        set(gca, 'YTickLabel',[])
                        %title({['size-bucket = ',num2str(optOB_out.results(n).sb)],['slackness = ',num2str(optOB_out.results(n).sln)]})
                        
                % Highlight the bins
                if expandedBuckets
                    highlightROIs( optResults(n).expandedBuckets', max(matrixss(:)) ,'color','none','edgeColor','k')
                else
                    highlightROIs( optResults(n).binsWithPeaks', max(matrixss(:)) ,'color','r')
                end
                hold off
        end
    end
    linkaxes(ax(:),'xy');
        
    % Add a title
        sttl = suptitle('Results for opt_bucket Parameter Optimization','interpreter','none');
    
    % Add one x label across the bottom
        suplabel('Chemical Shift (ppm)','x');

    % Add the slackness params across the top
        for i = 1:length(sn)
            %xlabel(ax(1,i),['slackness = ',num2str(iv(i))],'visible','on','FontWeight','bold');
            title(ax(1,i),['slackness = ',num2str(sn(i))],'visible','on','FontWeight','bold');
        end
    % Add the bucket sizes across up left side
        for j = 1:length(bs)
            ylabel(ax(j,1),['bucket size = ',num2str(bs(j))],'visible','on','FontWeight','bold');
        end
    
    % Do the optimal result selection (if desired)
        answer = questdlg('Select Optimal Buckets?',... % question
                    'Selection Menu',...    % title
                    'Click this when you''re ready to click a set of bins to select them',... % button 1
                    'I''m not selecting bins, I just want to look this time around',... % button 2
                    'I''m not selecting bins, I just want to look this time around');   % default
           
            % Handle response
            
                switch answer
                    case 'Click this when you''re ready to click a set of bins to select them'
                        % Choose a result 
                            % Make selection by clicking on the plot
                                waitforbuttonpress; % pause till click, then set gca
                                selectedPlot = gca; % get axis handle

                            % Update the figure to indicate that the plot's been selected
                                set(selectedPlot,'Box',1)
                                set(selectedPlot,'LineWidth',4)   

                            % Use axis positions to ID the subplot of interest

                                % Get all axes positions

                                    positions = zeros(numel(ax),4); 
                                    posMap = reshape(1:numel(ax),length(bs),length(sn))';
                                    for i = 1:numel(ax)
                                         positions(posMap(i),:) = ax(i).Position;
                                    end

                                % Find the plot that matches the position of selectedPlot

                                    [~,selectedAxInd] = ismember(selectedPlot.Position,positions,'rows');                    

                            % Set that as the optOB_out result

                                buckets.optimization.optInd = selectedAxInd;
                                buckets.optimized = optResults(selectedAxInd);
                                % Other info?
                                
                            % Modify the title
%                                 delete(sttl)
%                                 suptitle(['opt_bucket Parameters selected: Slackness: ',num2str(optOB_out.results(optOB_out.optInd).sln),', Bucket Size: ',num2str(optOB_out.results(optOB_out.optInd).sb)],'interpreter','none');
                            % Save the figure
                            
                                fig = gca;
                                
                            % Save  
                                if saveFig
                                % Make figure name
                                    buckets.optimization.plotOptBuckets.selectionFigure =  ['opt_bucket optimum selected- sln ',num2str(buckets.optimized.sln),', sb ',num2str(buckets.optimized.sb),', timestamp ',num2str(now)];
                                    
                                    msgbox({;...%  Figure saved as ''',bucketStruct.refineBuckets.figure,''''
                                            [];...
                                            ['Please wait, saving figure as '];...
                                            [];...
                                            [buckets.optimization.plotOptBuckets.selectionFigure];...
                                            [];...
                                            ['in'];...  
                                            [];...
                                            [cd()];...
                                            [];...
                                            ['All data will be recorded in buckets.optimized (returned by this function).']})
                                    saveas(fig,[buckets.optimization.plotOptBuckets.selectionFigure,'.fig'])                                
                                    fprintf(['\n\tFigure saved as \n\t\t''',buckets.optimization.plotOptBuckets.selectionFigure,'''\n','\tin\n\t\t','''',cd(),'''.\n'])
                                else
                                    buckets.optimization.plotOptBuckets.selectionFigure =  'NA';
                                    msgbox('Figure was not saved. Results saved in buckets.optimized.');
                                end
                                
                                fprintf('\n\tResults saved in buckets.optimized.\n\n')
                                
                                %close(fig)
                                               
                            
                    case 'I''m not selecting bins, I just want to look this time around'
                        msgbox(['Make sure you rerun plotOptBucket_optResult() to regenerate this figure and pick a set of bins,',...
                            ' or set   optOB_out.optInd   manually if you intend to use them in downstream functions.'])

                end
        
            
end