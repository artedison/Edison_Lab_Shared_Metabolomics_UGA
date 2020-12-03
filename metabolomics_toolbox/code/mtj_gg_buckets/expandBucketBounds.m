function [bucketStruct] = expandBucketBounds(bucketStruct,X,ppm,varargin)

%% expandBinBounds
%
%     Takes buckets from opt_bucket pipeline (optimize_optBucket -> optOB_out 
%     (+ peaks -> filterBuckets_Peaks_opt), and expands the bucket bounds to the
%     lowest point between each bucket (i.e. bounds = valleys). 
%
%     Potential Issues:  - huge bins possible
%                        - overlapping or perfectly adjacent bins 
% 
%     Inputs:
%         buckets     buckets to be expanded to lowest valley between bins
%         X           spectral matrix
%         ppm         ppm vector
%         varargin    only supports flags:
%
%                       'plotResult' (default is don't plot)
% 
%     Outputs: 
%
%         bucketStruct  modified with field for expanded bin bounds (ppm) 
%
%     Usage: 
%
%           [buckets] = expandBucketBounds(buckets,matrix,ppm,'plotResult');
%
%   MTJ 2020
%

    %% Handle varargin

        if ~isempty(varargin)
            if ismember('plotResult',varargin)                          
                plotResult = true;
            else
                plotResult = false;
            end
        end
        
    %% Expand the bins
    
        buckets = bucketStruct.optimized.binsWithPeaks;
        
        % Lowest point in the representative spectrum between two bins
        
            % Get the lowest point in the max spectrum of each non-bin region

                betweenBuckets = [buckets(2:end,2), buckets(1:end-1,1)];

                        % figure,
                        % hold on
                        % plot(ppm,X)
                        %     set(gca,'XDir','reverse')
                        %     set(gca, 'YTickLabel',[])
                        %     highlightROIs( betweenBuckets', max(X(:)) ,'color','k')
                        %     highlightROIs( bins', max(X(:)) ,'color','none','edgeColor','k')

                breakPoints = zeros(size(betweenBuckets));

                for i = 1:size(betweenBuckets,1)
                    pinds = fillRegions(matchPPMs(betweenBuckets(i,:),ppm));
                    reg = pinds{:};
                        %figure,
                        %plot(ppm(reg),max(    X(:,reg),[],1));
                    [breakPoints(i,1),breakPoints(i,2)] = min(  max(    X(:,reg),[],1)   );
                    breakPoints(i,2) = reg(breakPoints(i,2));
                end

%                         figure,
%                             hold on
%                             plot(ppm,X)
%                             % Plot the peak points for the max of the spectra
%                                 scatter(ppm(breakPoints(:,2)),breakPoints(:,1),'o','k')
%                                 set(gca,'XDir','reverse')
%                                 xlabel('Chemical Shift (ppm)')
%                                 ylabel('Signal Intensity')
%                                 set(gca, 'YTickLabel',[])
%                                 % Plot the bins
%                                 highlightROIs( bins', max(X(:)) ,'color','r')
%                         
%                         hold off
                
        % Set new bin bounds to breakpoints
                
                newBuckets = buckets;
                bps = breakPoints;
                newBuckets(1:end-1,1) = ppm(bps(:,2))';
                newBuckets(2:end  ,2) = ppm(bps(:,2))';
                
    %% Plot expanded bin bounds (if desired)
    
        if plotResult
            
                figure,
                    hold on
                    plot(ppm,X)

            % Plot the peak points for the max of the spectra
                %scatter(ppm(breakPoints(:,2)),breakPoints(:,1),'o','k')
                    set(gca,'XDir','reverse')
                    xlabel('Chemical Shift (ppm)')
                    ylabel('Signal Intensity')
                    set(gca, 'YTickLabel',[])
                    title('Expanded Buckets - Lowest points method')

                % Draw the new bin bounds
                highlightROIs( newBuckets', max(X(:)) ,'color','none','edgeColor','k')

                hold off
        end
        
    %% Try expansion based on (smoothed) low slope threshold
        % (maybe in future iterations)
    
    %% Store result
    
        bucketStruct.optimized.expandedBuckets = newBuckets;
        
        
end