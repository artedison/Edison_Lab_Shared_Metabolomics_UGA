function [newBounds] = expandToValleys(bounds,X,ppm)


% MTJ 2020
        %buckets = bucketStruct.optimized.binsWithPeaks;
        
%             if all(bounds(:,1)<bounds(:,2))
%                 
%             end
                
        % Lowest point in the representative spectrum between two bins
        
            % Get the lowest point in the max spectrum of each non-bin region
            
                    
                    betweenBuckets = [bounds(1:end-1,2),bounds(2:end,1)];
                    
%                 if all(bounds(:,1)<bounds(:,2))
%                     
%                 else
%                     betweenBuckets = [bounds(2:end,2), bounds(1:end-1,1)];
%                 end

%                         figure,
%                         hold on
%                         plot(ppm,X)
%                             set(gca,'XDir','reverse')
%                             set(gca, 'YTickLabel',[])
%                             highlightROIs( betweenBuckets', max(X(:)) ,'color','k')
%                             highlightROIs( bins', max(X(:)) ,'color','none','edgeColor','k')

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
                
                newBounds = bounds;
                bps = breakPoints;
                newBounds(2:end,1) = ppm(bps(:,2))';
                newBounds(1:end-1  ,2) = ppm(bps(:,2))';


end