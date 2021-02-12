function [regions,peaksPerBin] = guessRegions(matrix,ppm,peaksPerRegion_target,thresh)


            
	%% Number of Peaks (very strict) ~ equal for all
    
    if ~exist('thresh','var')
        thresh = 0.8;
        response = 'initiate';
        while ~isempty(response)
            [peaks,shifts,~] = Peakpick1D_noFigures(matrix,ppm,'max',thresh,'Complex');
            f = gcf;
            title('Close figure to proceed')
            uiwait(f,60)
            response = inputdlg(['Last peakpick threshold setting: ',num2str(thresh),' . New threshold (cancel to accept last one):']);
            if ~isempty(response)
                thresh = str2num(response{:});
            end
        end
    else
        [peaks,shifts,~] = Peakpick1D_noFigures(matrix,ppm,'max',thresh,'Complex');
    end   
    
    %% Generate initial bins by number of peaks

        % Evenly space the bin dividers among the indices
            binDividers = shifts(1:peaksPerRegion_target:end);
                binDividers(end) = max(shifts);     % make sure the end gets included with the remainder (not exact peak #)
                
        % Build the initial bounds
            bins = zeros(2,length(binDividers)-1);
            bins(1,:) = binDividers(1:end-1);
            bins(2,:) = binDividers(2:end);
        
        % Shrink the bins (to allow expansion to natural valley points)
%                         avgBinSize = max(ppm)/length(binDividers);
            binSizes = bins(2,:)-bins(1,:);
            bins(1,2:end) = bins(1,2:end) + 0.25 * binSizes(2:end);
            bins(2,1:end-1) = bins(2,1:end-1) - 0.25 * binSizes(1:end-1);
        %figure,hold on,plotr(ppm,max(matrix,[],1)),highlightROIs(bins,max(matrix(:)), 'edgeColor','k'  )

    %% Expand Bounds to valleys
        
        [regions] = expandToValleys(bins',matrix,ppm)';
        
        
%         figure,hold on,plotr(ppm,max(matrix,[],1));
%                        highlightROIs(regions,max(matrix(:)), 'edgeColor','k'  );
%                        plot(shifts,max(peaks,[],1),'k*')
        
        
        [~,~,~,~,peaksPerBin] = filterBins_peaks(regions',shifts);
        
    %% Generate initial bins with equal spacing
    
% %         nbins = floor(ppm(end)/regionsize);
% %         bounds = zeros(2,nbins);
% %         bounds(1,:) = ppm(  matchPPMs(  linspace(ppm(1),ppm(end)-regionsize*0.75,nbins), ppm  )  );
% %         bounds(2,:) = ppm(  matchPPMs(  linspace(ppm(1)+regionsize*0.75,ppm(end),nbins), ppm  )  );
% %         
% %         figure,hold on,plotr(ppm,max(matrix,[],1)),highlightROIs(bounds,max(matrix(:)))
%     % Expand Bounds to valleys
%         
%         [regions] = expandToValleys(bins',matrix,ppm)';
%         
%         
% %         figure,hold on,plotr(ppm,max(matrix,[],1));
% %                        highlightROIs(regions,max(matrix(:)), 'edgeColor','k'  );
% %                        plot(shifts,max(peaks,[],1),'k*')
% %         
% %         
%             [binsWithPeaks,matchingPeaks,binFilt,peakFilt,peaksPerBin] = filterBins_peaks(regions',shifts);
%             
%             %figure,plot(peaksPerBin(binFilt))
%     
%             regions = regions(:,binFilt)';
%             
%             %figure,plot(ppm(matchPPMs(shifts,ppm)),find(shifts))

                    
    %% Plot em            
        
%         figure,hold on,plotr(ppm,max(matrix,[],1));
%                        highlightROIs(regions,max(matrix(:)), 'edgeColor','k'  );
%                        plot(shifts,max(peaks,[],1),'k*')
                       
        [~,regions] = refineBuckets(matrix,ppm,regions','expandedBuckets');
        % Regions are sorted above ^
        close(gcf)
        regions = regions';
        
        figure,hold on,plotr(ppm,max(matrix,[],1));
                       highlightROIs(regions,max(matrix(:)), 'edgeColor','k'  );
                       plot(shifts,max(peaks,[],1),'k*')
        
        
                                          
end