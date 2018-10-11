function [output] = ridgeTracing_PeakPick1D(matrix,ppm,times,ROI,smoothingFilter,peakPickThreshold)

%% Setup, filter, and peak pick 
    ROIinds = matchPPMs(ROI(1),ppm):matchPPMs(ROI(2),ppm);
    windowPPMs = ppm(ROIinds(1:end));
    window = matrix(:,ROIinds);
    % Smooth the window to make peakpicking easier:
        %window = filter2(fspecial('average',[1,5]),window);
        eval(sprintf('window = %s;',smoothingFilter))
        %window = imgaussfilt(window,[1,15]);
set(0,'DefaultFigureVisible','off')

row = cell(1,size(window,1));
col = cell(1,size(window,1));

for i = 1:size(window,1)
    [~,col{i}] = Peakpick1D_noFigures(window(i,:),1:size(window,2),1,peakPickThreshold,'Simple');
    row{i} = repmat(i,1,length(col{i}));
end

peakInds = sub2ind(size(window),cell2mat(row),cell2mat(col));
set(0,'DefaultFigureVisible','on')
    
    %Surface plot of window
        figure, hold on
            surf(ppm(ROIinds(1:end)),1:size(window,1),window,'FaceColor','Interp');         % Plot the smoothed spectra as surface
            scatter3(windowPPMs(cell2mat(col)),cell2mat(row),window(peakInds),'.r')         % Plot the picked peaks   
                %s.CData = S./repmat(1.5*std(S,[],1),size(S,1),1);
                shading interp
                set(gca,'xdir','reverse') 
                xlabel('ppm')
                zlabel('Intensity')
                ylabel('Time (h)')  
                set(gcf, 'InvertHardCopy', 'off');
                set(gca,'fontsize',20)
                set(gca,'box','off')
                title('Gaussian-smoothed, peak-picked spectra')               

           
           % Save a structure to pass to ridgeTracing_clusterPeaks()      
            output.matrix = matrix;
            output.ppm = ppm;
            output.times = times;
            output.ROI = ROI;
            output.ROIinds = ROIinds;
            output.smoothingFilter = smoothingFilter;
            output.peakPickThreshold = peakPickThreshold;
            output.window = window;
            output.windowPPMs = windowPPMs;
            output.cols = cell2mat(col);
            output.rows = cell2mat(row);
            output.peakInds = peakInds;
            
            
end