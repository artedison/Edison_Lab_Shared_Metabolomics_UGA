function stackSpectra_mirror_2(matrix1,matrix2,currentppm1,currentppm2,times1,times2,horzshift,vertshift,baselineSpacing,subTitles,plotTitle)
%% stackSpectra_mirror
%{
    This function plots spectra sequentially so that trends can be seen for
    peaks across samples instead of simply overlaying them. The plot is
    functionally like plot(), but uses area() instead so that the area
    under each spectrum blocks what lies behind it. This is critical for
    interpretability. Exporting this in a vectorized way, however, takes 
    forever if possible at all.

%}
%% Usage
%     matrix = Xcollapsed_1h1d;                         % a spectral matrix (each row is a spectrum)
%     matrix2 = Xcollapsed_1h1d;                        % another spectral matrix (each row is a spectrum)
%     currentppm = ppm_1h1d;                            % a ppm vector for matrix 1
%     currentppm2 = ppm_1h1d;                           % the ppm vector for matrix 2
%     horzshift = -0.001;                               % the horizontal shift (in ppm) between sequential spectra. 
%     vertshift = 0.01;                                 % the vertical distance between spectra.
%     baselineSpacing = 0.0001;                         % controls spacing between the two samples' spectra
%         plotTitle = {'Comparison of...'};             % Plot title
%                 %uisetcolor();                        % this is helpful for getting the desired colors, example results stored below 
%                 yellow = [1.0000    0.7490    0.0510];
%                 red = [0.800000000000000,0,0]; 
%         subTitles = {'sample1','sample2';...          % text boxes to label the upper and lower halves
%                      yellow,                    red;...
%                      8,                          8;...
%                      0.2,                        0.1};
%     
%     stackSpectra_mirror(matrix,matrix2,currentppm,currentppm2,horzshift,vertshift,baselineSpacing,subTitles,plotTitle)

% MJ 2018

%% Setup stuff

    baseline = 0;  % center on zero
    baseShift1 = baselineSpacing * max(max([matrix1(1:end),matrix2(1:end)]));   % shift the first matrix up 
    baseShift2 = 0-baselineSpacing * max(max([matrix1(1:end),matrix2(1:end)])); % shift the second matrix down

    vertshift2 = vertshift*times2;
    vertshift1 = vertshift*times1;
    
    horzshift2 = horzshift*times2;
    horzshift1 = horzshift*times1;
    
%% Plot from the bottom up for matrix2 (below)    
    % Plot the second matrix, inverted, underneath
        matrix2 = baseShift2-( matrix2-min( matrix2(1:end) ) ); % invert the zero-bottomed matrix
    
    figure('PaperType','<custom>','PaperSize',[24 24],'Color',[1 1 1]),hold on,
    
    for i = fliplr(1:size(matrix2,1))
        adjRow = matrix2(i,:)-vertshift2(i);% * i;      
%         ar = area(currentppm2 - horzshift2(i) * i,adjRow,'BaseValue',baseline);
        ar = area(currentppm2 - horzshift2(i),adjRow,'BaseValue',baseline);
        lineColor = ar.FaceColor;
        ar.FaceColor = 'w';
        ar.EdgeColor = lineColor;
    end    
        %set(gca,'XAxisLocation','top');
    
    
%% Plot matrix1 as usual    
        
matrix1 = flipud(matrix1);
matrix1 = matrix1-min( matrix1(1:end) ) + vertshift*size(matrix1,1) + baseShift1; % get the zero-bottomed matrix
% Adjust ppm to match the first spectrum:
%     currentppm1 = currentppm1 - horzshift1 * size(matrix1,1);
    currentppm1 = currentppm1 - horzshift1(end);

    for i = 1:size(matrix1,1)
        adjRow = matrix1(i,:)-vertshift1(i);% * i;
%         ar = area(currentppm1 + horzshift(i) * i,adjRow,'BaseValue',baseline);
        ar = area(currentppm1 + horzshift1(i),adjRow,'BaseValue',baseline);
        lineColor = ar.FaceColor;
        ar.FaceColor = 'w';
        ar.EdgeColor = lineColor;
    end
    
    
    set(gca,'XDir','reverse')
    if isempty(plotTitle)
        title('Stack Plot')
    else 
        title(plotTitle)
    end
    xlabel('Chemical Shift (ppm)')
    ylabel('Signal Intensity')
    ax1 = gca;                   
    ax1.YAxis.Visible = 'off';
    set(gcf, 'InvertHardCopy', 'off');
    set(gca,'fontsize',20)
    if ~isempty(subTitles)
        text(subTitles{3,1},subTitles{4,1}*max(matrix1(end,:)+vertshift*size(matrix1,1)),subTitles{1,1},'FontSize',14,'Color',subTitles{2,1})
        text(subTitles{3,2},subTitles{4,2}*min(matrix2(end,:)-vertshift*size(matrix1,1)),subTitles{1,2},'FontSize',14,'Color',subTitles{2,2})
    end
end