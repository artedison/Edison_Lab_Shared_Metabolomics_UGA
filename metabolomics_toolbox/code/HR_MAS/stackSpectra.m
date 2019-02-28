function stackSpectra(matrix,currentppm,horzshift,vertshift,plotTitle)
%% stackSpectra
%{
    This function plots spectra sequentially so that trends can be seen for
    peaks across samples instead of simply overlaying them. The plot is
    functionally like plot(), but uses area() instead so that the area
    under each spectrum blocks what lies behind it. This is critical for
    interpretability. Exporting this in a vectorized way, however, takes 
    forever if possible at all.

%}
%% Usage
%     matrix = X;                           % your spectral matrix (each row is a spectrum)
%     currentppm = ppm;                     % ppm vector
%     horzshift = 0.01;                     % the horizontal shift (in ppm) between sequential spectra. 
%     vertshift = 0.005 * max(max(matrix)); % the vertical distance between spectra. I tend to scale to the spectal maximum.
%     plotTitle = '';                       % a string 
%
%     stackSpectra(matrix,currentppm,horzshift,vertshift,plotTitle)

% MJ 2018
%% 
%vertshift = vertshift*max(matrix(1:end)); % MJ edit 27SEP2018. Not
%necessary to scale vertshift if data are normalized. It causes issues. 
matrix = flipud(matrix);
%baseline = mean(matrix(1:end,:))-std(matrix(1:end,:)) - vertshift * size(matrix,1); % get rid of sides
baseline = mean(matrix(end,:))-std(matrix(end,:)) - vertshift * size(matrix,1); % get rid of sides
% Adjust ppm to match the first spectrum:
currentppm = currentppm - horzshift * size(matrix,1);

    figure('PaperType','<custom>','PaperSize',[24 24],'Color',[1 1 1]),hold on,
    for i = 1:size(matrix,1)
        adjRow = matrix(i,:)-vertshift * i;
        adjRow(adjRow<baseline) = baseline;
        ar = area(currentppm + horzshift * i,adjRow,'BaseValue',baseline);%,'EdgeColor','flat');
        lineColor = ar.FaceColor;
        ar.FaceColor = 'w';
        ar.EdgeColor = lineColor;
    end
    
    
    
    set(gca,'XDir','reverse')
    if isempty(plotTitle)
        title('Stack Plot')
    else 
        title(plotTitle,'Interpreter','none')
    end
    xlabel('Chemical Shift (ppm)')
    ylabel('Signal Intensity')
    ax1 = gca;                   % gca = get current axis
    ax1.YAxis.Visible = 'off';
    set(gcf, 'InvertHardCopy', 'off');
    set(gca,'fontsize',20)
    set(gca,'box','off')
end