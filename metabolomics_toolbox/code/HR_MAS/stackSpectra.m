function stackSpectra(matrix,currentppm,horzshift,vertshift,plotTitle,varargin)
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
%     vertshift = 0.005 * max(max(matrix)); % the vertical distance between spectra. I tend to scale to the spectral maximum. Usually positive. 
%     plotTitle = '';                       % a string 
%
%     stackSpectra(matrix,currentppm,horzshift,vertshift,plotTitle)

% MJ 2018
%% Parse varargins

    if ~isempty(varargin)
        % Name-value pairs
            ind = find(strcmp(varargin,'colors'));
            if ~isempty(ind)
                colors = varargin{ind(1)+1};
                colors.rgb = flipud(colors.rgb);
            end
    end

%% Make the plot    
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
            ar(i) = area(currentppm + horzshift * i,adjRow,'BaseValue',baseline);%,'EdgeColor','flat');
            lineColor = ar(i).FaceColor;
            ar(i).FaceColor = 'w';
            ar(i).EdgeColor = lineColor;
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
 
    
% Change line colors 
    
    if exist('colors','var')
        for i = fliplr(1:length(ar))
            ar(i).EdgeColor = colors.rgb(i,:);
        end
        addReasonableLegend(table2cell(colors.categories(:,1)),colors.colorList)    
    end
    
     
%% Add white triangles to cover undesirable regions
% 
%     % Get the coordinates:
%         adjMat = matrix - (  (  (1:size(matrix,1))  ) * vertshift)';
%         
%         % Calculate the endpoints of the ppm vect (beginning and end x
%         % coords)
%         
%         %figure,plotr(currentppm,adjMat)
%             bufferPoints = round(0.002*length(currentppm));
% 
%         % Calculate the endpoints of the ppm vect (beginning and end x coords)
%             adjPPM(:,1) = currentppm(1+bufferPoints) + (1:size(matrix,1)) * horzshift;
%             adjPPM(:,2) = currentppm(end-bufferPoints) + (1:size(matrix,1)) * horzshift;
%             
%             
%             
%             
%             
%             adjPPM(1,:)    
%             adjPPM(end,:)
%     
%             plot(adjPPM(1,1), adjMat(1,1),'*','MarkerSize',10)
%             plot(adjPPM(end,1), adjMat(end,1),'*','MarkerSize',10)
%             
%         %    
%             m = vertshift/horzshift;
%             x = currentppm(end)+size(matrix,1)*horzshift;
%             b = 0;
%             
%         upperLeft = [innerLeft, m*x + b  ];
%         
%         upperLeft = [innerLeft, m*x + b  ];
%         plot(upperLeft,'*','MarkerSize',15)
%         
%         
%         lowerLeft = innerLeft
%         
%         upperRight
%         lowerRight
% 
%             
%     % Triangle
%         pgon = polyshape([0 0 1 3], [0 3 3 0]);
%         plot(pgon)
% 
%     % Base Rectangle
    
end