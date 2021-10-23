function [noise,p] = stackSpectra(matrix,ppm,horzshift,vertshift,plotTitle,varargin)
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

    noise = quantile(sort(abs(matrix(:))),0.25); % useful later
    whiteShapes = true;
    continuousColors = false;
    timeVect = (1:size(matrix,1))';
    
    if ~isempty(varargin)
        % Name-value pairs
        
        % If we have colors specified
            ind = find(strcmp(varargin,'colors'));
            if ~isempty(ind)
                colors = varargin{ind(1)+1};
                %colors.rgb = flipud(colors.rgb);
                if ~isfield(colors,'categories') % figure out if it's a colorCategories or customColormap struct
                    continuousColors = true;
                end
%             % If we then have indices specified for those colors
%                 ind = find(strcmp(varargin,'colorInds'));
%                 if ~isempty(ind)
%                     colorInds = varargin{ind(1)+1};
%                     colors.rgb = flipud(colors.rgb(colorInds,:)); % subset them, then flip to match matrix
%                 end  
        % If we are plotting based on a time vector (should come after
        % plotInds
        
            ind = find(strcmp(varargin,'timeVect'));
            if ~isempty(ind)
                %tvflag = 1;
                timeVect = varargin{ind(1)+1};
                if size(timeVect,2) > size(timeVect,1) % make sure it's a column
                    timeVect = timeVect';
                end
            end
                
            % If we then have plot indices specified
                ind = find(strcmp(varargin,'plotSubset'));
                if ~isempty(ind)
                    plotInds = varargin{ind(1)+1};               
                    matrix = matrix(plotInds,:);                 % don't flip this yet
                    timeVect = timeVect(plotInds);
                    if exist('colors','var')
                        colors.rgb = flipud(colors.rgb(plotInds,:)); % subset them, then flip to match matrix
                    end
                else
                    if exist('colors','var')
                        colors.rgb = flipud(colors.rgb);             % flip to match matrix
                    end
                end
            end

            if any(strcmp(varargin,'autoVert'))
                vertshift = vertshift * noise;
                %fprintf(['\n\n\tVertshift -> noise multiple mode. Estimated noise level : ',num2str(noise),'\n'])
                timeVect = (1:size(matrix,1))'; % override timeVect (only way)
            end     

            
            if any(strcmp(varargin,'noWhiteShapes'))
                whiteShapes = false;
            end     
            
            
    end

    p = reportParams('exclude',{'matrix','ppm','timeVect','colors'});
    
%% Do some calculations ahead of time

    % Vertshift and horzshift need to be related 
    
    matrix = flipud(matrix);
%%    
    hshiftvect = horzshift * timeVect; 
    
    shiftedppms = flipud(repmat(ppm,length(timeVect),1) - repmat(hshiftvect,1,length(ppm)));
    
    % Adjust ppm to match the first spectrum in shiftedmat
    
        %ppm = ppm + horzshift * max(timeVect);
        
    % All vertical adjustments made here
    
        vshiftvect = flipud(vertshift * timeVect); 
        
        shiftedmat = matrix + vshiftvect; 
        
            baseline = mean(shiftedmat(end,:))-std(shiftedmat(end,:)) - vertshift * length(timeVect); % get rid of sides
            
        shiftedmat(shiftedmat<baseline) = baseline;

    
%% Make the plot       
    % Spectra are plotted from the back to the front, and top to bottom
    % (i.e. the last row of the matrix is plotted highest and first, and
    % the rest are plotted on top of it. area() is used to mask spectra
    % farther back, although it is heavy on the graphics. 
    
        figure('PaperType','<custom>','PaperSize',[24 24],'Color',[1 1 1]),hold on,
        
% Make adjustments during plotting
%         for i = 1:size(matrix,1)
%             adjRow = matrix(i,:)-vertshift * i;
%             adjRow(adjRow<baseline) = baseline;
%             ar(i) = area(currentppm + horzshift * i,adjRow,'BaseValue',baseline);%,'EdgeColor','flat');
%             lineColor = ar(i).FaceColor;
%             ar(i).FaceColor = 'w';
%             ar(i).EdgeColor = lineColor;
%         end
% Adjustments made as matrices
    
        for i = 1:size(matrix,1)
            ar(i) = area(shiftedppms(i,:),shiftedmat(i,:),'BaseValue',baseline);%,'EdgeColor','flat');
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
 
    
% Update line colors (if provided)
    
    if exist('colors','var')
        if continuousColors
            for i = 1:length(ar) % *** must update them in the correct order so they lay on top of each other correctly
                ar(i).EdgeColor = colors.rgb(i,:); % colors.rgb was already flipped to match matrix
            end
            
            colorbar(gca,'Colormap',colors.cmap,...
                 'TickLabels',round(colors.setPts,1),...
                 'Ticks', (colors.setPts - min(colors.setPts))/range(colors.setPts),...
                 'FontSize',20);
             
        else
            for i = 1:length(ar) % *** must update them in the correct order so they lay on top of each other correctly
                ar(i).EdgeColor = colors.rgb(i,:); % colors.rgb was already flipped to match matrix
            end
            legendInds = unique(colors.inds_cat(plotInds));
            addReasonableLegend(table2cell(colors.categories(legendInds,1)),colors.colorList(legendInds,:))    % these were not flipped
        end
    end
    
%     set(gca,'xlim',[-3.21862204651590,-3.04023643628845])
%     set(gca,'ylim',[-126132930513.596,73262839879.1538])
     
%% Add white shapes to cover undesirable regions on edges and @ baseline

% NEED TO UPDATE so vert is calculated from a single point, not rows 1 or end
    if whiteShapes
        % Get the coordinates:

                bufferSpace = (shiftedppms(1,1) - shiftedppms(end,1) ) *1.1;

                % R shape (parallelogram on top of a rectangle):

                    % X coords
                        trap1(1,1) = shiftedppms(1,1) + bufferSpace;    % upper left
                        trap1(2,1) = shiftedppms(1,1) - bufferSpace;	% upper right      
                        trap1(3,1) = shiftedppms(end,1) - bufferSpace;  % middle right       
                        trap1(4,1) = shiftedppms(end,1) - bufferSpace;  % lower right
                        trap1(5,1) = shiftedppms(end,1) + bufferSpace;  % lower left
                        trap1(6,1) = shiftedppms(end,1) + bufferSpace;  % middle left

                    % Y coords
                        trap1(1,2) = shiftedmat(1,1)+vertshift*2;     % upper left
                        trap1(2,2) = shiftedmat(1,1)+vertshift*2;     % upper right      
                        trap1(3,2) = shiftedmat(end,1)-vertshift*2;   % middle right       
                        trap1(4,2) = baseline-vertshift;              % lower right
                        trap1(5,2) = baseline-vertshift;              % lower left
                        trap1(6,2) = shiftedmat(end,1)-vertshift*2;   % middle left

                    pgonR = fill(trap1(:,1),trap1(:,2),'w','EdgeColor','w','HandleVisibility','off');

                % L shape (parallelogram on top of a rectangle): 

                    % X coords
                        trap2(1,1) = shiftedppms(1,end) + bufferSpace;      % upper left
                        trap2(2,1) = shiftedppms(1,end) - bufferSpace;      % upper right      
                        trap2(3,1) = shiftedppms(end,end) - bufferSpace;    % middle right       
                        trap2(4,1) = shiftedppms(end,end) - bufferSpace;    % lower right       
                        trap2(5,1) = shiftedppms(end,end) + bufferSpace;    % lower left
                        trap2(6,1) = shiftedppms(end,end) + bufferSpace;    % middle left 

                    % Y coords
                        trap2(1,2) = shiftedmat(1,end)+vertshift*2;         % upper left      
                        trap2(2,2) = shiftedmat(1,end)+vertshift*2;         % upper right
                        trap2(3,2) = shiftedmat(end,end)-vertshift*2;       % middle right
                        trap2(4,2) = baseline-vertshift;                    % lower right
                        trap2(5,2) = baseline-vertshift;                    % lower left
                        trap2(6,2) = shiftedmat(end,end)-vertshift*2;       % middle left

                    pgonL = fill(trap2(:,1),trap2(:,2),'w','EdgeColor','w','HandleVisibility','off');                

            % Add the rectangle at the bottom
                % Get the coords
                    % x 

                        xl = get(gca,'xlim'); % matlab plots farther than I would, so we'll just cover the whole axis
                        rec(1,1) = xl(2);
                        rec(2,1) = xl(1);
                        rec(3,1) = rec(2,1);
                        rec(4,1) = rec(1,1);
                    % y     
                        rec(1,2) = baseline + vertshift;
                        rec(2,2) = baseline + vertshift;
                        rec(3,2) = baseline - vertshift;
                        rec(4,2) = baseline - vertshift;

                % Plot it:
                        fill(rec(:,1),rec(:,2),'w','EdgeColor','w','HandleVisibility','off');
    end           
   %% Other odds and ends            
        % Zoom in to the relevant plot region
            set(gca,'xlim',[min(shiftedppms(:)),max( shiftedppms(:) )])

        % Do away with the colored background
            set(gcf, 'InvertHardCopy', 'off');
            hold off
    
end