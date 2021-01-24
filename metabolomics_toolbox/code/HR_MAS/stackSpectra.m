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

    noise = quantile(sort(abs(matrix(:))),0.25); % useful later
            %figure,plot(1:length(sort(abs(matrix(:)))),sort(abs(matrix(:))))
    if ~isempty(varargin)
        % Name-value pairs
        
        % If we have colors specified
            ind = find(strcmp(varargin,'colors'));
            if ~isempty(ind)
                colors = varargin{ind(1)+1};
                %colors.rgb = flipud(colors.rgb);
                
            % If we then have indices specified for those colors
                ind = find(strcmp(varargin,'colorInds'));
                if ~isempty(ind)
                    colorInds = varargin{ind(1)+1};
                    colors.rgb = flipud(colors.rgb(colorInds,:)); % subset them, then flip to match matrix
                end  
                
            end
                        
            if ~isempty(strcmp(varargin,'autoVert'))
                vertshift = vertshift * noise;
                fprintf(['\n\n\tVertshift -> noise multiple mode. Estimated noise level : ',num2str(noise),'\n'])
            end     
            
    end

%% Make the plot    

    matrix = flipud(matrix);
    baseline = mean(matrix(end,:))-std(matrix(end,:)) - vertshift * size(matrix,1); % get rid of sides

    hshiftvect = (horzshift * (1:size(matrix,1)))';
    shiftedppms = flipud(repmat(currentppm,size(matrix,1),1) - repmat(hshiftvect,1,length(currentppm)));

    vshiftvect = vertshift * (1:size(matrix,1))';
    shiftedmat = matrix - vshiftvect;
    
    % Adjust ppm to match the first spectrum:
        currentppm = currentppm - horzshift * size(matrix,1);
    
%%        
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
 
    
% Update line colors (if provided)
    
    if exist('colors','var')
        for i = 1:length(ar) % *** must update them in the correct order so they lay on top of each other correctly
            ar(i).EdgeColor = colors.rgb(i,:); % colors.rgb was already flipped to match matrix
        end
        %addReasonableLegend(table2cell(colors.categories(:,1)),colors.colorList)    % these were not flipped
    end
    
%     set(gca,'xlim',[-3.21862204651590,-3.04023643628845])
%     set(gca,'ylim',[-126132930513.596,73262839879.1538])
     
%% Add white triangles to cover undesirable regions
% NEED TO UPDATE so vert is calculated from a single point, not rows 1 or end

    % Get the coordinates:
%         adjMat = matrix - (  (  (1:size(matrix,1))  ) * vertshift)';
        
        % Calculate the endpoints of the ppm vect (beginning and end x
        % coords)
        
        %figure,plotr(currentppm,adjMat)
            %bufferPoints = round(0.002*length(currentppm));
            bufferSpace = (shiftedppms(1,1) - shiftedppms(end,1) ) *1.1;
            
        % Calculate the endpoints of the ppm vect (beginning and end x coords)
%             adjPPM(:,1) = currentppm(1+bufferPoints) + (1:size(matrix,1)) * horzshift;
%             adjPPM(:,2) = currentppm(end-bufferPoints) + (1:size(matrix,1)) * horzshift;
%             adjPPM(:,1) = shiftedppms(:,1);
%             adjPPM(:,2) = shiftedppms(:,end);
            
            % R pentagon: 
            
                % X coords
                    trap1(1,1) = shiftedppms(1,1) + bufferSpace;     % upper middle
                    trap1(2,1) = shiftedppms(1,1) - bufferSpace;	% upper left      
                    trap1(3,1) = shiftedppms(end,1) - bufferSpace;   % lower left       
                    trap1(4,1) = shiftedppms(end,1) - bufferSpace;      % lower right
                    trap1(5,1) = shiftedppms(end,1) + bufferSpace;      % upper right 
                    trap1(6,1) = shiftedppms(end,1) + bufferSpace;  % dummy
                     
                % Y coords
                    trap1(1,2) = shiftedmat(1,1)+vertshift*2;       % upper middle
                    trap1(2,2) = shiftedmat(1,1)+vertshift*2;     % upper left
                    trap1(3,2) = shiftedmat(end,1)-vertshift*2;                % lower left
                    trap1(4,2) = baseline-vertshift;                % lower right
                    trap1(5,2) = baseline-vertshift;       % upper right
                    trap1(6,2) = shiftedmat(end,1)-vertshift*2;       % dummy
                    
                pgonR = fill(trap1(:,1),trap1(:,2),'w','EdgeColor','w','HandleVisibility','off');
%                 for p = 1:size(trap1,1)
%                     p = 6
%                     plot(trap1(p,1), trap1(p,2),'*','MarkerSize',10,'HandleVisibility','off');
%                 end                
                
%%
            % L pentagon: 
            
                % X coords
                    trap2(1,1) = shiftedppms(1,end) + bufferSpace;
                    trap2(2,1) = shiftedppms(1,end) - bufferSpace;	% upper middle      
                    trap2(3,1) = shiftedppms(end,end) - bufferSpace;   % lower left       
                    trap2(4,1) = shiftedppms(end,end) - bufferSpace;   % lower left       
                    trap2(5,1) = shiftedppms(end,end) + bufferSpace;      % lower right
                    trap2(6,1) = shiftedppms(end,end) + bufferSpace;      % upper right 

                % Y coords
%                     trap2(1,2) = shiftedmat(1,end)+vertshift*2;       % upper right
                    trap2(1,2) = shiftedmat(1,end)+vertshift*2;       % upper right
                    trap2(2,2) = shiftedmat(1,end)+vertshift*2;     % upper middle
                    trap2(3,2) = shiftedmat(end,end)-vertshift*2;   % lower left
                    trap2(4,2) = baseline-vertshift;                % lower right
                    trap2(5,2) = baseline-vertshift;       % upper right
                    trap2(6,2) = shiftedmat(end,end)-vertshift*2;       % upper right
                    
                pgonL = fill(trap2(:,1),trap2(:,2),'w','EdgeColor','w','HandleVisibility','off');                

%                 for p = 1:size(trap1,1)
%                     plot(trap2(p,1), trap2(p,2),'*','MarkerSize',10,'HandleVisibility','off');
%                 end
%                 delete(marks(p))

% Add the rectangle at the bottom
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
    
                    fill(rec(:,1),rec(:,2),'w','EdgeColor','w','HandleVisibility','off');
                    
                    set(gca,'xlim',[min(shiftedppms(:)),max( shiftedppms(:) )])
                    
    
end