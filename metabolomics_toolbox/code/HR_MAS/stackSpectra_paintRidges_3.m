function [h] = stackSpectra_paintRidges_3(matrix,currentppm,horzshift,vertshift,plotTitle,peaks,pointSize,plotNames,ridgeColor,ridgeLift)

    matrix = flipud(matrix);
    % Adjust ppm to match the first spectrum:
    currentppm = currentppm - horzshift * size(matrix,1);
    vertshift = vertshift*max(matrix(1:end));
    baseline = mean(matrix(end,:))-std(matrix(end,:)) - vertshift * size(matrix,1); % get rid of sides

    ridgeInts = {peaks.RidgeIntensities};
    ridgePPMs = {peaks.Ridges};
    
    % Get the Ridge points for each row
        ridgeRowInds = {peaks.RowInds};
        ridgeNumbers = cell(1,length(ridgeRowInds));
        % Set up a mapping from each ridge point to its ridge index (for
        % coloring)
%         for i = 1:length(ridgeRowInds)
%             ridgeNumbers{i} = repmat(i,numel(ridgeRowInds{i}),1)'; % only using the numel
%         end
        % Set up a mapping from each ridge point to its group number (for
        % coloring), where groups are defined based on annotation so that 
        % ridges in the same compound are the same color.
        if ~isempty({peaks.CompoundNames})
            try
                annotations = regexprep({peaks.CompoundNames},'( )\d',''); 
                [compoundList,ic,groupNumber] = unique(annotations); 
            catch        
                annotations = peaks.CompoundNames; 
                compoundList = 1:length(annotations);
                ic = compoundList;
                groupNumber = compoundList;  
            end
        end            
        for i = 1:length(ridgeRowInds)
            ridgeNumbers{i} = repmat(i,numel(ridgeRowInds{i}),1)'; % only using the numel
        end
                           
        ridgeRowInds = size(matrix,1)-[ridgeRowInds{:}]+1; % invert the points (different from flipping)        
        %ridgeRowInds = [ridgeRowInds{:}];
        allpoints_ridgePos = [ridgePPMs{:}] - horzshift * size(matrix,1);
        allpoints_ridgeInt = [ridgeInts{:}];
        
        allpoints_ridgePos_shifted = allpoints_ridgePos + horzshift * ridgeRowInds;
        allpoints_ridgeInt_shifted = allpoints_ridgeInt - vertshift * ridgeRowInds;% +  ridgeLift * max(matrix(1:end));
        
    figure('PaperType','<custom>','PaperSize',[24 24],'Color',[1 1 1]),hold on,
    
    % Set up the colors to correspond with ridges (not including empty
    % rows)    
        %cmap = lines(length(ridgeRowInds));
    % Set up the colors to correspond with ridges (including empty
    % rows) 
            cmap = lines(length(ridgeInts));
    if strcmp(ridgeColor,'compound')
        for i = 1:length(groupNumber)
            ridgeNumbers{i} = repmat(groupNumber(i),numel(ridgeInts{i}),1)'; % only using the numel
        end      
        %ridgeNumbers = ;
        cmap = lines(length(compoundList));
            cmap([4,6],:) = cmap(flip([4,6]),:);
        
    end
    if strcmp(ridgeColor,'provided')
        for i = 1:length(groupNumber)
            ridgeNumbers{i} = repmat(groupNumber(i),numel(ridgeInts{i}),1)'; % only using the numel
        end      
        %ridgeNumbers = ;
        cmap = cell2mat({peaks(ic).ridgeColors}');      
    end
    
    ridgeNumbers = [ridgeNumbers{:}]; % cast as 1D vector
    
    % Plot the spectra over the scatters, spectrum by spectrum
    for i = 1:size(matrix,1)
        adjRow = matrix(i,:)-vertshift * i;
        adjRow(adjRow<baseline) = baseline;
        ar = area(currentppm + horzshift * i,adjRow,'BaseValue',baseline);%,'EdgeColor','flat');
        lineColor = ar.FaceColor;
        ar.FaceColor = 'w';
        ar.EdgeColor = lineColor;
        % Find the ridgepoints on this row
            thisRowRidgePoints = find(ridgeRowInds == i);
        % Plot them:
            scatter(allpoints_ridgePos_shifted(thisRowRidgePoints),allpoints_ridgeInt_shifted(thisRowRidgePoints),pointSize,cmap(ridgeNumbers(thisRowRidgePoints),:),'filled');
            
        %find([ridgeRowInds{:}] == i)
    end
    if strcmp(plotNames,'plotNames')
        if ~isempty(find(~cellfun(@isempty,{peaks.CompoundNames})))
            % Require both a name and data:
            for i = find(and(~cellfun(@isempty,{peaks.RidgeIntensities}),~cellfun(@isempty,{peaks.CompoundNames})))
                xpos = ridgePPMs{i}(1);
                    % remember that matrix is flipud'd
                ypos = matrix(end,matchPPMs(xpos - horzshift * size(matrix,1),currentppm)) - vertshift * size(matrix,1); 
                text(xpos,ypos-0.1,peaks(i).CompoundNames,'Rotation',-45,'FontSize',20)
            end
        end
    end
%     if ~iscell(ridgeColor)
%         rowInds = cellfun(@fliplr,{peaks.RowInds},'UniformOutput',0);
%     
%         allpoints_ridgePos = [ridgePPMs{:}] - horzshift * size(matrix,1);
%         allpoints_ridgeInt = [ridgeInts{:}];
%         allpoints_rowInd = [rowInds{:}];
% 
%         allpoints_ridgePos_shifted = allpoints_ridgePos + horzshift * allpoints_rowInd;
%         allpoints_ridgeInt_shifted = allpoints_ridgeInt - vertshift * allpoints_rowInd;% +  ridgeLift * max(matrix(1:end));
% 
%         h = scatter(allpoints_ridgePos_shifted,allpoints_ridgeInt_shifted,10,'MarkerEdgeColor','none','MarkerFaceColor',ridgeColor);
%     else        
%         cmap = lines(length(ridgePPMs));
%         for i = 1:length(ridgePPMs) 
%             rowInds = fliplr(peaks(i).RowInds);
%             
%             allpoints_ridgePos = [ridgePPMs{i}] - horzshift * size(matrix,1);
%             allpoints_ridgeInt = [ridgeInts{i}]+ridgeLift;
%             allpoints_rowInd = rowInds;
%             allpoints_ridgePos_shifted = allpoints_ridgePos + horzshift * allpoints_rowInd;
%             allpoints_ridgeInt_shifted = allpoints_ridgeInt - vertshift * allpoints_rowInd;% +  ridgeLift * max(matrix(1:end));
%             scatter(allpoints_ridgePos_shifted,allpoints_ridgeInt_shifted,10,'MarkerEdgeColor','none','MarkerFaceColor',cmap(i,:));
%         end
%     end
    

    set(gca,'XDir','reverse')
    if isempty(plotTitle)
        title('Stack Plot')
    else 
        title(plotTitle)
    end
    xlabel('Chemical Shift (ppm)')
    ylabel('Signal Intensity')
    ax1 = gca;                   % gca = get current axis
    ax1.YAxis.Visible = 'off';
    set(gcf, 'InvertHardCopy', 'off');
    set(gca,'fontsize',20)
    set(gca,'box','off')    

end