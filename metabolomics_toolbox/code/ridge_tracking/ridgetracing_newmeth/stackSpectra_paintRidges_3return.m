function [h] = stackSpectra_paintRidges_3return(matrix,currentppm,horzshift,vertshift,plotTitle,peaks,pointSize,nameflag,ridgeColor,ridgeLift,thickpair)
%% this script plot stackspectra based on the nmr data
%% argument:
%%% matrix: the nmr data matrix
%%% currentppm: ppm vector
%%% horzshift: horizontal shift for spectra plot
%%% vertshift: vertical shift for spectra plot
%%% plotTitle: the title of the figure
%%% peaks: the peak/ridge table
%%% pointSize: font size
%%% ridgeColor: ridge color not used yet
%%% ridgeLift: not used yet
%%% nameflag: whether print ridge name on the spectra
%%% thickpair: [rowind linewidth] to change linewidth of specific line
%% return: h the figure handler
%% by MJ
%% modified by YW (11/28/2018)
%% add line thickness pair with rowindex to change thickness of selected line

%%test%%%%
% matrix=mathere;
% currentppm=ppmhere;
% horzshift=horzshift;
% vertshift=vertshift;
% plotTitle=plotTitle;
% peaks=peakshere;
% pointSize=10;
% nameflag=false;
%%%%%%%%%%

if ~exist('nameflag', 'var')
  nameflag=false;
end
if ~exist('thickpair', 'var')
  thickpair=[];
end
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
% Set up a mapping from each point to its ridge index (for
% coloring)
for i = 1:length(ridgeRowInds)
    ridgeNumbers{i} = repmat(i,numel(ridgeRowInds{i}),1)'; % only using the numel
end
ridgeRowInds = size(matrix,1)-[ridgeRowInds{:}]+1; % invert the points (different from flipping)
ridgeNumbers = [ridgeNumbers{:}];
%ridgeRowInds = [ridgeRowInds{:}];
allpoints_ridgePos = [ridgePPMs{:}] - horzshift * size(matrix,1);
allpoints_ridgeInt = [ridgeInts{:}];
allpoints_ridgePos_shifted = allpoints_ridgePos + horzshift * ridgeRowInds;
allpoints_ridgeInt_shifted = allpoints_ridgeInt - vertshift * ridgeRowInds;% +  ridgeLift * max(matrix(1:end));
h=figure('PaperType','<custom>','PaperSize',[24 24],'Color',[1 1 1]);hold on;
  set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')');
    % Set up the colors to correspond with ridges (not including empty
    % rows)
        %cmap = lines(length(ridgeRowInds));
    % Set up the colors to correspond with ridges (including empty
    % rows)
        cmap = lines(length(ridgeInts));
    %             annotations = regexprep({peaks.CompoundNames},'( )\d','');
    %             compoundList = unique(annotations,'stable');
    %             cmap = lines(length(compoundList));

    % Plot the spectra over the scatters, spectrum by spectrum
    for i = 1:size(matrix,1)
        adjRow = matrix(i,:)-vertshift * i;
        adjRow(adjRow<baseline) = baseline;
        %% the line index need also to be upside down as well
        if length(thickpair)~=0 && (size(matrix,1)-thickpair(1)+1)==i
          ar = area(currentppm + horzshift * i,adjRow,'BaseValue',baseline,'LineWidth',thickpair(2));%,'EdgeColor','flat');
        else
          ar = area(currentppm + horzshift * i,adjRow,'BaseValue',baseline);%,'EdgeColor','flat');
        end
        lineColor = ar.FaceColor;
        ar.FaceColor = 'w';
        ar.EdgeColor = lineColor;
        % Find the ridgepoints on this row
            thisRowRidgePoints = find(ridgeRowInds == i);
        % Plot them:
            scatter(allpoints_ridgePos_shifted(thisRowRidgePoints),allpoints_ridgeInt_shifted(thisRowRidgePoints),pointSize,cmap(ridgeNumbers(thisRowRidgePoints),:),'filled');
        %find([ridgeRowInds{:}] == i)
    end
    if ~isempty(find(~cellfun(@isempty,{peaks.CompoundNames})))&&nameflag
        % Require both a name and data:
        for i = find(and(~cellfun(@isempty,{peaks.RidgeIntensities}),~cellfun(@isempty,{peaks.CompoundNames})))
            xpos = ridgePPMs{i}(1);
                % remember that matrix is flipud'd
            ypos = matrix(end,matchPPMs(xpos - horzshift * size(matrix,1),currentppm)) - vertshift * size(matrix,1);
            text(xpos,ypos-0.1,peaks(i).CompoundNames,'Rotation',-45,'FontSize',20);
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
    set(gca,'XDir','reverse');
    if isempty(plotTitle)
        title(''); %% give the possibility of giving no title. for illustratro purpose
    else
        title(plotTitle);
    end
    xlabel('Chemical Shift (ppm)');
    ylabel('Signal Intensity');
    ax1 = gca;                   % gca = get current axis
    ax1.YAxis.Visible = 'off';
    set(gcf, 'InvertHardCopy', 'off');
    set(gca,'fontsize',20);
    set(gca,'box','off');
end
