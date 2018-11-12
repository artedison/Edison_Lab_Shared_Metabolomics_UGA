%% STEP 4: Plotting Compound Relative Intensity over Time

%save('compounds_22JUN2018.mat')
%load('compounds_22JUN2018.mat')
%save('compounds_26SEP2018.mat')
% load('compounds_26SEP2018.mat')
load('compounds.mat')
load('sampleData.mat')%the dataset we provided
%% Plot the means of the scaled ridges by condition
% mkdir('AllPlots_updated_26SEP2018_scaled')
% cd('AllPlots_updated_26SEP2018_scaled')
mkdir('AllPlots_updated_scaled')
cd('AllPlots_updated_scaled')

% Plot the trajectories of the different compounds as a function of time
%     mkdir('compounds_plots_relativeQuant')
%     cd('compounds_plots_relativeQuant')
plotList=unique({compounds.Name});
for c = 1:length(plotList)
    figure('PaperType','<custom>','PaperSize',[8 6],'Color',[1 1 1]),hold on,
    % Get the relevant rows
    currentRidges = compounds(find(strcmp({compounds.Name},plotList{c})));
    [conditions,conditionInds] = unique({compounds.Condition});
    for s=1:length(currentRidges)
        h(s) = plot(currentRidges(s).AverageTrajectory_times,currentRidges(s).AverageTrajectory_intensities_scaled,...
            ...%'DisplayName',[currentRidges(s).Condition,' replicate'],...
            'Color',currentRidges(s).PlotColor,...
            'LineWidth',1,...
            'Marker','o',...
            'MarkerSize',3);
            set(h(s),'DisplayName',[currentRidges(s).Condition,' replicate'])
    end
    %legend([h(conditionInds(1)),h(conditionInds(2))], conditions{1}, conditions{2})
    %xlabel('Time (h)')
    %ylabel({'Scaled Ridge Intensity','(time-wise mean)'})
    %title(['Peak Traces in an ','\itN. crassa\rm\bf culture '])
    %title([plotList{c},' (lag)'])
    title(plotList{c})
    set(gcf, 'InvertHardCopy', 'off');
    set(gca,'fontsize',40)
    set(gca,'XLim',[0,13])
        allpoints = {currentRidges.AverageTrajectory_intensities_scaled};
        maxValue = max([allpoints{:}]);
        minValue = min([allpoints{:}]);

    set(gca,'YLim',[minValue,maxValue])
    title(plotList{c})
    % Get rid of annoying whitespace on the outsides
    fig = gca;
    InSet = get(fig, 'TightInset');
    set(fig, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
    % Programmatically save as .fig and .pdf by title
    saveas(fig,fig.Title.String)
    print(['-f',num2str(get(gcf,'Number'))],fig.Title.String,'-dpdf')

end
cd ..

%% Plot the means of the shifted ridges by condition
% mkdir('AllPlots_updated_26SEP2018_shifted')
% cd('AllPlots_updated_26SEP2018_shifted')
mkdir('AllPlots_updated_shifted')
cd('AllPlots_updated_shifted')
% Plot the trajectories of the different compounds as a function of time
%     mkdir('compounds_plots_relativeQuant')
%     cd('compounds_plots_relativeQuant')
plotList = unique({compounds.Name});

for c = 1:length(plotList)
    figure('PaperType','<custom>','PaperSize',[8 6],'Color',[1 1 1]),hold on,
    % Get the relevant rows
    currentRidges = compounds(find(strcmp({compounds.Name},plotList{c})));
    [conditions,conditionInds] = unique({compounds.Condition});
    for s = 1:length(currentRidges)
        h(s) = plot(currentRidges(s).AverageTrajectory_times,currentRidges(s).AverageTrajectory_intensities_shifted,...
            ...%'DisplayName',[currentRidges(s).Condition,' replicate'],...
            'Color',currentRidges(s).PlotColor,...
            'LineWidth',1,...
            'Marker','o',...
            'MarkerSize',3);
            set(h(s),'DisplayName',[currentRidges(s).Condition,' replicate'])
    end
    %legend([h(conditionInds(1)),h(conditionInds(2))], conditions{1}, conditions{2})
    xlabel('Time (h)')
    ylabel('Scaled Ridge Intensity (time-wise mean)')
            %title(['Peak Traces in an ','\itN. crassa\rm\bf culture '])
            %title([plotList{c},' (lag)'])
    title(plotList{c})
    set(gcf, 'InvertHardCopy', 'off');
    set(gca,'fontsize',40)
    set(gca,'XLim',[0,13])
        allpoints = {currentRidges.AverageTrajectory_intensities_shifted};
        maxValue = max([allpoints{:}]);
        minValue = min([allpoints{:}]);
    set(gca,'YLim',[minValue,maxValue])
    title(plotList{c})
    % Get rid of annoying whitespace on the outsides
    fig = gca;
    InSet = get(fig, 'TightInset');
    set(fig, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
    % Programmatically save as .fig and .pdf by title
    %             saveas(fig,fig.Title.String)
    %             print(['-f',num2str(get(gcf,'Number'))],fig.Title.String,'-dpdf')
end
cd ..

%% Plot the means of the raw ridges by condition
% mkdir('AllPlots_updated_26SEP2018_raw')
% cd('AllPlots_updated_26SEP2018_raw')
mkdir('AllPlots_updated_raw')
cd('AllPlots_updated_raw')
% Make the colors
blue = [0,.3,1];
yellow = [1    0.7490    0.0510];%
red = [0.8,0,0];
for c = 1:length(compounds)
    if strcmp(compounds(c).Condition,'aerobic')
        compounds(c).PlotColor = red;
    else
        compounds(c).PlotColor = blue;
    end
end

% Plot the trajectories of the different compounds as a function of time
%     mkdir('compounds_plots_relativeQuant')
%     cd('compounds_plots_relativeQuant')

plotList = unique({compounds.Name});
for c = 1:length(plotList)
    figure('PaperType','<custom>','PaperSize',[8 6],'Color',[1 1 1]),hold on,
    % Get the relevant rows
    currentRidges = compounds(find(strcmp({compounds.Name},plotList{c})));
    [conditions,conditionInds] = unique({compounds.Condition});
    for s = 1:length(currentRidges)
        h(s) = plot(currentRidges(s).AverageTrajectory_times,currentRidges(s).AverageTrajectory_intensities_raw,...
            ...%'DisplayName',[currentRidges(s).Condition,' replicate'],...
            'Color',currentRidges(s).PlotColor,...
            'LineWidth',1,...
            'Marker','o',...
            'MarkerSize',3);
            set(h(s),'DisplayName',[currentRidges(s).Condition,' replicate'])
    end
    %legend([h(conditionInds(1)),h(conditionInds(2))], conditions{1}, conditions{2})
    xlabel('Time (h)')
    ylabel('Scaled Ridge Intensity (time-wise mean)')
    %title(['Peak Traces in an ','\itN. crassa\rm\bf culture '])
    %title([plotList{c},' (lag)'])
    title(plotList{c})
    set(gcf, 'InvertHardCopy', 'off');
    set(gca,'fontsize',40)
    set(gca,'XLim',[0,13])
        allpoints = {currentRidges.AverageTrajectory_intensities_raw};
        maxValue = max([allpoints{:}]);
        minValue = min([allpoints{:}]);
    set(gca,'YLim',[minValue,maxValue])
    title(plotList{c})
    % Get rid of annoying whitespace on the outsides
    fig = gca;
    InSet = get(fig, 'TightInset');
    set(fig, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
    % Programmatically save as .fig and .pdf by title
    %             saveas(fig,fig.Title.String)
    %             print(['-f',num2str(get(gcf,'Number'))],fig.Title.String,'-dpdf')
end
cd ..

%%
clearvars -except Samples2 Samples compounds sampleData samples

%% Plot the raw ridges for each replicate of a compound (color by condition):
% mkdir('compounds_plots_allRidges')
% cd('compounds_plots_allRidges')
plotList = unique({compounds.Name});
for c = 1:length(plotList)
    figure('PaperType','<custom>','PaperSize',[8 6],'Color',[1 1 1]),hold on,
    % Get the relevant rows
    currentRidges = compounds(find(strcmp({compounds.Name},plotList{c})));
    [conditions,conditionInds] = unique({compounds.Condition});
    for s = 1:length(currentRidges)
        for r = 1:length(currentRidges(s).Intensities)
            h = plot(currentRidges(s).Times{r,1},currentRidges(s).Intensities{r,1},...
                ...%'DisplayName',[currentRidges(s).Condition,' replicate'],...
                'Color',currentRidges(s).PlotColor,...
                'LineWidth',1,...
                'Marker','o',...
                'MarkerSize',3);
                set(h,'DisplayName',[currentRidges(s).Condition,' replicate'])
                %clear('h')
        end
    end

    %legend([h(conditionInds(1)),h(conditionInds(2))], conditions{1}, conditions{2})
    xlabel('Time (h)')
    ylabel('Raw Ridge Intensity')
    %title(['Peak Traces in an ','\itN. crassa\rm\bf culture '])
    %title([plotList{c},' (lag)'])
    set(gcf, 'InvertHardCopy', 'off');
    set(gca,'fontsize',40)
    set(gca,'XLim',[0,13])
        allpoints = {currentRidges.Intensities};
        allpoints = [allpoints{:}];
        maxValue = max([allpoints{:}]);
        minValue = min([allpoints{:}]);
    set(gca,'YLim',[minValue,maxValue])
    title([plotList{c},' - all ridges by condition'])
    % Get rid of annoying whitespace on the outsides
    fig = gca;
    InSet = get(fig, 'TightInset');
    set(fig, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
    % Programmatically save as .fig and .pdf by title
    %             saveas(fig,fig.Title.String)
    %             print(['-f',num2str(get(gcf,'Number'))],fig.Title.String,'-dpdf')
    %
end

%% Plot the raw ridges for each replicate of a compound (color by replicate):
%  mkdir('compounds_plots_allRidges_byReplicate')
%  cd('compounds_plots_allRidges_byReplicate')
mkdir('Figure3_allridges_raw')
cd('Figure3_allridges_raw')
plotList = unique({compounds.Name});
cmap = lines(6);
for c = 1:length(plotList)
    figure('PaperType','<custom>','PaperSize',[8 6],'Color',[1 1 1]),hold on,
    % Get the relevant rows
    currentRidges = compounds(find(strcmp({compounds.Name},plotList{c})));
    %currentRidges = compounds(find(     and(strcmp({compounds.Name},plotList{c}) , [compounds.SampleNumber]==1)    )); % all ridges for one rep
    [conditions,conditionInds] = unique({compounds.Condition});
    for s = 1:length(currentRidges)
        for r = 1:length(currentRidges(s).Intensities)
            h = plot(currentRidges(s).Times{r,1},currentRidges(s).Intensities{r},...
                ...%'DisplayName',[currentRidges(s).Condition,' replicate'],...
                'Color',cmap(s,:),...
                ...%'Color',cmap(6,:),...
                'LineWidth',1,...
                'Marker','o',...
                'MarkerSize',3);
                set(h,'DisplayName',[currentRidges(s).Condition,' replicate'])
                %clear('h')
            % Plot the means
            % h = plot(currentRidges(s).AverageTrajectory_times,currentRidges(s).AverageTrajectory_intensities_raw,...
            %                     ...%'DisplayName',[currentRidges(s).Condition,' replicate'],...
            % ...%'Color',red,...
            % 'Color',cmap(4,:),...
            % 'LineWidth',3,...
            % 'MarkerSize',3);
            % set(h,'DisplayName',[currentRidges(s).Condition,' replicate'])                      %clear('h')
        end
    end
    %legend([h(conditionInds(1)),h(conditionInds(2))], conditions{1}, conditions{2})
    % xlabel('Time (h)')
    % ylabel('Raw Ridge Intensity')
    %   title(['Peak Traces in an ','\itN. crassa\rm\bf culture '])
    %   title([plotList{c},' (lag)'])
    set(gcf, 'InvertHardCopy', 'off');
    set(gca,'fontsize',40)
    set(gca,'XLim',[0,13])
        allpoints = {currentRidges.Intensities};
        allpoints = [allpoints{:}];
        maxValue = max([allpoints{:}]);
        minValue = min([allpoints{:}]);
    set(gca,'YLim',[minValue,maxValue])
    %title([plotList{c},' - all ridges by replicate'])
    title(plotList{c})
    % Get rid of annoying whitespace on the outsides
    fig = gca;
    InSet = get(fig, 'TightInset');
    set(fig, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
    % Programmatically save as .fig and .pdf by title
    saveas(fig,fig.Title.String)
    print(['-f',num2str(get(gcf,'Number'))],fig.Title.String,'-dpdf')

end

%% Plot the shifted ridges for each replicate of a compound (color by replicate):
%  mkdir('compounds_plots_allRidges_byReplicate')
%  cd('compounds_plots_allRidges_byReplicate')
plotList = unique({compounds.Name});
cmap = lines(6);
for c = 1:length(plotList)
    figure('PaperType','<custom>','PaperSize',[8 6],'Color',[1 1 1]),hold on,
    % Get the relevant rows
    currentRidges = compounds(find(strcmp({compounds.Name},plotList{c})));
    [conditions,conditionInds] = unique({compounds.Condition});

    for s = 1:length(currentRidges)
        for r = 1:length(currentRidges(s).Intensities)
            h = plot(currentRidges(s).Times{r,1},currentRidges(s).trajectory_ScaledRidgeIntensities_shifted{r},...
                ...%'DisplayName',[currentRidges(s).Condition,' replicate'],...
                'Color',cmap(s,:),...
                'LineWidth',1,...
                'Marker','o',...
                'MarkerSize',3);
                set(h,'DisplayName',[currentRidges(s).Condition,' replicate'])
                %clear('h')
        end
    end

    %legend([h(conditionInds(1)),h(conditionInds(2))], conditions{1}, conditions{2})
    xlabel('Time (h)')
    ylabel('Shifted Ridge Intensity')
    %title(['Peak Traces in an ','\itN. crassa\rm\bf culture '])
    %title([plotList{c},' (lag)'])
    set(gcf, 'InvertHardCopy', 'off');
    set(gca,'fontsize',40)
    set(gca,'XLim',[0,13])
        allpoints = {currentRidges.trajectory_ScaledRidgeIntensities_shifted};
        allpoints = [allpoints{:}];
        maxValue = max([allpoints{:}]);
        minValue = min([allpoints{:}]);

    set(gca,'YLim',[minValue,maxValue])
    title([plotList{c},' - all ridges by replicate'])
    % Get rid of annoying whitespace on the outsides
    fig = gca;
    InSet = get(fig, 'TightInset');
    set(fig, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
    % Programmatically save as .fig and .pdf by title
    % saveas(fig,fig.Title.String)
    % print(['-f',num2str(get(gcf,'Number'))],fig.Title.String,'-dpdf')
end

%% Plot the scaled ridges for each replicate of a compound (color by replicate):
%  mkdir('compounds_plots_allRidges_byReplicate')
%  cd('compounds_plots_allRidges_byReplicate')
mkdir('Figure3_allridges_scaled')
cd('Figure3_allridges_scaled')
plotList = unique({compounds.Name});
cmap = lines(6);
for c = 1:length(plotList)
    figure('PaperType','<custom>','PaperSize',[8 6],'Color',[1 1 1]),hold on,
    % Get the relevant rows
    currentRidges = compounds(find(  strcmp({compounds.Name},plotList{c})   )); % all ridges all reps
    %currentRidges = compounds(find(     and(strcmp({compounds.Name},plotList{c}) , [compounds.SampleNumber]==1)    )); % all ridges for one rep
    [conditions,conditionInds] = unique({compounds.Condition});
    % G1P issue with sample 1, ridge 4
    for s = 1:length(currentRidges)
        for r = 1:length(currentRidges(s).trajectory_ScaledRidgeIntensities_scaled)
            h = plot(currentRidges(s).Times{r,1},currentRidges(s).trajectory_ScaledRidgeIntensities_scaled{r},...
                ...%'DisplayName',[currentRidges(s).Condition,' replicate'],...
                'Color',cmap(s,:),...
                ...%'Color',cmap(6,:),...
                'LineWidth',1,...
                'Marker','o',...
                'MarkerSize',3);
                set(h,'DisplayName',[currentRidges(s).Condition,' replicate'])
            % Plot the means
            % h = plot(currentRidges(s).AverageTrajectory_times,currentRidges(s).AverageTrajectory_intensities_scaled,...
            %     ...%'DisplayName',[currentRidges(s).Condition,' replicate'],...
            %     ...%'Color',red,...
            %     'Color',cmap(4,:),...
            %     'LineWidth',3,...
            %     'MarkerSize',3);
            %     set(h,'DisplayName',[currentRidges(s).Condition,' replicate'])                      %clear('h')
        end
    end
    %legend([h(conditionInds(1)),h(conditionInds(2))], conditions{1}, conditions{2})
    % xlabel('Time (h)')
    % ylabel('Scaled Ridge Intensity (time-wise mean)')
    %     title(['Peak Traces in an ','\itN. crassa\rm\bf culture '])
    %     title([plotList{c},' (lag)'])
    set(gcf, 'InvertHardCopy', 'off');
    set(gca,'fontsize',40)
    set(gca,'XLim',[0,13])
        allpoints = {currentRidges.trajectory_ScaledRidgeIntensities_scaled};
        allpoints = [allpoints{:}];
        maxValue = max([allpoints{:}]);
        minValue = min([allpoints{:}]);

    set(gca,'YLim',[minValue,maxValue])
    % title([plotList{c},' - all ridges by replicate'])
    title(plotList{c})

    % Get rid of annoying whitespace on the outsides
    fig = gca;
    InSet = get(fig, 'TightInset');
    set(fig, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
    % Programmatically save as .fig and .pdf by title
    saveas(fig,fig.Title.String)
    print(['-f',num2str(get(gcf,'Number'))],fig.Title.String,'-dpdf')

end

%% Change properties of existing axes without calling them forward:
thisAx = findall(0,'type','axes');
for i = 1:length(thisAx)
   set(thisAx(i),'fontsize',20)
end

reorderFigures('lowToHigh')

%% Make a stackspectra plot from sampleData and compounds
%% Supplement: Plot DSS Stability over time
samples = [1:3,7:9];
% Use Raw Intensity
figure, hold on
for i = samples
    vect = sampleData(i).normFactors;
    %vect = vect - mean(vect);
    %vect = vect/max(vect);
    plot(sampleData(i).timesR_1h1d, vect,'DisplayName',[num2str(i),' (',sampleData(i).condition,')'])
end
legend
xlabel('Time (h)')
ylabel('Raw Intensity')
title(['DSS intensity (Full Resolution)'])
% Use scaled, normalized
figure, hold on
for i = samples
    vect = sampleData(i).normFactors;
    vect = vect - mean(vect);
    vect = vect/max(vect);
    plot(sampleData(i).timesR_1h1d, vect,'DisplayName',[num2str(i),' (',sampleData(i).condition,')'])
end
legend
xlabel('Time (h)')
ylabel('Scaled, Normalized Intensity')
title(['DSS intensity (Full Resolution)'])

%% Supplement: Plot several baselines over time

%ROI = [0.5,0.55];
ROI = [0.6];
%ROI = [6.4];
%ROI = [8.8];

% Use Raw intensity
figure, hold on
for i = samples
    newROIinds = matchPPMs(ROI,sampleData(i).ppmR_1h1d);
    %vect = sum(sampleData(i).XN_1h1d(:,newROIinds(1):newROIinds(2))'); % AUC
    vect = sampleData(i).XR_1h1d(:,newROIinds(1)); % intensity at a point
    %vect = vect - mean(vect);
    %vect = vect/max(vect);

    plot(sampleData(i).timesR_1h1d, vect, 'DisplayName',[num2str(i),' (',sampleData(i).condition,')'])
end
legend
ylabel('Raw Intensity')
xlabel('Time (h)')
%title(['Baseline (Full Resolution) at ',num2str(ROI(1)),'-',num2str(ROI(2)),'ppm'])
title(['Baseline (Full Resolution) at ',num2str(ROI(1)),'ppm'])

% Use scaled, normalized
figure, hold on
for i = samples
    newROIinds = matchPPMs(ROI,sampleData(i).ppmR_1h1d);
        %vect = sum(sampleData(i).XN_1h1d(:,newROIinds(1):newROIinds(2))'); % AUC
    vect = sampleData(i).XR_1h1d(:,newROIinds(1)); % intensity at a point
    vect = vect - mean(vect);
    vect = vect/max(vect);

    plot(sampleData(i).timesR_1h1d, vect, 'DisplayName',[num2str(i),' (',sampleData(i).condition,')'])
end
legend
ylabel('Scaled, Normalized Intensity')
xlabel('Time (h)')
%title(['Baseline (Full Resolution) at ',num2str(ROI(1)),'-',num2str(ROI(2)),'ppm'])
title(['Baseline (Full Resolution) at ',num2str(ROI(1)),'ppm'])

%% Calculate correlations between two compounds
cmpd1 = 'citrate';
cmpd2 = 'choline';
cmpd1_inds = find(strcmp({compounds.Name},cmpd1));
cmpd2_inds = find(strcmp({compounds.Name},cmpd2));
% Get the data, times, and repNumber for each compound
cmpd1_trajectories = {compounds(cmpd1_inds).AverageTrajectory_intensities_scaled};
cmpd1_times = compounds(cmpd1_inds).AverageTrajectory_times;
cmpd1_repNumber = [compounds(cmpd1_inds).SampleNumber];

cmpd2_trajectories = {compounds(cmpd2_inds).AverageTrajectory_intensities_scaled};
cmpd2_times = compounds(cmpd2_inds).AverageTrajectory_times;
cmpd2_repNumber = [compounds(cmpd2_inds).SampleNumber];

% Which samples to compare with those of cmpd1?
[~,~,order2] = intersect(cmpd1_repNumber , cmpd2_repNumber);

% For taking derivative, need to adjust the time vector (use values in
% between timepoints)
cmpd1_times = cmpd1_times - (cmpd1_times(end)-cmpd1_times(1))/length(cmpd1_times);
cmpd1_times = cmpd1_times(2:end);
% cmpd1_times = cmpd1_times - (cmpd1_times(end)-cmpd1_times(1))/length(cmpd1_times);
% cmpd1_times = cmpd1_times(2:end);
cmpd2_times = cmpd2_times - (cmpd2_times(end)-cmpd2_times(1))/length(cmpd2_times);
cmpd2_times = cmpd2_times(2:end);
% cmpd2_times = cmpd2_times - (cmpd2_times(end)-cmpd2_times(1))/length(cmpd2_times);
% cmpd2_times = cmpd2_times(2:end);

% Calculate the correlations between two compounds of the same sample, for each sample
for i = 1:length(order2)
    % Average trajectories
    traj1 = cmpd1_trajectories{i}';
    traj2 = cmpd2_trajectories{order2(i)}';
    % Mean-center and normalize
    traj1 = (traj1-mean(traj1))/max(traj1);
    traj2 = (traj2-mean(traj2))/max(traj2);
    % Take the derivative
    traj1 = diff(traj1);
    traj2 = diff(traj2);

    [correlations(i,2),correlations(i,3)] = corr(traj1,traj2,'Type','Pearson');
    correlations(i,1) = cmpd1_repNumber(i);

    figure,hold on
        plot(cmpd1_times,traj1,'r','DisplayName',cmpd1)
        plot(cmpd2_times,traj2,'b','DisplayName',cmpd2)
        title([cmpd1,' rep ',num2str(cmpd1_repNumber(i)),' - and - ',cmpd2,' rep ',num2str(cmpd1_repNumber(order2(i)))])
        legend
        hold off
    figure,hold on
        plot(traj1,traj2,'*k')
        title([cmpd1,' rep ',num2str(cmpd1_repNumber(i)),' - vs - ',cmpd2,' rep ',num2str(cmpd1_repNumber(order2(i)))])
end
