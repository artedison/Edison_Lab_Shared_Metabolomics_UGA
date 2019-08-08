%% this script is to trace ridges from experimental data
%%%ridges that have been traced in previous publication will be repeated here
%%%you can look into details of the previous publication in:
%%% Judge Michael T., Wu Yue, Tayyari Fariba, Hattori Ayuna, Glushka John, Ito Takahiro, Arnold Jonathan, Edison Arthur S. 2019. Continuous in vivo Metabolism by NMR. Frontiers in Molecular Biosciences. Vol. 6
%% the tracing process is semi-automatic
%% it is based on preset ppm region.
%% this script select complex/representative region in the real spectral as there are possible too many region to deal with
set(0,'DefaultFigureVisible','on');
close all;
clear all;
comp='/Users/yuewu/';%the computer user location
workdir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge_tracing_manuscript/analysis_res/result/experiment_quantification_manul_compare/'];%%the working folder
cd(workdir);
load('../../data/sampleData.mat');
path=workdir;
%% visual parameter
horzshift=0.001;
vertshift=1E-3;
%% stack plot show
samp_i=1;
stackSpectra(sampleData(samp_i).Xcollapsed_1h1d,sampleData(samp_i).ppm_1h1d,horzshift,1E-2,'N. crassa NMR metabolic fingerprint over time')
%% selected region to work on
% for ridi=172:178
%   ppmvec=Samples2(1).adjustedRidges(ridi).ppms;
%   [min(ppmvec) max(ppmvec)]
% end

%% [range] name
%% [2.2 2.22] 'unknown-179'
%% [0.97 0.99] 'valine 1' 'valine 2'
%% [1.022 1.042]  'valine 3' 'valine 4'
%% [0.992 1.01] 'isoleucine 1' 'isoleucine 2'
%% [1.45 1.49] 'alanine 1' 'alanine 2'
%% [1.31 1.34] 'lactate 1' 'lactate 2'
%% [1.61 1.67] 'arginine 1' 'arginine 2' 'arginine 3' 'arginine 4' 'arginine 5' 'arginine 6' 'arginine 7' 'arginine 8' # not the second shoulder ridge
%% [5.785 5.81] 'uracil 1' 'uracil 2'
%% [5.45 5.5] 'glucose-1-phosphate 1' 'glucose-1-phosphate 2' 'glucose-1-phosphate 3' 'glucose-1-phosphate 4'
%% [8.25 8.36] 'adenosine 3' 'adenosine 4'
%% [6.87 6.91] 'tyrosine 1' 'tyrosine 2'
%% [7.4 7.43] 'phenylalanine 1' 'phenylalanine 2'
%% [7.85 7.875] 'uridine 1' 'uridine 2'
%% [3.15 3.2] 'choline'
%% [5.17 5.195] 'trehalose 1' 'trehalose 2'
%% [7.95 8.05] 'guanosine' (left)
%% [2.55 2.8] 'citrate 1' 'citrate 2' 'citrate 3' 'citrate 4'
%% [2.4 2.55] 'succinate'
%% [2.32 2.42] 'glutamate' the left side peak of the six peak
%% [6.5 6.6] 'fumarate'
%% [8.42 8.45] 'formate'
%% [5.2 5.26] 'glucose 7' 'glucose 8' [4.6 4.67] 'glucose 1' 'glucose 2' [3.87 3.92] 'glucose 3' 'glucose 4'  'glucose 5' 'glucose 6'
%% [3.62 3.68] 'ethanol 1' 'ethanol 2' 'ethanol 3' 'ethanol 4' [1.15 1.2] 'ethanol 5' 'ethanol 6' 'ethanol 7'
regionsele=[2.2 2.22; 0.97 0.99; 1.022 1.042; 0.992 1.01; 1.45 1.49; 1.31 1.34; 1.61 1.67; 5.785 5.81; 5.45 5.5; 8.25 8.36; 6.87 6.91; 7.405 7.43; 7.85 7.875; 3.15 3.2; 5.17 5.2; 7.95 8.05; 2.55 2.8; 2.4 2.55; 2.32 2.42; 6.5 6.6; 8.42 8.45; 5.2 5.26; 4.6 4.67; 3.87 3.92; 3.62 3.68; 1.15 1.2];
compoundList={'unknown-179','valine','isoleucine','alanine','lactate','arginine','uracil','glucose-1-phosphate','adenosine','tyrosine','phenylalanine','uridine','choline','trehalose','guanosine','citrate','succinate','glutamate','fumarate','formate','glucose','ethanol'};
samples=[1,2,3,7,8,9];
sampleKey={'aerobic','aerobic','aerobic','anaerobic','anaerobic','anaerobic'};
% preparation of data structure
Sample=[];
Sample(1).ridges(1).parameters=[];
Sample(1).ridges(1).result=[];
for i = 2:length(samples)
    Sample(i).ridges(1)=Sample(1).ridges(1);
end
defaultinputstart=struct();
defaultinputstart.compd='unknown';
defaultinputstart.quan='N';
%%%%%%%% test run %%%%%%%%%%%%%%
thredseg=10; %the main tune parameter default 10
maxaddon=1;%% some times need changes for wavy peaks in intensity, default 1
samp_i=1;
i=22;
regionhere=regionsele(i,:);
data=sampleData(samp_i);
mat=data.Xcollapsed_1h1d;
ppm=data.ppm_1h1d;
time=data.timesCollapsed_1h1d;
[returndata]=ridgetrace_power2_ext(mat,ppm,time,regionhere,path,thredseg,maxaddon,'defaultinput',defaultinputstart);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the production run
thredseglist=repmat(10,[1 size(regionsele,1)]);% default 10
thredseglist(17)=15;%citrate
thredseglist(18)=15;%succinate
maxaddonlist=repmat(1,[1 size(regionsele,1)]);%default 1
maxaddonlist(6)=2;%Lactate&Threonine
maxaddonlist(15)=2;%trehalose
maxaddonlist(17)=4;%citrate
maxaddonlist(19)=6;%glutamate
maxaddonlist(25)=4;%ethanol
maxaddonlist(26)=3;%ethanol
defaultinput=defaultinputstart;
for i = 1:size(regionsele,1)
  % showfigtitle=[];
  for samp_i_i=1:length(samples)
    % if samp_i~=1
    %   showfig=openfig(strcat(path,showfigtitle,'.surf.experiment.manual.fig'));
    % end
    samp_i=samples(samp_i_i);
    data=sampleData(samp_i);
    mat=data.Xcollapsed_1h1d;
    ppm=data.ppm_1h1d;
    time=data.timesCollapsed_1h1d;
    disp(['sample ' num2str(samp_i) ' region ' num2str(i)]);
    regionhere=regionsele(i,:);
    plotTitle=[num2str(regionhere(1)),'.',num2str(regionhere(2)),'.',num2str(samp_i),'.testplot'];
    [returndata]=ridgetrace_power2_ext(mat,ppm,time,regionhere,path,thredseglist(i),maxaddonlist(i),'defaultinput',defaultinput);
    defaultinput.compd=returndata.names{1};
    defaultinput.quan=returndata.quantifyvec{1};
    fig=gcf;
    saveas(fig,strcat(path,plotTitle,'.surf.experiment.manual.fig'));
    if samp_i==1
      showfigtitle=plotTitle;
    end
    close(fig);
    result=returndata.result;
    ridnames=returndata.names;
    quanvec=returndata.quantifyvec;
    groups=result(:,5);
    %% store as a struct table
    tempstorerids=[];
    tempstorerids.ridges(1).parameters=[];
    tempstorerids.ridges(1).result=[];
    for group=unique(groups,'stable')'
      groupind=find(result(:,5)==group);
      temptab=result(groupind,:);
      temptab(:,5)=[];
      resdata=struct();
      resdata.linearind=temptab(:,1);
      resdata.colind=temptab(:,2);
      resdata.rowind=temptab(:,3);
      resdata.intensity=temptab(:,4);
      resdata.names=ridnames(groupind);
      resdata.quanvec=quanvec(groupind);
      resdata.time=temptab(:,5);
      resdata.ppm=temptab(:,6);
      res=struct();
      res.parameters=returndata.para;
      res.result=resdata;
      Sample(samp_i_i).ridges=[Sample(samp_i_i).ridges res];
      tempstorerids.ridges=[tempstorerids.ridges res];
    end
    %% the structure for ploting
    peakshere=struct();
    ridgenumbers=1:length(tempstorerids.ridges);
    for scali=1:length(ridgenumbers)
      if ~isempty(tempstorerids.ridges(ridgenumbers(scali)).result)
        % peakshere(scali).Ridges=ppm(tempstorerids.ridges(ridgenumbers(scali)).result.colind);
        peakshere(scali).Ridges=tempstorerids.ridges(ridgenumbers(scali)).result.ppm';
        peakshere(scali).RowInds=tempstorerids.ridges(ridgenumbers(scali)).result.rowind';
        peakshere(scali).RidgeIntensities=tempstorerids.ridges(ridgenumbers(scali)).result.intensity';
        peakshere(scali).CompoundNames=tempstorerids.ridges(ridgenumbers(scali)).result.names;
        peakshere(scali).quantifiable=tempstorerids.ridges(ridgenumbers(scali)).result.quanvec;
      else
        peakshere(scali).Ridges=[];
        peakshere(scali).RowInds=[];
        peakshere(scali).RidgeIntensities=[];
        peakshere(scali).CompoundNames=[];
        peakshere(scali).quantifiable=[];
      end
    end
    reg=matchPPMs(regionhere,ppm);
    ind=reg(1):reg(2);
    mathere=mat(:,ind);
    ppmhere=ppm(ind);
    fig=stackSpectra_paintRidges_3return(mathere,ppmhere,horzshift,0.01,plotTitle,peakshere,10);
    saveas(fig,strcat(path,plotTitle,'.scatter.experiment.manual.fig'));
    close(fig);
  end
  % close(showfig);
end
%%plotting
for samp_i_i=1:length(samples)
  samp_i=samples(samp_i_i);
  data=sampleData(samp_i);
  mat=data.Xcollapsed_1h1d;
  ppm=data.ppm_1h1d;
  peaks=struct();
  ridgenumbers=1:length(Sample(samp_i_i).ridges);
  for i=1:length(ridgenumbers)
    if ~isempty(Sample(samp_i_i).ridges(ridgenumbers(i)).result)
      % peaks(i).Ridges=ppm(Sample(samp_i_i).ridges(ridgenumbers(i)).result.colind);
      peaks(i).Ridges=Sample(samp_i_i).ridges(ridgenumbers(i)).result.ppm';
      peaks(i).RowInds=Sample(samp_i_i).ridges(ridgenumbers(i)).result.rowind';
      peaks(i).RidgeIntensities=Sample(samp_i_i).ridges(ridgenumbers(i)).result.intensity';
      peaks(i).CompoundNames=Sample(samp_i_i).ridges(ridgenumbers(i)).result.names;
      peaks(i).quantifiable=Sample(samp_i_i).ridges(ridgenumbers(i)).result.quanvec;
    else
      peaks(i).Ridges=[];
      peaks(i).RowInds=[];
      peaks(i).RidgeIntensities=[];
      peaks(i).CompoundNames=[];
      peaks(i).quantifiable=[];
    end
  end
  plotTitle=[num2str(samp_i) '.' 'region'];
  fig=stackSpectra_paintRidges_3return(mat,ppm,0.0,0.02,plotTitle,peaks,10);
  saveas(fig,strcat(path,plotTitle,'.scatter.experiment.manual.fig'));
  close(fig);
end
save([workdir 'tracing.newmeth.experiment.manual.mat'],'Sample');

%% conctruct compound structure for plotting
load([workdir 'tracing.newmeth.experiment.manual.mat']);
lags=[4   0.013194444 1
        5   0.011805556 2
        6   0.011805556 3
        10  0.018055556 7
        11  0.01875 8
        12  0.015277778 9];
lags(:,2)=lags(:,2)*24; % convert days to hours
% For each sample correct time
for s=1:length(Sample)
  for r=2:length(Sample(s).ridges)
    Sample(s).ridges(r).result.time=Sample(s).ridges(r).result.time+lags(s,2);
  end
end
[compounds]=combineRidges_new(Sample,compoundList,samples,sampleKey);
save([workdir 'mapped.newmeth.experiment.manual.mat'],'compounds');

%% Plot the trajectories of the different compounds as a function of time
%%% to compare with published paper
load([workdir 'tracing.newmeth.experiment.manual.mat']);
mkdir('AllPlots_updated_scaled')
cd('AllPlots_updated_scaled')
blue=[0,.3,1];
yellow=[1,0.7490,0.0510];%
red=[0.8,0,0];
% Add to the compounds structure
colorlist={red,red,red,blue,blue,blue};
conditionlist={'aerobic','aerobic','aerobic','anaerobic','anaerobic','anaerobic'};
for c=1:length(compoundList)
    figure('PaperType','<custom>','PaperSize',[8 6],'Color',[1 1 1]),hold on,
    % Get the relevant rows
    totalintensity=[];
    for s=1:length(Sample)
        currentRidges=compounds(find(strcmp({compounds.Name},compoundList{c})&[compounds.SampleNumber]==samples(s)),:);
        if length(currentRidges)==0
          continue;
        end
        times=currentRidges.AverageTrajectory_times;
        intensity=currentRidges.AverageTrajectory_intensities_scaled;
        [timesort index]=sort(times);
        intensitysort=intensity(index);
        h(s)=plot(timesort,intensitysort,...
            ...%'DisplayName',[currentRidges(s).Condition,' replicate'],...
            'Color',colorlist{s},...
            'LineWidth',1,...
            'Marker','o',...
            'MarkerSize',3);
            set(h(s),'DisplayName',[conditionlist{s},' replicate'])
        totalintensity=[totalintensity intensity];
    end
    %legend([h(conditionInds(1)),h(conditionInds(2))], conditions{1}, conditions{2})
    %xlabel('Time (h)')
    %ylabel({'Scaled Ridge Intensity','(time-wise mean)'})
    %title(['Peak Traces in an ','\itN. crassa\rm\bf culture '])
    %title([compoundList{c},' (lag)'])
    % title(compoundList{c})
    set(gcf,'InvertHardCopy', 'off');
    set(gca,'fontsize',40)
    set(gca,'XLim',[0,13])
    maxValue=max(totalintensity(:));
    minValue=min(totalintensity(:));
    set(gca,'YLim',[minValue,maxValue])
    title(compoundList{c})
    % Get rid of annoying whitespace on the outsides
    fig=gca;
    InSet=get(fig,'TightInset');
    set(fig,'Position',[InSet(1:2),1-InSet(1)-InSet(3),1-InSet(2)-InSet(4)])
    % Programmatically save as .fig and .pdf by title
    saveas(fig,fig.Title.String)
    print(['-f',num2str(get(gcf,'Number'))],fig.Title.String,'-dpdf')
end
cd ..

%%residue plot
compoundslist=struct();
load([workdir 'mapped.newmeth.experiment.manual.mat']);
compoundslist.newmeth=compounds;
load('/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/archive/test.code/analysis/multi_sample_data/compounds.mat');%%%if you want to run this line of code, please read our previous publication and download corresponding files
compoundslist.oldmeth=compounds;
mkdir('residue_plot');
cd('residue_plot');
colorlist={red,red,red,blue,blue,blue};
conditionlist={'aerobic','aerobic','aerobic','anaerobic','anaerobic','anaerobic'};
for c=1:length(compoundList)
    figure('PaperType','<custom>','PaperSize',[8 6],'Color',[1 1 1]),hold on,
    % Get the relevant rows
    totalintensity=[];
    totalresidue=[];
    for s=1:length(Sample)
        currentRidgeslist=struct();
        fieldnamevec=fieldnames(compoundslist);
        for fieldnam_i=1:length(fieldnamevec)
          fieldnam=fieldnamevec{fieldnam_i};
          localdata=compoundslist.(fieldnam);
          currentRidgeslist.(fieldnam)=localdata(find(strcmp({compoundslist.(fieldnam).Name},compoundList{c})&[compoundslist.(fieldnam).SampleNumber]==samples(s)),:);
        end
        if length(currentRidgeslist.newmeth)==0||length(currentRidgeslist.oldmeth)==0
          continue;
        end
        [times,iold,inew]=intersect(currentRidgeslist.oldmeth.AverageTrajectory_times,currentRidgeslist.newmeth.AverageTrajectory_times);
        intensityold=currentRidgeslist.oldmeth.AverageTrajectory_intensities_scaled(iold);
        intensitynew=currentRidgeslist.newmeth.AverageTrajectory_intensities_scaled(inew);
        residues=intensitynew-intensityold;
        [timesort index]=sort(times);
        residuessort=residues(index);
        h(s)=plot(timesort,residuessort,...
            ...%'DisplayName',[currentRidges(s).Condition,' replicate'],...
            'Color',colorlist{s},...
            'LineWidth',1,...
            'Marker','o',...
            'MarkerSize',3);
            set(h(s),'DisplayName',[conditionlist{s},' replicate'])
        totalintensity=[totalintensity intensityold];
        totalresidue=[totalresidue residuessort];
    end
    set(gcf,'InvertHardCopy', 'off');
    set(gca,'fontsize',40)
    set(gca,'XLim',[0,13])
    maxValue=max(totalintensity(:));
    minValue=min(totalresidue(:));
    set(gca,'YLim',[minValue,maxValue])
    title([compoundList{c} ' residues'])
    % Get rid of annoying whitespace on the outsides
    fig=gca;
    InSet=get(fig,'TightInset');
    set(fig,'Position',[InSet(1:2),1-InSet(1)-InSet(3),1-InSet(2)-InSet(4)])
    % Programmatically save as .fig and .pdf by title
    saveas(fig,fig.Title.String)
    print(['-f',num2str(get(gcf,'Number'))],fig.Title.String,'-dpdf')
end
cd ..

%%residue plot
compoundslist=struct();
load([workdir 'tracing.newmeth.experiment.manual.mat']);
load([workdir 'mapped.newmeth.experiment.manual.mat']);
compoundslist.newmeth=compounds;
load('/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/archive/test.code/analysis/multi_sample_data/compounds.mat');%%%if you want to run this line of code, please read our previous publication and download corresponding files
compoundslist.oldmeth=compounds;
mkdir('residue_plot_2');
cd('residue_plot_2');
blue=[0,.3,1];
yellow=[1,0.7490,0.0510];%
red=[0.8,0,0];
colorlist={red,red,red,blue,blue,blue};
conditionlist={'aerobic','aerobic','aerobic','anaerobic','anaerobic','anaerobic'};
totalintvec=[];
for c=1:length(compoundList)
    % Get the relevant rows
    totalintensity=[];
    totalresidue=[];
    for s=1:length(Sample)
        currentRidgeslist=struct();
        fieldnamevec=fieldnames(compoundslist);
        for fieldnam_i=1:length(fieldnamevec)
          fieldnam=fieldnamevec{fieldnam_i};
          localdata=compoundslist.(fieldnam);
          currentRidgeslist.(fieldnam)=localdata(find(strcmp({compoundslist.(fieldnam).Name},compoundList{c})&[compoundslist.(fieldnam).SampleNumber]==samples(s)),:);
        end
        if length(currentRidgeslist.newmeth)==0||length(currentRidgeslist.oldmeth)==0
          continue;
        end
        [times,iold,inew]=intersect(currentRidgeslist.oldmeth.AverageTrajectory_times,currentRidgeslist.newmeth.AverageTrajectory_times);
        intensityold=currentRidgeslist.oldmeth.AverageTrajectory_intensities_scaled(iold);
        intensitynew=currentRidgeslist.newmeth.AverageTrajectory_intensities_scaled(inew);
        residues=intensitynew-intensityold;
        [timesort index]=sort(times);
        residuessort=residues(index);
        totalintvec=[totalintvec residuessort];
    end
end
for c=1:length(compoundList)
    figure('PaperType','<custom>','PaperSize',[8 6],'Color',[1 1 1]),hold on,
    % Get the relevant rows
    totalintensity=[];
    totalresidue=[];
    for s=1:length(Sample)
        currentRidgeslist=struct();
        fieldnamevec=fieldnames(compoundslist);
        for fieldnam_i=1:length(fieldnamevec)
          fieldnam=fieldnamevec{fieldnam_i};
          localdata=compoundslist.(fieldnam);
          currentRidgeslist.(fieldnam)=localdata(find(strcmp({compoundslist.(fieldnam).Name},compoundList{c})&[compoundslist.(fieldnam).SampleNumber]==samples(s)),:);
        end
        if length(currentRidgeslist.newmeth)==0||length(currentRidgeslist.oldmeth)==0
          continue;
        end
        [times,iold,inew]=intersect(currentRidgeslist.oldmeth.AverageTrajectory_times,currentRidgeslist.newmeth.AverageTrajectory_times);
        intensityold=currentRidgeslist.oldmeth.AverageTrajectory_intensities_scaled(iold);
        intensitynew=currentRidgeslist.newmeth.AverageTrajectory_intensities_scaled(inew);
        residues=intensitynew-intensityold;
        [timesort index]=sort(times);
        residuessort=residues(index);
        h(s)=plot(timesort,residuessort,...
            ...%'DisplayName',[currentRidges(s).Condition,' replicate'],...
            'Color',colorlist{s},...
            'LineWidth',1,...
            'Marker','o',...
            'MarkerSize',3);
            set(h(s),'DisplayName',[conditionlist{s},' replicate'])
        totalintensity=[totalintensity intensityold];
        totalresidue=[totalresidue residuessort];
    end
    set(gcf,'InvertHardCopy', 'off');
    set(gca,'fontsize',40)
    set(gca,'XLim',[0,13])
    maxValue=max(totalintvec(:))+1;
    minValue=min(totalintvec(:));
    set(gca,'YLim',[minValue,maxValue])
    title([compoundList{c} ' residues 2'])
    % Get rid of annoying whitespace on the outsides
    fig=gca;
    InSet=get(fig,'TightInset');
    set(fig,'Position',[InSet(1:2),1-InSet(1)-InSet(3),1-InSet(2)-InSet(4)])
    % Programmatically save as .fig and .pdf by title
    saveas(fig,fig.Title.String)
    print(['-f',num2str(get(gcf,'Number'))],fig.Title.String,'-dpdf')
end
cd ..

%%%test and check
%%%%spectra visulization and compare between simulated data and real data
% samp_i=3;
% mat=matrixstr{samp_i};
% plotr(ppm,mat(1,:))

%%check scatter plotting result
% dirfig=workdir;
% cd(dirfig);
% i=14;
% regionhere=regionsele(i,:);
% posivec={'northwest' 'northeast'};
% for sample_i_i=1:length(sample)
%   sample_i=sample(sample_i_i);
%   file=[dirfig num2str(regionhere(1)) '.' num2str(regionhere(2)) '.' num2str(sample_i) '.testplot.scatter.experiment.manual.fig']
%   fig=openfig(file);
%   movegui(fig,posivec{sample_i_i});
% end
% close all;
