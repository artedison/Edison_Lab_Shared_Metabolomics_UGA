%% this script is to trace ridges from experimental data
%% the tracing process is semi-automatic
%% it is based on preset ppm region.
%% this script select complex/representative region in the real spectral as there are possible too many region to deal with
close all;
clear all;
%%%%%%%%%%%%%%%%%%%%path and data loading%%%%%%%%%%%%%%%%%%%%
%%the user can cd to their own directory and load sampleData.mat from corresponding location
% comp='/Users/yuewu/';%the computer user location
% workdir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge_tracing_manuscript/analysis_res/result/experiment_quantification_manul/'];%%the working folder
% cd(workdir);
% load('../../data/sampleData.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path=pwd;
sample=1;%the plotting sample%[1 4];
%% visual parameter
horzshift= 0;
vertshift= 1E-3;
%% stack plot show
samp_i=1;
stackSpectra(sampleData(samp_i).Xcollapsed_1h1d,sampleData(samp_i).ppm_1h1d,horzshift,vertshift,'N. crassa NMR metabolic fingerprint over time')%%visualize time series spectra
%% selected region to work on
% allregionsele=[-0.05 0.05; 0.6 0.7; 0.8 0.9; 0.9 0.955; 0.955 0.97; 0.97 0.99; 0.992 1.05; 1.05 1.06; 1.06 1.07; 1.075 1.1; 1.1 1.135; 1.135 1.144; 1.144 1.15; 1.15 1.2; 1.25 1.3; 1.3 1.35; 1.35 1.4; 1.46 1.49; 1.5 1.6; 1.6 1.625; 1.625 1.67; 1.67 1.7; 1.7 1.8; 1.85 1.92; 1.90 2.0; 2.0 2.09; 2.07 2.17; 2.2 2.26; 2.3 2.365; 2.33 2.44; 2.4 2.6; 2.5 2.9; 2.79 2.83; 2.83 2.9; 2.9 2.95; 2.97 3.03; 3.03 3.1; 3.1 3.15; 3.15 3.3; 3.3 3.36; 3.36 3.43; 3.43 3.5; 3.5 3.54; 3.545 3.585; 3.585 3.62; 3.62 3.685; 3.68 3.8; 3.8 3.855; 3.855 3.875; 3.875 3.91; 3.92 3.96; 3.96 4.0; 4.0 4.08; 4.08 4.155; 4.06 4.19; 4.19 4.24; 4.24 4.29; 4.28 4.35; 4.38 4.41; 4.41 4.45; 4.48 4.54; 4.6 4.68; 5.15 5.2; 5.2 5.25; 5.25 5.35; 5.36 5.41; 5.46 5.5; 5.78 5.81; 5.86 5.92; 5.95 6.0; 6.05 6.08; 6.04 6.25; 6.51 6.55; 6.87 6.9; 7.15 7.2; 7.25 7.4; 7.4 7.44; 7.5 7.55; 7.8 8.0; 8.18 8.24; 8.25 8.36; 8.4 8.45; 8.5 8.62; 9.65 9.7];
regionsele=[1.3 1.35; 2.33 2.44; 3.15 3.3; 3.68 3.8];%the plotting region
% [1.1 1.135; 1.3 1.35; 1.85 1.92; 2.00 2.11; 2.33 2.44; 3.1 3.15; 3.15 3.3; 3.68 3.8; 4.08 4.155; 4.28 4.35; 4.48 4.54; 5.36 5.41; 6.05 6.08; 6.51 6.55];
%% complex region to compare [1.1 1.135] [1.3 1.35] [2.00 2.11] [2.33 2.44]
%% [3.1 3.15] [3.15 3.3] [3.68 3.8]
% preparation of data structure
Sample=[];
Sample(1).ridges(1).parameters=[];
Sample(1).ridges(1).result=[];
for i = 2:length(sample)
    Sample(i).ridges(1)=Sample(1).ridges(1);
end
%%%%%%%% test run %%%%%%%%%%%%%%
thredseg=10; %the main tune parameter default 10
maxaddon=6;%% some times need changes for wavy peaks in intensity, default 1
samp_i=1;
i=2;
regionhere=regionsele(i,:);
data=sampleData(samp_i);
mat=data.Xcollapsed_1h1d;
ppm=data.ppm_1h1d;
time=data.timesCollapsed_1h1d;
[returndata]=ridgetrace_power2_ext(mat,ppm,time,regionhere,path,thredseg,maxaddon);%%ridge tracking function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the production run
thredseglist=repmat(10,[1 size(regionsele,1)]);% default 10
thredseglist(3)=5;
thredseglist(4)=2;
% thredseglist(4)=5;
% thredseglist(7)=5;
% thredseglist(8)=2;
% thredseglist(9)=5;
maxaddonlist=repmat(1,[1 size(regionsele,1)]);%default 1
maxaddonlist([1])=2;
maxaddonlist([2])=6;
% maxaddonlist([2])=2;
% maxaddonlist([5])=6;
% maxaddonlist([6])=3;
for i = 1:size(regionsele,1)
  % showfigtitle=[];
  for samp_i_i=1:length(sample)
    % if samp_i~=1
    %   showfig=openfig(strcat(path,showfigtitle,'.surf.experiment.manual.fig'));
    % end
    samp_i=sample(samp_i_i);
    data=sampleData(samp_i);
    mat=data.Xcollapsed_1h1d;
    ppm=data.ppm_1h1d;
    time=data.timesCollapsed_1h1d;
    disp(['sample ' num2str(samp_i) ' region ' num2str(i)]);
    regionhere=regionsele(i,:);
    plotTitle=[num2str(regionhere(1)),'.',num2str(regionhere(2)),'.',num2str(samp_i),'.testplot'];
    [returndata]=ridgetrace_power2_ext(mat,ppm,time,regionhere,path,thredseglist(i),maxaddonlist(i));
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
for samp_i_i=1:length(sample)
  samp_i=sample(samp_i_i);
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
save([workdir 'tracing.newmeth.experiment.manual.mat'],'Sample')

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

%%test plot for valine region [1.3 1.35] 0.000,0.02 [2.33 2.44] 0.0004,0.04 [3.15 3.3] 0.0004,0.04 [3.68 3.8] 0.0004,0.04
% close all;
% clear all;
% ppmregion=[3.68 3.8];%
% comp='/Users/yuewu/';%the computer user location
% workdir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge_tracing_manuscript/analysis_res/result/experiment_quantification_manul/'];%%the working folder
% cd(workdir);
% % load([workdir 'tracing.newmeth.experiment.manual.mat']);
% % load('../../data/sampleData.mat');
% i=2;
% j=4;
% compname='region2';
% peakcell=Sample(i).ridges;
% matchcomppeak=length(peakcell);%number of peaks
% peaks=struct();
% ridgenumbers=2:length(peakcell);
% data=sampleData(j);
% mat=data.Xcollapsed_1h1d;
% ppm=data.ppm_1h1d;
% jm=1;
% peakshere=struct();
% for j=ridgenumbers
%   temppeak=peakcell(j).result;
%   ppmvec=temppeak.ppm';
%   if min(ppmvec)<ppmregion(1)|max(ppmvec)>ppmregion(2)
%     continue;
%   end
%   if ~isempty(temppeak)
%     peakshere(jm).Ridges=temppeak.ppm';
%     peakshere(jm).RowInds=temppeak.rowind';
%     peakshere(jm).RidgeIntensities=temppeak.intensity';
%     peakshere(jm).CompoundNames=temppeak.names;
%     peakshere(jm).quantifiable=temppeak.quanvec;
%   else
%     peakshere(jm).Ridges=[];
%     peakshere(jm).RowInds=[];
%     peakshere(jm).RidgeIntensities=[];
%     peakshere(jm).CompoundNames=[];
%     peakshere(jm).quantifiable=[];
%   end
%   jm=jm+1;
% end
% ppmrang=matchPPMs(ppmregion,ppm);
% ppmind=ppmrang(1):ppmrang(2);
% mathere=mat(:,ppmind);
% ppmrere=ppm(ppmind);
% namestr=[compname '.' num2str(ppmregion(1)) '.' num2str(ppmregion(2))];
% fig=stackSpectra_paintRidges_3return(mathere,ppmrere,0.0004,0.04,namestr,peakshere,10);
% saveas(fig,strcat(workdir,[namestr '.fig']));
% print(fig,[namestr '.pdf'],'-dpdf','-fillpage','-r2400','-painters');%'-painters','-fillpage',
% close(fig);
