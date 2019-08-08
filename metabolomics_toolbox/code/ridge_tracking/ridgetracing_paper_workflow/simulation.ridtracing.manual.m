%% This script is to trace ridges from simulated data
%% The tracing process is semi-automatic
%% It is based on preset ppm region.
%%% For ridge tracing based on simulated data, return informaiton on compound name and quantifiablity is not relevant/used
%%% This is because ridge match is not done in this script but by match_and_output.m script for 'annotation'.
%%% This design is for automation purpose.
%%% instead, in ridge tracing in real data, the return informaiton in compound name and quantifiablity is informative as the user input that information based on their annotation
close all;
clear all;
comp='/Users/yuewu/';%the computer user location
workdir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge_tracing_manuscript/analysis_res/result/simulated_quantification_manual/'];%%the working folder
cd(workdir);
load('../simulateddata/simulated.timeseries.ph.complex.more.mat');
path=workdir;
sample=1:4;
%% visual parameter
horzshift= -0.0025;
vertshift= 1E-3;
%% visual check the sample spectra to select region to work on
stackSpectra(matrixstr{4},ppm,0.0,1,'N. crassa NMR metabolic fingerprint over time')
%% selected region to work on
regionsele=[-0.03 0.03; 0.6 0.65; 0.8 1.1; 1.15 1.185; 1.29 1.37; 1.4 1.6; 1.64 1.8; 1.9 2.1; ...
            2.22 2.31; 2.35 2.4; 2.88 2.94; 3.18 3.25; 3.3 3.5; 3.5 3.7; 3.7 3.85; ...
            3.87 4.0; 4.1 4.15; 4.18 4.28; 4.3 4.4; 4.75 4.95; 5.8 6; 7.8 7.9; 8.3 8.6; 8.8 9.1];
%complex regions: *[1.64 1.8] *[1.9 2.1] *[2.22 2.31]
% *[3.3 3.5] *[3.7 3.85]
%  *[4.18 4.28] *[4.75 4.95]
% preparation of data structure
Sample=[];
Sample(1).ridges(1).parameters=[];
Sample(1).ridges(1).result=[];
for i = 2:length(sample)
    Sample(i).ridges(1)=Sample(1).ridges(1);
end
defaultinputstart=struct();
defaultinputstart.compd='unknown';
defaultinputstart.quan='N';
%%%%%%%% test run %%%%%%%%%%%%%%
%%%%%the two main tuning parameter
thredseg=10; %% the more curvy the ridges the higher this value
maxaddon=1; %% if find out ridges not traced, increase this number to the number of ridges in the region
%%%%
samp_i=3;
i=24;
regionhere=regionsele(i,:);
mat=matrixstr{samp_i};
[returndata]=ridgetrace_power2_ext(mat,ppm,times',regionhere,path,thredseg,maxaddon,'defaultinput',defaultinputstart);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the production run
thredseglist=repmat(10,[1 size(regionsele,1)]);% default 10
thredseglist(8)=50;%[1.9 2.1]
thredseglist(13)=50;%[3.3 3.5]
thredseglist(18)=50;%[4.18 4.28]
thredseglist(20)=50;%[4.75 4.95]
thredseglist(23)=50;%[8.3 8.6]
maxaddonlist=repmat(1,[1 size(regionsele,1)]);%default 1
maxaddonlist(2)=3;%[0.6 0.65]
maxaddonlist(7)=10;%[1.64 1.8]
maxaddonlist(9)=7;%[2.22 2.31]
maxaddonlist(11)=3;%[2.88 2.94]
defaultinput=defaultinputstart;
for i = 1:size(regionsele,1)
  for samp_i=sample
    mat=matrixstr{samp_i};
    disp(['sample ' num2str(samp_i) ' region ' num2str(i)]);
    regionhere=regionsele(i,:);
    plotTitle=[num2str(regionhere(1)),'.',num2str(regionhere(2)),'.',num2str(samp_i),'.testplot'];
    [returndata]=ridgetrace_power2_ext(mat,ppm,times',regionhere,path,thredseglist(i),maxaddonlist(i),'defaultinput',defaultinput);
    defaultinput.compd=returndata.names{1};
    defaultinput.quan=returndata.quantifyvec{1};
    fig=gcf;
    saveas(fig,strcat(path,plotTitle,'.surf.manual.fig'));
    close(fig);
    result=returndata.result;
    ridnames=returndata.names;
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
      resdata.time=temptab(:,5);
      resdata.ppm=temptab(:,6);
      res=struct();
      res.parameters=returndata.para;
      res.result=resdata;
      Sample(samp_i).ridges=[Sample(samp_i).ridges res];
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
      else
        peakshere(scali).Ridges=[];
        peakshere(scali).RowInds=[];
        peakshere(scali).RidgeIntensities=[];
        peakshere(scali).CompoundNames=[];
      end
    end
    reg=matchPPMs(regionhere,ppm);
    ind=reg(1):reg(2);
    mathere=mat(:,ind);
    ppmhere=ppm(ind);
    fig=stackSpectra_paintRidges_3return(mathere,ppmhere,0.0,0.01,plotTitle,peakshere,10);
    saveas(fig,strcat(path,plotTitle,'.scatter.manual.fig'));
    close(fig);
  end
end
%%plotting
for samp_i=sample
  mat=matrixstr{samp_i};
  peaks=struct();
  ridgenumbers=1:length(Sample(samp_i).ridges);
  for i=1:length(ridgenumbers)
    if ~isempty(Sample(samp_i).ridges(ridgenumbers(i)).result)
      % peaks(i).Ridges=ppm(Sample(samp_i).ridges(ridgenumbers(i)).result.colind);
      peaks(i).Ridges=Sample(samp_i).ridges(ridgenumbers(i)).result.ppm';
      peaks(i).RowInds=Sample(samp_i).ridges(ridgenumbers(i)).result.rowind';
      peaks(i).RidgeIntensities=Sample(samp_i).ridges(ridgenumbers(i)).result.intensity';
      peaks(i).CompoundNames=Sample(samp_i).ridges(ridgenumbers(i)).result.names;
    else
      peaks(i).Ridges=[];
      peaks(i).RowInds=[];
      peaks(i).RidgeIntensities=[];
      peaks(i).CompoundNames=[];
    end
  end
  plotTitle=[num2str(samp_i) '.' 'region'];
  fig=stackSpectra_paintRidges_3return(mat,ppm,0.0,0.02,plotTitle,peaks,10);
  saveas(fig,strcat(path,plotTitle,'.scatter.manual.fig'));
  close(fig);
end
save([workdir '/tracing.newmeth.simulated.manual.mat'],'Sample')

%%%test and check
%%%%spectra visulization and compare between simulated data and real data
% samp_i=3;
% mat=matrixstr{samp_i};
% plotr(ppm,mat(1,:))

%%check scatter plotting result
% dirfig=workdir;
% cd(dirfig);
% i=1;
% regionhere=regionsele(i,:);
% posivec={'northwest' 'southwest' 'northeast' 'southeast'};
% for j=1:4
%   file=[dirfig num2str(regionhere(1)) '.' num2str(regionhere(2)) '.' num2str(j) '.testplot.scatter.manual.fig']
%   fig=openfig(file);
%   movegui(fig,posivec{j});
% end
% close all;
