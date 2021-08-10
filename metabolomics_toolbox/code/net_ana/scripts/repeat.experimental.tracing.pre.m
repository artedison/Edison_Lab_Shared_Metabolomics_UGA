%% this script is to trace ridges from experimental data
%% this tracking aimed to cover most trackable ridges of the 6 experimental datasets.
%% the tracing process is semi-automatic (need manual input)
%%
%% The data for reproduction is provided starting from [Combined tracking data]
%% However you are welcome to try the tracking procedure and report to Yue.Wu@uga.edu if you have any questions. Tracking the whole ~300 ridges for 6 replicates will take at least one day. For more information and tutorial on ridge tracking, please refer to
%% 1. the library: https://github.com/artedison/Edison_Lab_Shared_Metabolomics_UGA/tree/master/metabolomics_toolbox/code/ridge_tracking
%% 2. the publication:
% Yue Wu, Michael T Judge, Jonathan Arnold, Suchendra M Bhandarkar, Arthur S Edison, RTExtract: time-series NMR spectra quantification based on 3D surface ridge tracking, Bioinformatics, , btaa631, https://doi.org/10.1093/bioinformatics/btaa631
set(0,'DefaultFigureVisible','on');
close all;
clear all;
% please going through and change the following paths as needed
comp='/Users/yuewu/';%the computer user location
pardir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge.net/result_reprod/'];
% OR
% pardir=[pwd() '/'];
workdir=[pardir 'result/ridgetracking/'];%%the working folder
datadir=[pardir 'data/'];
continue_dir=[pardir 'result_data/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(workdir);
load([datadir 'unshared/sampleData.mat']);%please refer to the orginal CIVM-NMR paper for the data here
path=workdir;
%% visual parameter
horzshift=0;
vertshift=1E-3;
%% stack plot show
samp_i=1;
stackSpectra(sampleData(samp_i).Xcollapsed_1h1d,sampleData(samp_i).ppm_1h1d,horzshift,1E-2,'N. crassa NMR metabolic fingerprint over time');

%%the regions that have been tracked in previous result is ignored here
regionsele=[-0.05 0.07; 0.8 0.9; 0.9 0.965; 1.05 1.1; 1.1 1.135; 1.135 1.16; 1.2 1.25; 1.25 1.3; 1.35 1.4; 1.5 1.6; 1.67 1.7; 1.7 1.76; 1.85 1.875; 1.875 1.93; 1.93 2.0; 2.0 2.10; 2.09 2.17; 2.22 2.26; 2.32 2.34; 2.32 2.43; 2.4 2.46; 2.55 2.62; 2.79 2.828; 2.828 2.9; 2.9 2.95; 2.95 3.03; 3.03 3.1; 3.1 3.15; 3.2 3.28; 3.28 3.36; 3.36 3.43; 3.43 3.5; 3.5 3.54; 3.545 3.585; 3.585 3.63; 3.675 3.685; 3.685 3.7; 3.7 3.8; 3.8 3.855; 3.855 3.878; 3.92 3.96; 3.96 4.0; 4.0 4.08; 4.08 4.155; 4.06 4.19; 4.19 4.24; 4.24 4.29; 4.28 4.35; 4.35 4.38; 4.38 4.45; 4.48 4.54; 4.645 4.68; 5.05 5.1; 5.24 5.35; 5.36 5.41; 5.86 5.92; 5.95 6.0; 6.05 6.08; 6.04 6.25;  7.15 7.2; 7.25 7.4; 7.5 7.55; 7.873 7.983; 8.18 8.24; 8.45 8.62; 9.65 9.7];

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
thredseg=1; %the main tune parameter default 10
maxaddon=1;%% some times need changes for wavy peaks in intensity, default 1
for i=1:size(regionsele,1)
  for samp_i=[1 7]
    regionhere=regionsele(i,:);
    data=sampleData(samp_i);
    mat=data.Xcollapsed_1h1d;
    ppm=data.ppm_1h1d;
    time=data.timesCollapsed_1h1d;
    [returndata]=ridgetrace_power2_ext(mat,ppm,time,regionhere,path,thredseg,maxaddon,'defaultinput',defaultinputstart);%(28:52,:)
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the production run (need manual input)
thredseglist=repmat(1,[1 size(regionsele,1)]);%
thredseglist([10,17,23,24,28,30,43,44,47,48,49,51,53,55,57,58,61])=5;
thredseglist([20,45,64,59])=10;
thredseglist([65])=15;
maxaddonlist=repmat(1,[1 size(regionsele,1)]);%
maxaddonlist([26,28,30])=8;%
maxaddonlist([46,50,55])=5;%
timerang_mod=repmat([1 52],[size(regionsele,1),1]);%as some region need to be zoomed in to viausally check
timerang_mod(21,:)=[20,52]; %[2.4 2.46]
timerang_mod(22,:)=[20 52]; %[2.55 2.62]
timerang_mod(19,:)=[20 52]; %[2.32 2.34]
timerang_mod(52,:)=[20 52]; %[4.645 4.68]

defaultinput=defaultinputstart;
%%not tracking: the left most glutamate peak; the overlapped peak around 7.36. one peak around 4.29
%% for [4.645 4.68] sample 7-9 there might be trouble
% load([workdir 'tracing.newmeth.experiment.manual.mat'])
for i = 1:size(regionsele,1)
  % showfigtitle=[];
  for samp_i_i=1:length(samples)
    % if samp_i~=1
    %   showfig=openfig(strcat(path,showfigtitle,'.surf.experiment.manual.fig'));
    % end
    samp_i=samples(samp_i_i);
    data=sampleData(samp_i);
    loctimerang=timerang_mod(i,:);
    mat=data.Xcollapsed_1h1d;
    ppm=data.ppm_1h1d;
    time=data.timesCollapsed_1h1d;
    loctimeind=loctimerang(1):loctimerang(2);
    mat=mat(loctimeind,:);
    time=time(loctimeind,:);
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
      resdata.rowind=temptab(:,3)+loctimerang(1)-1;
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
save([workdir 'tracing.newmeth.experiment.manual.mat'],'Sample');

%%time subregion refining information update
% loading tracking data
% load([workdir 'tracing.newmeth.experiment.manual.mat'])
Sample2=Sample
subtime_ind=find(timerang_mod(:,1)~=1);
sub_ppmrange_tab=regionsele(subtime_ind,:);
subtimerange_tab=timerang_mod(subtime_ind,:);
for samp_i_i=1:length(samples)
  samp_i=samples(samp_i_i);
  data=sampleData(samp_i);
  mat=data.Xcollapsed_1h1d;
  wholesize=size(mat);
  ridgeslist=Sample2(samp_i_i).ridges;
  for ridge_i=2:size(ridgeslist,2)
    ridregion=ridgeslist(ridge_i).parameters.region;
    matchind=find(sub_ppmrange_tab(:,1)==ridregion(1)&sub_ppmrange_tab(:,2)==ridregion(2));
    if(length(matchind)~=0)
      ridgeslist(ridge_i).parameters.timregion_sub=subtimerange_tab(matchind,:);
      locridind=ridgeslist(ridge_i).result.rowind;
      matchind
      % if matchind~=4
      %   locridind=locridind+subtimerange_tab(matchind,1)-1;
      %   ridgeslist(ridge_i).result.rowind=locridind;
      % end
      ridgeslist(ridge_i).result.linearind=sub2ind(wholesize,locridind,ridgeslist(ridge_i).result.colind);
    end
  end
  Sample2(samp_i_i).ridges=ridgeslist;
end
save([workdir 'tracing.newmeth.experiment.manual.refined.mat'],'Sample2');

% loading saved data for continuing [Combined tracking data]
load([continue_dir 'tracing.newmeth.experiment.manual.refined.mat']);
%%plotting spectra
for samp_i_i=1:length(samples)
  samp_i=samples(samp_i_i);
  data=sampleData(samp_i);
  mat=data.Xcollapsed_1h1d;
  ppm=data.ppm_1h1d;
  peaks=struct();
  ridgenumbers=1:length(Sample2(samp_i_i).ridges);
  for i=1:length(ridgenumbers)
    if ~isempty(Sample2(samp_i_i).ridges(ridgenumbers(i)).result)
      % peaks(i).Ridges=ppm(Sample(samp_i_i).ridges(ridgenumbers(i)).result.colind);
      peaks(i).Ridges=Sample2(samp_i_i).ridges(ridgenumbers(i)).result.ppm';
      peaks(i).RowInds=Sample2(samp_i_i).ridges(ridgenumbers(i)).result.rowind';
      peaks(i).RidgeIntensities=Sample2(samp_i_i).ridges(ridgenumbers(i)).result.intensity';
      peaks(i).CompoundNames=Sample2(samp_i_i).ridges(ridgenumbers(i)).result.names;
      peaks(i).quantifiable=Sample2(samp_i_i).ridges(ridgenumbers(i)).result.quanvec;
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


%%combined with tracked data
load([datadir 'shared/tracing.newmeth.experiment.manual.mat'])
Sample_newrid=Sample2;
Sample_oldrid=Sample;
Sample_complete_rid=Sample_newrid;
for i=1:length(Sample_newrid)
  addtab=Sample_oldrid(i).ridges;
  addtab=addtab(:,2:size(addtab,2));
  Sample_complete_rid(i).ridges=[Sample_complete_rid(i).ridges addtab];
end
save([workdir 'tracing.newmeth.experiment.manual.completerid.mat'],'Sample_complete_rid');

%% checking plot
for samp_i_i=1:length(samples)
  samp_i=samples(samp_i_i);
  data=sampleData(samp_i);
  mat=data.Xcollapsed_1h1d;
  ppm=data.ppm_1h1d;
  peaks=struct();
  ridgenumbers=1:length(Sample_complete_rid(samp_i_i).ridges);
  for i=1:length(ridgenumbers)
    if ~isempty(Sample_complete_rid(samp_i_i).ridges(ridgenumbers(i)).result)
      % peaks(i).Ridges=ppm(Sample(samp_i_i).ridges(ridgenumbers(i)).result.colind);
      peaks(i).Ridges=Sample_complete_rid(samp_i_i).ridges(ridgenumbers(i)).result.ppm';
      peaks(i).RowInds=Sample_complete_rid(samp_i_i).ridges(ridgenumbers(i)).result.rowind';
      peaks(i).RidgeIntensities=Sample_complete_rid(samp_i_i).ridges(ridgenumbers(i)).result.intensity';
      peaks(i).CompoundNames=Sample_complete_rid(samp_i_i).ridges(ridgenumbers(i)).result.names;
      peaks(i).quantifiable=Sample_complete_rid(samp_i_i).ridges(ridgenumbers(i)).result.quanvec;
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
  saveas(fig,strcat(path,plotTitle,'.scatter.experiment.manual.comprid.fig'));
  close(fig);
end
