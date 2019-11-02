%%this function will simulate spectra with different snr level and test the performance of ridge tracing
%%the peak will be a triplet moving in chemical shift
close all;
clear all;
comp='/Users/yuewu/';%the computer user location
workdir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge_tracing_manuscript/analysis_res/result/snr_check/'];%%the working folder
cd(workdir);
noiselevelabsolute=[5 50 100 500];
lambda=-0.00035;
totalnoiseregion=[-0.5 -0.2];
sampsize=50;
peakref_unknown=2.39;
concrange=[0 2000];% the whole range
times=1:sampsize;
%%simulation
load('../../data/spectral.real.list.mat');
load('../../data/spectral.complex.more.mat');
% matrix=zeros(sampsize,length(ppm));
peaksetting=struct('timevec',[1 sampsize],'ppmvec',[0.7 0.85],'polynlevel',1,'direction',0);
[matrix ppmrealvec]=poly_mov_spec_syn(strdata.unkownref,ppm,peakref_unknown,sampsize,peaksetting.timevec,peaksetting.ppmvec,peaksetting.polynlevel,peaksetting.direction,concrange(2));
matrixstrlist={};
matrixstrnoiselist={};
matsize=size(matrix);
noiserang=matchPPMs(totalnoiseregion,ppm);
noiseppmind=noiserang(1):noiserang(2);
rng(1);
for i=1:length(noiselevelabsolute)
  sigma=noiselevelabsolute(i);
  matrixnoise=matrix+random('Normal',0,sigma,matsize(1),matsize(2));
  matrixnoiseexp=zeros(matsize(1),matsize(2));
  for matrowi=1:matsize(1)
    matvechere=matrixnoise(matrowi,:);
    Y=hilbert(matvechere);
    X1=ifft(Y);
    timesresolu=1:1:length(matvechere);% the first point should be fixed to 1. normalization
    X1=X1.*exp(lambda*(length(matvechere)-timesresolu));
    y1=fft(X1);
    new_matvechere=real(y1);
    %%shift the intercept of spectral
    new_matvechere=new_matvechere-mean(new_matvechere(noiseppmind));
    matrixnoiseexp(matrowi,:)=new_matvechere;
  end
  matrixstrnoiselist{i}=matrixnoise;
  matrixstrlist{i}=matrixnoiseexp;
end
save('noisespec.mat','matrixstrlist','matrixstrnoiselist')
%%ridge tracing
load('noisespec.mat')
sample=1:length(noiselevelabsolute);
%% visual parameter
horzshift= -0.0025;
vertshift= 1E-3;
%% visual check the sample spectra to select region to work on
samp_i=2;
stackSpectra(matrixstrlist{samp_i},ppm,0.0,10,['noiselevel' num2str(noiselevelabsolute(samp_i))])
%% selected region to work on
regionsele=[0.68 0.88];
Sample=[];
Sample(1).ridges(1).parameters=[];
Sample(1).ridges(1).result=[];
for i = 2:length(sample)
    Sample(i).ridges(1)=Sample(1).ridges(1);
end
%%%%%%%% test run %%%%%%%%%%%%%%
%%%%%the two main tuning parameter
thredseg=40; %% the more curvy the ridges the higher this value
maxaddon=1; %% if find out ridges not traced, increase this number to the number of ridges in the region
%%%%
samp_i=4;
regionhere=regionsele;
mat=matrixstrlist{samp_i};
[returndata]=ridgetrace_power2_ext(mat,ppm,times',regionhere,workdir,thredseg,maxaddon);
for samp_i=sample
  mat=matrixstrlist{samp_i};
  disp(['sample ' num2str(samp_i)]);
  plotTitle=[num2str(regionsele(1)),'.',num2str(regionsele(2)),'.',num2str(samp_i),'.testplot'];
  [returndata]=ridgetrace_power2_ext(mat,ppm,times',regionhere,workdir,thredseg,maxaddon);
  fig=gcf;
  saveas(fig,strcat(workdir,plotTitle,'.surf.automatic.fig'));
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
  fig=stackSpectra_paintRidges_3return(mathere,ppmhere,-0.0025,0.01,plotTitle,peakshere,10);
  saveas(fig,strcat(workdir,plotTitle,'.scatter.automatic.fig'));
  close(fig);
end
save('noiseridge_trace.mat','Sample');
%%plotting
for samp_i=sample
  mat=matrixstrlist{samp_i};
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
  saveas(fig,strcat(workdir,plotTitle,'.scatter.automat.fig'));
  close(fig);
end
%%smaller region plot
smallregionppm=[0.4 1.0];
for samp_i=sample
  mat=matrixstrlist{samp_i};
  reg=matchPPMs(smallregionppm,ppm);
  ind=reg(1):reg(2);
  mathere=mat(:,ind);
  ppmhere=ppm(ind);
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
  fig=stackSpectra_paintRidges_3return(mathere,ppmhere,-0.00025,0.02,plotTitle,peaks,10);
  saveas(fig,strcat(workdir,plotTitle,'.scatter.automat.smallregion.fig'));
  % fig.PaperOrientation='landscape';
  % fig.Position=[0 500 1600 400];
  % pos=get(fig,'Position');
  set(fig,'PaperPositionMode','Auto','PaperUnits','inches','PaperSize',[7.777777777777778,5.833333333333332]);%%%get this value from figure object
  % set(fig,'PaperPositionMode','auto');
  print(fig,[plotTitle,'.scatter.automat.smallregion.pdf'],'-dpdf','-r2400');%'-painters','-fillpage',
  saveas(fig,['simulated.spectra.complexity.' num2str(4) '.selected.fig']);
  close(fig);
end
%%practical SNR
load("noiseridge_trace.mat")
snrlist=[];
for samp_i=sample
  mathere=matrixstrlist{samp_i};
  noisereg=matchPPMs(noiserang,ppm);
  noisemat=mathere(:,noisereg(1):noisereg(2));
  %% signal to noise ratio
  % SNR=max(mathere(:))/std(noisemat(:));
  SNR=mean(max(mathere,[],2)./std(noisemat,0,2));
  snrlist=[snrlist SNR];
end
save('noisemat_snr.mat','snrlist');
