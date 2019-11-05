% this script is to match and output the ridges traced from simulation data
%% compound peak list will be mapped with traced ridges
%%% For each compound peak the closet ridge will be selected as the match ridge. The distance is normalized by ridge length
%% the recovery rates for peaks and compounds will also be calculated
%%% a ridge will be mapped to at most one peak in one compound and one peak in one compound will be matched with at most one ridges
%% the result will be transformed to RData
%% total MSE on ppm direction will also be calculated
%%% as difficult region is selected in experiment spectra, and these regions are often not annotated even by hand, this script does not apply for experiment data
close all;
clear all;
comp='/Users/yuewu/';
workdir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge_tracing_manuscript/analysis_res/result/peakmatching_simulated/'];
% addpath([comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/code/simulatespectra'])
cd(workdir);
load('../simulateddata/simulated.timeseries.ph.complex.more.mat');
datasources={'manual'}; % {'manual' 'automatic'};
load('../../data/spectral.real.list.withshift.mat');%% the true value
horzshift= 0.0;%-0.0025;
vertshift= 0.02;
scalmsethrehold=0.5;
compnameslist=struct();
compnameslist.full=fieldnames(strdatalist);
compnameslist.simplify={'ethanol' 'Valine' 'Acetate' 'Glycerol' 'Choline' 'alanine' 'Uridine' 'Leucine' 'Formate','DSS'};
list_match=struct();% the matched ridge information {ridges}
list_eval=struct();% the evaluation inforation: recovery rates {recovery_peak, recovery_molecule}
ssestr=struct();%store sse information 'sum of squared error'
%% this section will go into each sample by each ridge tracing method
%% it will map ridges to compound
%% it will also collect statistics information(sse,mse,and criteria for performance evaluation)
for datasource_i=1:length(datasources)
  datasource=datasources{datasource_i};
  load(['../simulated_quantification_' datasource '/' 'tracing.newmeth.simulated.' datasource '.mat']);
  matchstatpeak={};
  matchstatcompd={};
  list_eval.(datasource)={};
  ssestr.(datasource)={};
  matchstore={};%the matched ridge information for each data source
  for i=1:length(Sample)
    ssestr.(datasource){i}=0;
    rids=Sample(i).ridges;
    matchstatpeak{i}=0;
    matchstatcompd{i}=0;
    matchstore{i}=struct();
    npeak=0;%total number of peaks for the sample
    msevec=[];%store mse of each match
    tempstorecompd=struct();%store peak matches for each ridge
    %%this loop will find matches for each compound peaks
    if i==1 || i==2
      compnames=compnameslist.simplify;
    else
      compnames=compnameslist.full;
    end
    for compnamei=1:length(compnames) %for each compound
      compname=compnames{compnamei};
      comppeaklist=strdatalist.(compname);
      matchstore{i}.(compname)={};
      npeakcompd=0;%number of total peak for a compound
      %% match each compound peak in the peaklist with a ridge
      for t=1:size(comppeaklist,1)% for each peak of one compound
        npeakcompd=npeakcompd+1;
        if size(comppeaklist,2)==2
          ppmcandidate=comppeaklist(t,1);
        else
          ppmcandidate=comppeaklist;
        end
        peakcandid=ridgematch(ppmcandidate,rids);
        candidind=peakcandid.ind;
        peakcandid.compd=compname;
        peakcandid.ppm=ppmcandidate(1); %%if the peak shift, there will be multiple ppm, and the first one will be used as the name tag.
        if peakcandid.scalmse<scalmsethrehold
          ridid=['N' num2str(candidind)];
          if ~isfield(tempstorecompd,ridid)
            tempstorecompd.(ridid)={};
          end
          tempstorecompd.(ridid)=[tempstorecompd.(ridid) peakcandid];
        end
      end
      npeak=npeak+npeakcompd;
    end
    %% reduce match between peaklist and ridges to only the closest peak in peak list
    ridids=fieldnames(tempstorecompd);
    for ridi=1:length(ridids)
      ridid=ridids{ridi};
      tempstorerid=tempstorecompd.(ridid);
      lentempstorerid=length(tempstorerid);
      if lentempstorerid>1
        mseridvec=[];%%store mse for matches
        for tempstoreridi=1:lentempstorerid
          mseridvec=[mseridvec tempstorerid{tempstoreridi}.scalmse];
        end
        [temp minind]=min(mseridvec);
        tempstorerid=tempstorerid{minind};
      else
        tempstorerid=tempstorerid{1};
      end
      matchstatpeak{i}=matchstatpeak{i}+1;%%number of peak match for each sample
      % iind=length(usedrid_samecompd)+1;
      compname=tempstorerid.compd;
      tempstorerid.candidate.match=tempstorerid.ppm;
      tempstorerid.candidate.mse=tempstorerid.scalmse;
      matchstore{i}.(compname)=[matchstore{i}.(compname) tempstorerid.candidate];
      ssestr.(datasource){i}=ssestr.(datasource){i}+tempstorerid.sse;%%sum mse for each tracing
      msevec=[msevec tempstorerid.scalmse];
    end
    %%statistics for each compound and plot for each compound
    for compnamei=1:length(compnames)
      compname=compnames{compnamei};
      peakcell=matchstore{i}.(compname);
      if length(peakcell)>0
        matchstatcompd{i}=matchstatcompd{i}+1;%%number of molecule match for each sample
      else
        continue;
      end
      matchcomppeak=length(peakcell);%number of peak match for individual compound
      peaks=struct();
      ridgenumbers=1:length(peakcell);
      for j=ridgenumbers
        temppeak=peakcell{j}.result;
        if ~isempty(temppeak)
          peaks(j).Ridges=temppeak.ppm';
          peaks(j).RowInds=temppeak.rowind';
          peaks(j).RidgeIntensities=temppeak.intensity';
          peaks(j).CompoundNames=temppeak.names;
        else
          peaks(j).Ridges=[];
          peaks(j).RowInds=[];
          peaks(j).RidgeIntensities=[];
          peaks(j).CompoundNames=[];
        end
      end
      mathere=matrixstr{i};
      plotTitle=['check.' datasource '.' num2str(i) '.' compname 'nmatch.' num2str(matchcomppeak/npeakcompd) '.horz.' num2str(horzshift) '.'];%'.mse.' num2str(msestr.(datasource){i})
      fig=stackSpectra_paintRidges_3return(mathere,ppm,horzshift,vertshift,plotTitle,peaks,10);
      saveas(fig,strcat(workdir,plotTitle,'.scatter.fig'));
      close(fig);
    end
    % normalized by DSS {not used here}
    % dssintensity=matchstore{i}.('DSS'){1}.result.intensity;%% the location of DSS peak;
    % for compnamei=1:length(compnames)
    %   compname=compnames{compnamei};
    %   ridlist=matchstore{i}.(compname);
    %   for t=1:length(ridlist);
    %     matchstore{i}.(compname){t}.result.intensity=ridlist{t}.result.intensity./dssintensity;
    %   end
    % end
    tempstr=struct();
    tempstr.compd=matchstatcompd{i}/length(compnames);%%recall for compound
    tempstr.peak=matchstatpeak{i}/npeak;%%recall for peaks
    tempstr.npeak=npeak;%% positive condition
    tempstr.predposi=length(rids)-1;%%prediction positve
    tempstr.trueposi=matchstatpeak{i};%%true positive
    tempstr.sse=ssestr.(datasource){i};%%total sse
    tempstr.mserid=ssestr.(datasource){i}/tempstr.predposi; %%averaged mse for each mapped peak
    list_eval.(datasource){i}=tempstr;
    titlehis=['histogram.' datasource '.' num2str(i) '.mse'];
    fig=figure(), hold on
        histogram(msevec,30);
        title(titlehis);
    saveas(fig,strcat(workdir,titlehis,'.fig'));
    close(fig);
  end
  list_match.(datasource)=matchstore;
end
save([workdir 'eval.ridtracing.mat'],'list_eval','list_match');

%%addon performance evaluation calculation
%% only for automatic run of simulated data
% evaldata=list_eval.automatic;
% ml_evaltab=zeros([4 4]);
% for i=1:length(evaldata)
%   evaldataexp=evaldata{i};
%   ml_evaltab(i,1)=evaldataexp.trueposi/tempstr.predposi;%precision
%   ml_evaltab(i,2)=(tempstr.predposi-evaldataexp.trueposi)/tempstr.predposi;%FDR: False discover rates
%   ml_evaltab(i,3)=evaldataexp.trueposi/tempstr.npeak;%Recall
%   ml_evaltab(i,4)=(tempstr.npeak-evaldataexp.trueposi)/tempstr.npeak;%FNR: False negative rates
% end
% ml_eval=array2table(ml_evaltab,'VariableNames',{'precision','FDR','Recall','FNR'});
% datasource={'simualtion1','simualtion2','simualtion3','simualtion4'}';
% ml_eval=addvars(ml_eval,datasource,'Before','precision');
% save([workdir 'mleval.ridtracing.automatic.mat'],'ml_eval');
% writetable(ml_eval,'mleval.ridtracing.automatic.txt','Delimiter','\t');

%%check compound peaks
dirfig=workdir;
cd(dirfig);
load([workdir 'eval.ridtracing.mat'])
load('../simulateddata/simulated.timeseries.ph.complex.more.mat');
load('../../data/spectral.real.list.withshift.mat');%% the true value
i=18;
compnames=[fieldnames(strdatalist)];
compd=compnames{i};
% sampvec=[1 2 3 4];%[3 4];%[1 1 2 2 3 3 4 4];%%[3 3 4 4];%%
sampvec=[3 4];
% datasourcevec={'manual' 'manual' 'manual' 'manual'}; %%{'manual' 'manual'} ;%%{'manual' 'automatic' 'manual' 'automatic'};
datasourcevec={'manual' 'manual'};
posivec={'northwest' 'southwest'}; % {'northwest' 'southwest' 'north' 'south'}; % {'northwest' 'southwest' 'north' 'south' 'northeast' 'southeast' 'northeast' 'southeast'};%%{'northwest' 'southwest' 'north' 'south'};
% posivec={'northwest' 'southwest' 'north' 'south'};
for j=1:2%%1:8 %%
  pattern=['check.' datasourcevec{j} '.' num2str(sampvec(j)) '.' compd 'nmatch.*.horz.0..scatter.fig'];
  % pattern=['check.' datasourcevec{j} '.' num2str(sampvec(j)) '.' compd 'nmatch.*.mse.*.scatter.fig'];
  file=ls(pattern);
  file=regexprep(file,'[\n\r]+','')
  fig=openfig(file);
  movegui(fig,posivec{j});
end
close all;

%%evaluation of spectral complexity
%% attention, ridge need to be traced to feed into this
%%%%experiment spectra
clear all;
load(['../../data/sampleData.mat']);
% load(['../experiment_quantification_manul/archive/tracing.newmeth.experiment.manual.mat']);
load(['../experiment_quantification_manul_compare/tracing.newmeth.experiment.manual.mat'])
temptab1=zeros([2 5]);
sample=[1 4];%the index of the chosen sample in ridges data
sampleid=[1 7];%the index of the chosen sample in spectra data
% regionsele=[1.1 1.135; 1.3 1.35; 1.85 1.92; 2.00 2.11; 2.33 2.44; 3.1 3.15; 3.15 3.3; 3.68 3.8; 4.08 4.155; 4.28 4.35; 4.48 4.54; 5.36 5.41; 6.05 6.08; 6.51 6.55];
regionsele=[2.2 2.22; 0.97 0.99; 1.022 1.042; 0.992 1.01; 1.45 1.49; 1.31 1.34; 1.61 1.67; 5.785 5.81; 5.45 5.5; 8.25 8.36; 6.87 6.91; 7.405 7.43; 7.85 7.875; 3.15 3.2; 5.17 5.2; 7.95 8.05; 2.55 2.8; 2.4 2.55; 2.32 2.42; 6.5 6.6; 8.42 8.45; 5.2 5.26; 4.6 4.67; 3.87 3.92; 3.62 3.68; 1.15 1.2];
for samp_i_i=1:length(sample)
  samp_i=sample(samp_i_i);
  samp_isp=sampleid(samp_i_i);
  data=sampleData(samp_isp);
  mat=data.Xcollapsed_1h1d;
  ppm=data.ppm_1h1d;
  time=data.timesCollapsed_1h1d;
  % ridges=Sample(samp_i_i);
  ridges=Sample(samp_i);
  indrag=matchPPMs([-0.1,0.1],ppm);
  dssvec=max(mat(:,indrag(1):indrag(2)),[],2);
  objstr=complex_obj_func(mat,ppm,time',ridges,[-0.4 -0.2],10,[-0.4 10],regionsele,0.05,dssvec);
  temptab1(samp_i_i,1)=objstr.SNR;
  temptab1(samp_i_i,2)=objstr.Objcomppm;
  temptab1(samp_i_i,3)=objstr.Objcomdyn;
  temptab1(samp_i_i,4)=objstr.Objcomrange;
  temptab1(samp_i_i,5)=objstr.peakdensity;
end
%%%%simulated spectra
load(['../simulateddata/simulated.timeseries.ph.complex.more.mat']);
load(['../simulated_quantification_manual/tracing.newmeth.simulated.manual.mat']);
temptab2=zeros([4 5]);
regionsele=[-0.03 0.03; 0.6 0.65; 0.8 1.1; 1.15 1.185; 1.29 1.37; 1.4 1.6; 1.64 1.8; 1.9 2.1; ...
            2.22 2.31; 2.35 2.4; 2.88 2.94; 3.18 3.22; 3.3 3.5; 3.5 3.7; 3.7 3.85; ...
            3.87 4.0; 4.1 4.15; 4.18 4.28; 4.3 4.4; 4.75 4.95; 5.8 6; 7.8 7.9; 8.3 8.6; 8.8 9.1];
for samp_i=1:4
  ridges=Sample(samp_i);
  mat=matrixstr{samp_i};
  indrag=matchPPMs([-0.1,0.1],ppm);
  dssvec=max(mat(:,indrag(1):indrag(2)),[],2);
  objstr=complex_obj_func(mat,ppm,times,ridges,[-0.4 -0.2],10,[-0.4 10],regionsele,0.05,dssvec);
  temptab2(samp_i,1)=objstr.SNR;
  temptab2(samp_i,2)=objstr.Objcomppm;
  temptab2(samp_i,3)=objstr.Objcomdyn;
  temptab2(samp_i,4)=objstr.Objcomrange;
  temptab2(samp_i,5)=objstr.peakdensity;
end
comptab=array2table([temptab1; temptab2],'VariableNames',{'SNR','Objcomppm','Objcomdyn','Objcomrange','peakdensity'});
datasource={'experiment4','experiemnt10','simualtion1','simualtion2','simualtion3','simualtion4'}';
comptab=addvars(comptab,datasource,'Before','SNR');
save('spectra_complexity_tab.mat','comptab');
writetable(comptab,'spectra_complexity_tab.txt','Delimiter','\t');


%%test plot for valine region [2.22 2.3]
ppmregion=[2.22 2.3];
datasource='manual';
comp='/Users/yuewu/';
workdir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge_tracing_manuscript/analysis_res/result/peakmatching_simulated/'];
% addpath([comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/code/simulatespectra'])
cd(workdir);
load(['../simulated_quantification_' datasource '/' 'tracing.newmeth.simulated.' datasource '.mat']);
load([workdir 'eval.ridtracing.mat']);
load('../simulateddata/simulated.timeseries.ph.complex.more.mat');
i=4;
compname='Valine';
peakcell=list_match.(datasource){i}.(compname);
matchcomppeak=length(peakcell);%number of peak match for individual compound
peaks=struct();
ridgenumbers=1:length(peakcell);
mat=matrixstr{i};
jm=1;
for j=ridgenumbers
  temppeak=peakcell{j}.result;
  ppmvec=temppeak.ppm';
  if min(ppmvec)<ppmregion(1)|max(ppmvec)>ppmregion(2)
    continue;
  end
  if ~isempty(temppeak)
    peaks(jm).Ridges=temppeak.ppm';
    peaks(jm).RowInds=temppeak.rowind';
    peaks(jm).RidgeIntensities=temppeak.intensity';
    peaks(jm).CompoundNames=temppeak.names;
  else
    peaks(jm).Ridges=[];
    peaks(jm).RowInds=[];
    peaks(jm).RidgeIntensities=[];
    peaks(jm).CompoundNames=[];
  end
  jm=jm+1;
end
ppmrang=matchPPMs(ppmregion,ppm);
ppmind=ppmrang(1):ppmrang(2);
mathere=mat(:,ppmind);
ppmrere=ppm(ppmind);
fig=stackSpectra_paintRidges_3return(mathere,ppmrere,-0.0002,0.1,'valine2.22.2.3.pdf',peaks,10);
saveas(fig,strcat(workdir,'valine2.22.2.3.fig'));
print(fig,['valine2.22.2.3.pdf'],'-dpdf','-fillpage','-r2400');%'-painters','-fillpage',
close(fig);

%%% the table of ppm difference for simulated data (RMSD)
%%%% manual and automatic
load('../../data/spectral.real.list.withshift.mat')
load([workdir 'eval.ridtracing.mat'])
list_restab=struct();
sources={'manual'}% 'automatic'};
compds=fieldnames(strdatalist);
rowlen=0;
compdcellvec={};
ppmcellvec={};
for compd_i=1:length(compds)
  compd=compds{compd_i};
  ppmtab=strdatalist.(compd);
  ppms=ppmtab(:,1);
  rowlen=rowlen+length(ppms);
  compdsrep=cell(1,length(ppms));
  compdsrep(:)={compd};
  compdcellvec=[compdcellvec compdsrep];
  ppmscell=cell(1,length(ppms));
  for i=1:length(ppms)
    ppmscell{i}=num2str(ppms(i));
  end
  ppmcellvec=[ppmcellvec ppmscell];
end
compd_rids=strcat(compdcellvec,'_',ppmcellvec);
compduniq=unique(compdcellvec,'stable');
for source=sources
  source=source{1};
  data_source=list_match.(source);
  temptab=repmat(-1,rowlen,4);%%-1 indicate no matches there
  for sampi=1:length(data_source)
    data_samp=data_source{1,sampi};
    for compdele_i=1:length(compduniq)
      compdele=compduniq{compdele_i};
      localfieldnames=fieldnames(data_samp);
      if ~ismember(compdele,localfieldnames)
        continue;
      end
      data_compd=data_samp.(compdele);
      ppmtab=strdatalist.(compdele);
      for rid_i=1:length(data_compd)
        data_rid=data_compd{1,rid_i};
        nameind=find(strcmp(compd_rids,[compdele,'_',num2str(data_rid.match)]));
        temptab(nameind,sampi)=sqrt(data_rid.mse);
      end
    end
  end
  temptab2=array2table(temptab,'VariableNames',{'sample1','sample2','sample3','sample4'});
  temptab2=addvars(temptab2,compd_rids','Before','sample1');
  writetable(temptab2,[source '_peak_match_differences.txt'],'Delimiter','\t');
  list_restab.(source)=temptab2;
end
save('peak_match_differences','list_restab');
