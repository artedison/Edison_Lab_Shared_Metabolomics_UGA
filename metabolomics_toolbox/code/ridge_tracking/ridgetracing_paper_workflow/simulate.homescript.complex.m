%% this script use homemade script to simulate spectra with pH changes (ppm shift)
%% this script will simulate spectral with different complex level (compound number and dynamics in intensity) and noise level
%% ref: Modelling the acid/base 1H NMR chemical shift limits of metabolites in human urine
%% ref: Bayesian estimation of the number of protonation sites for urinary metabolites from NMR spectroscopic data
%% ref: The NMR Chemical Shift pH Measurement Revisited: Analysis of Error and Modeling of a pH Dependent Reference
close all;
clear all;
%%setting up paths for directory
comp='/Users/yuewu/';%the computer user location
workdir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge_tracing_manuscript/analysis_res/' 'result/simulateddata/'];%%the working folder
cd(workdir);
%%parameter setting
% rng(1);%%random number seed
complexlevel=[1 2];
phrag=[4.0 6.0];%8.0
sampsize=50;
concrange=[0 2000];% the whole range
noiselevelabsolute=[0.5 5];%noise level
pH_ref=7.4;
pHpara=struct();
pHpara.Acetate=[4.578 2.089 1.910]; %pKa, acid ppm limit, base ppm limit
pHpara.Formate=[3.555 8.234 8.448];
lambda=-0.00035;% the exponential window function for fid, need to be modified according to inspection of result spectra
totalnoiseregion=[-0.5 -0.2];%% the region of no signal to shift spectra
ref=['DSS'];
peakref_unknown=2.39;
compd_ref_unkown='unkownref';
intensityscal=10;% the division for highest intensity for a subset of compound for the complex high condition 2
times=1:sampsize;
seqpH=linspace(phrag(1),phrag(2),sampsize);
strdataannotname={'ethanol' 'Valine' 'Acetate' 'Glycerol' 'Choline' 'alanine' 'Uridine' 'Leucine' 'Formate'  'Butanol' 'Caffeine' 'Serine' 'Purine' 'unkownref' 'DSS'};%%when update this, it is a good idea to keep DSS at end to ensure DSS included in the spectral as reference
unknown_paralist=struct();%%parameters for added synthized spectra
unknown_paralist.unknown1=struct('timevec',[1 sampsize],'ppmvec',[1.95 1.95],'polynlevel',1,'direction',0);%%set points and polynomial power
unknown_paralist.unknown2=struct('timevec',[1 sampsize],'ppmvec',[4.2 4.26],'polynlevel',1,'direction',0);
unknown_paralist.unknown3=struct('timevec',[1 25 sampsize],'ppmvec',[4.85 4.9 4.85],'polynlevel',2,'direction',0);
unknown_paralist.unknown4=struct('timevec',[1 sampsize],'ppmvec',[3.35 3.45],'polynlevel',1,'direction',-1);

%%don't run this block if you don't have corresponding gissmo library in the correct location
%%the data from gissmo will be provided as .mat file
%%spectral from library
% gissmodir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/data/gissmo/'];
% dataloc='/simulation_1/B0s/sim_600MHz';%% the relative location in gissmo folder
% peaklistloc='/simulation_1/peaks/sim_600MHz_peaks_standard.csv';
% ppmlimit=[-0.5 10];%the range of ppm that will be inlcuded in the data
% strdataannotbmrb={'bmse000297' 'bmse000052' 'bmse000191' 'bmse000184' 'bmse000285' 'bmse000028' 'bmse000158' 'bmse000042' 'bmse000203' 'bmse000447' 'bmse000206'  'bmse000048' 'bmse000454' 'unknown' 'bmse000795' }; %%ids for the compound in order
% strdata=struct();
% strdatalist=struct();
% ppm=[];
% %%read data from gissmo folder
% for i=1:length(strdataannotname)
%   %%spectra
%   name=strdataannotname{i};
%   filelocat=[gissmodir strdataannotbmrb{i} dataloc];
%   filelocat=ls([filelocat '*']);%accomodate files with different names
%   filelocat=regexprep(filelocat,'[\n\r]+','');
%   tab=readtable(filelocat);
%   tab=table2array(tab);
%   ragind=find(tab(:,1)>ppmlimit(1)&tab(:,1)<ppmlimit(2));
%   tab=tab(ragind,:);
%   strdata.(name)=tab(:,2);
%   ppm=tab(:,1)';
%   %%peaklist
%   if strcmp(strdataannotbmrb{i},'unknown') %%unknown are not included here
%     continue;
%   end
%   filelocat=[gissmodir strdataannotbmrb{i} peaklistloc];
%   tab=readtable(filelocat);
%   tab=table2array(tab);
%   strdatalist.(name)=tab;%%PPM, amplititude
% end
% save('../../data/spectral.complex.more.mat','strdata','ppm');
% save('../../data/spectral.real.list.mat','strdatalist');

%%add moving peaks to true peak list
load('../../data/spectral.real.list.mat')
refppmshift=struct();
refppmshift.Acetate=strdatalist.Acetate(1,1);%%it contains peaks list if peak didn't move and peak location for moving peaks(with one peak for these compound) or unkown compound
refppmshift.Formate=strdatalist.Formate(1,1);
strdatalist.Acetate=[];
strdatalist.Formate=[];
%% simulate spetra
load('../../data/spectral.complex.more.mat');
matrixstrlist={};
matrixstrnoiselist={};
truevallist={};
for compi=complexlevel
  if compi==1
    strdataannotnamehere=[strdataannotname(1:9) ref];
  else
    strdataannotnamehere=strdataannotname;
  end
  matrix=zeros(sampsize,length(ppm));
  truevallisttemp=struct();
  for j=1:length(strdataannotnamehere)
    compname=strdataannotnamehere{j};
    truevallisttemp.(compname)=[];
  end
  %%named spectra
  for i=1:sampsize
    spec_sum=zeros(length(ppm),1);
    for j=1:length(strdataannotnamehere)
      if compi==2
        if rem(j,3)==1
          seqconc=linspace(concrange(1),concrange(2),sampsize);
        elseif rem(j,3)==2
          seqconc=linspace(concrange(1),concrange(2)/intensityscal,sampsize);
        else
          seqconc=linspace(concrange(1),concrange(2),sampsize);
        end
      else
        seqconc=linspace(concrange(1),concrange(2),sampsize);
      end
      compname=strdataannotnamehere{j};
      spectral=strdata.(compname);
      %%ppm shift based on pH
      if ismember(compname,fieldnames(pHpara))
        parahere=pHpara.(compname);
        [spectral ppmnewshift]=ppmshift(spectral,parahere(3),parahere(2),parahere(1),seqpH(i),ppm,pH_ref);
        if compi==1
          newppmval=ppmnewshift+refppmshift.(compname);
          strdatalist.(compname)=[strdatalist.(compname) newppmval];
        end
      end
      %% the reference peak
      if ismember(compname,ref)%%DDS is not changed in intensity through time
        spec_sum=spec_sum+spectral*seqconc(end);
        truevallisttemp.(compname)=[truevallisttemp.(compname) seqconc(end)];
        continue;
      end
      %% different trend for different compound
      if rem(j,2)==1
        spec_sum=spec_sum+spectral*seqconc(i);
        truevallisttemp.(compname)=[truevallisttemp.(compname) seqconc(i)];
      else
        spec_sum=spec_sum+spectral*seqconc(sampsize-i+1);
        truevallisttemp.(compname)=[truevallisttemp.(compname) seqconc(sampsize-i+1)];
      end
    end
    matrix(i,:)=spec_sum';
  end
  %%add unknown spectral
  %%unknown compound spectra will be added to add complexity of the spectra, they will be adapted from succinate spectra
  if compi==2
    unknown_mat=zeros(sampsize,length(ppm));
    unknownnames=fieldnames(unknown_paralist);
    for unknownnamei=1:length(unknownnames)
      unknownname=unknownnames{unknownnamei};
      para=unknown_paralist.(unknownname);
      [mataddon ppmrealvec]=poly_mov_spec_syn(strdata.(compd_ref_unkown),ppm,peakref_unknown,sampsize,para.timevec,para.ppmvec,para.polynlevel,para.direction,concrange(2)/10);
      unknown_mat=unknown_mat+mataddon;
      strdatalist.(unknownname)=ppmrealvec;
    end
    matrix=matrix+unknown_mat;
  end
  %%noise add and line broadening
  matrixstrtemp={};
  matrixstrnoisetemp={};
  matsize=size(matrix);
  noiserang=matchPPMs(totalnoiseregion,ppm);
  noiseppmind=noiserang(1):noiserang(2);
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
    matrixstrnoisetemp{i}=matrixnoise;
    matrixstrtemp{i}=matrixnoiseexp;
  end
  matrixstrlist{compi}=matrixstrtemp;
  matrixstrnoiselist{compi}=matrixstrnoisetemp;
  truevallist{compi}=truevallisttemp;
end
%%reformat the data
%%% the order: simple_low_noise simple_high_noise complex_low_noise complex_high_noise
matrixstr={matrixstrlist{1}{1} matrixstrlist{1}{2} matrixstrlist{2}{1} matrixstrlist{2}{2}};
matrixstrnoise={matrixstrnoiselist{1}{1} matrixstrnoiselist{1}{2} matrixstrnoiselist{2}{1} matrixstrnoiselist{2}{2}};
save('simulated.timeseries.ph.complex.more.mat','matrixstr','matrixstrnoise','ppm','times','truevallist')
save('../../data/spectral.real.list.withshift.mat','strdatalist');


%%add unknown peaks into the true value
load('simulated.timeseries.ph.complex.more.mat')
changedpart=truevallist{2};
unknownnames=fieldnames(unknown_paralist);
for unknownnamei=1:length(unknownnames)
  unknownname=unknownnames{unknownnamei};
  para=unknown_paralist.(unknownname);
  refpeak=max(strdata.unkownref);
  direction=para.direction;
  conchigh=concrange(2)/10;
  if direction==1
    ratio=linspace(0,conchigh,sampsize);
  elseif direction== -1
    ratio=linspace(conchigh,0,sampsize);
  else%%direction==0
    ratio=repmat(conchigh,[1 sampsize]);
  end
  changedpart.(unknownname)=refpeak*ratio;
end
truevallist{2}=changedpart;
save('simulated.timeseries.ph.complex.more.addunknown.mat','matrixstr','matrixstrnoise','ppm','times','truevallist')
%% visualize check
%% plot a general spectra
plotr(ppm,matrixstr{4}(20,:));
plotr(ppm,matrixstrnoise{3}(10,:));
%%test plot of ppm shift region surface plot
matrix=matrixstr{4};
testregion=[1.8 2.2]%[1.8 2.2];%%[8.1 8.5];%%[1.8 2.2]
ppmrang=matchPPMs(testregion,ppm);
ppmind=ppmrang(1):ppmrang(2);
mattest=matrix(:,ppmind);
ppmtest=ppm(ppmind);
fig=figure(), hold on
    surf(ppmtest,times',mattest,'FaceColor','Interp');
    % surf(mattest);
    shading interp;
    set(gca,'xdir','reverse');
    ylabel('time');
    zlabel('intensity');
    title(['shiftregion']);
    xlabel('ppm');
%%test plot for crowded region
plotr(ppm,matrixstr{1}(20,:));
%%DSS region
testregion=[-0.5 0.5];
ppmrang=matchPPMs(testregion,ppm);
ppmind=ppmrang(1):ppmrang(2);
mattest=matrix(:,ppmind);
ppmtest=ppm(ppmind);
fig=figure(), hold on
    surf(ppmtest,times',mattest,'FaceColor','Interp');
    shading interp;
    set(gca,'xdir','reverse');
    ylabel('time');
    zlabel('intensity');
    title(['shiftregion']);
    xlabel('ppm');
%%stack plot of whole spectra
samp_i=4;
mat=matrixstr{samp_i};
horzshift=0.0;
vertshift=10;
stackSpectra(mat,ppm,horzshift,vertshift,'simulated N. crassa NMR metabolic fingerprint over time')


%%compare the experiment spectral and simulated spectral
%%%%experiment
load([comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge_tracing_manuscript/analysis_res/' 'data/sampleData.mat']);
fig1=figure();
plotr(sampleData(1).ppm_1h1d,sampleData(1).Xcollapsed_1h1d(25,:));
%%%%simulated
fig2=figure();
plotr(ppm,matrixstr{2}(10,:));

%%save figures
cd(workdir);
load('simulated.timeseries.ph.complex.more.mat');
load('../../data/spectral.real.list.withshift.mat');
%% visual parameter
for i=1:length(matrixstr)
  plotTitle=['simulated.spectra.complexity.' num2str(i)];
  stackSpectra(matrixstr{i},ppm,0.0,10,['N. crassa NMR metabolic fingerprint over time complexity ' num2str(i)]);
  fig=gcf;
  saveas(fig,strcat(workdir,plotTitle,'.fig'));
  close(fig);
end
for i=1:length(matrixstr)
  file=['simulated.spectra.complexity.' num2str(i) '.fig'];
  fig=openfig(file);
  print(fig,['simulated.spectra.complexity.' num2str(i) '.pdf'],'-fillpage','-dpdf','-r2400');%'-painters',
  close(fig);
end
%%specific region for the four simulated samples
ppmrag=[2.8 4.4];
reg=matchPPMs(ppmrag,ppm);
ind=reg(1):reg(2);
samprag=1:50;
for i=1:length(matrixstr)
  plotTitle=['simulated.spectra.complexity.' num2str(i)];
  mat=matrixstr{i};
  mathere=mat(samprag,ind);
  ppmhere=ppm(ind);
  stackSpectra(mathere,ppmhere,-0.001,50,['N. crassa NMR metabolic fingerprint over time complexity ' num2str(i)]);
  fig=gcf;
  saveas(fig,strcat(workdir,plotTitle,'_part2.fig'));
  close(fig);
end
for i=1:length(matrixstr)
  file=['simulated.spectra.complexity.' num2str(i) '_part2.fig'];
  fig=openfig(file);
  print(fig,['simulated.spectra.complexity.' num2str(i) '_part2.pdf'],'-fillpage','-dpdf','-r2400');%'-painters',
  close(fig);
end

%%plot specific region of complexity 4 figure
ppmregion=[0.5 3.0];
selerang=matchPPMs(ppmregion,ppm);
ppmind=selerang(1):selerang(2);
matselec=matrixstr{4}(:,ppmind);
ppmselec=ppm(ppmind);
stackSpectra(matselec,ppmselec,-0.0025,10,['N. crassa NMR metabolic fingerprint over time complexity ' num2str(4)]);
%%%%modified figure size to longer
fig=gcf;
% fig.PaperOrientation='landscape';
fig.Position=[0 500 1600 400];
pos=get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','inches','PaperSize',[22.22222222222222,5.555555555555555]);%%%get this value from figure object
% set(fig,'PaperPositionMode','auto');
print(fig,['simulated.spectra.complexity.' num2str(4) '.selected.pdf'],'-dpdf','-r2400');%'-painters','-fillpage',
saveas(fig,['simulated.spectra.complexity.' num2str(4) '.selected.fig']);
close(fig);
