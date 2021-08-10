%% The script fits smoothing curve for all features from ridges tracking by functional data analysis (FDA)
%% It calculate smoothed intensity and derivative for all peaks

set(0,'DefaultFigureVisible','on');
close all;
clear all;
%% the user will need to modify the path here for local run
comp='/Users/yuewu/';%the computer user location
pardir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge.net/result_reprod/'];
% OR
% pardir=[pwd() '/'];
workdir=[pardir 'result/smoothing/'];%%the working folder
datadir=[pardir 'data/'];%the data folder
continue_dir=[pardir 'result_data/'];
cd(workdir);
load([continue_dir 'tracing.newmeth.experiment.manual.completerid.mat']);%the ridge tracking result as input
load([datadir 'unshared/sampleData.mat']);%please refer to the orginal CIVM-NMR paper for the data here

ntime=52;%number of time points
%% construct data matrix for all replicates
samp_ind=[];%%the start index for each replicate on ymat_all
samp_loc=[];% the sample id and ridge id for each ridge
ymat_all=[];%%the whole data matrix
ppmmat_all=[];
timemat_all=[];
listloc={};%sample_i ppm
xvec=1:ntime;
matlist={};%intensity matrix for each replicate
ppmlist={};%ppm matrix for each replicate
timelist={};%time matrix for each replicate
namelist={};%name vector for each replicate
ridgerangelist={};%time range for each ridge in each replicate
samp_loc_list={};%the sample id and ridge id for replicate
%% reformat ridge data
for samp_i=1:length(Sample_complete_rid)
  sample_rep=Sample_complete_rid(samp_i);
  nridges=size(sample_rep.ridges,2)-1;
  ymat=NaN(ntime,nridges);
  ppmmat=NaN(ntime,nridges);
  timemat=NaN(ntime,nridges);
  samp_loc_here=NaN(nridges,2);
  namelist{samp_i}={};
  ridrange=NaN(nridges,2);
  for i=1:nridges
    locstr=sample_rep.ridges(i+1).result;
    timevec=locstr.rowind;
    intensity=locstr.intensity;
    ppmloc=locstr.ppm;
    % timeseq=locstr.time;
    %%reordering ridge information
    [timevec sortind]=sort(timevec);
    intensity=intensity(sortind);
    ppmloc=ppmloc(sortind);
    % timeseq=timeseq(sortind);
    indvec=matchPPMs(timevec',xvec);
    %%intensity information
    locintvec=NaN(1,ntime);
    locintvec(indvec)=intensity;
    rangetime=[min(timevec) max(timevec)];
    locintvec(1:(rangetime(1)-1))=locintvec(rangetime(1));
    locintvec((rangetime(2)+1):ntime)=locintvec(rangetime(2));
    % locintvec(find(isnan(locintvec)))=min(intensity);
    ymat(:,i)=locintvec;
    %%ppm information
    locppmvec=NaN(1,ntime);
    locppmvec(indvec)=ppmloc;
    locppmvec(1:(rangetime(1)-1))=locppmvec(rangetime(1));
    locppmvec((rangetime(2)+1):ntime)=locppmvec(rangetime(2));
    ppmmat(:,i)=locppmvec;
    % time information
    timemat(:,i)=sampleData(samp_i).timesCollapsed_1h1d;
    %
    namesvec=unique(locstr.names);
    namelist{samp_i}=[namelist{samp_i} namesvec(1)];
    ridrange(i,:)=rangetime;
    samp_loc_here(i,:)=[samp_i i];
    listloc=[listloc [samp_i min(locppmvec) max(locppmvec)]];
  end
  samp_ind=[samp_ind size(ymat_all,2)+1];
  ymat_all=[ymat_all ymat];
  ppmmat_all=[ppmmat_all ppmmat];
  timemat_all=[timemat_all timemat];
  samp_loc=[samp_loc; samp_loc_here];
  matlist{samp_i}=ymat;
  ppmlist{samp_i}=ppmmat;
  timelist{samp_i}=timemat;
  ridgerangelist{samp_i}=ridrange;
  samp_loc_list{samp_i}=samp_loc_here;
end
%%preprocessing data matrix
% shiftv=min(ymat_all(:))-0.000001;
% ymat_all_trans=log(ymat_all-shiftv);
ymat_all_trans=ymat_all;
%%smoothing
loglambda_vec= -4:0.25:4;
listres={};
smoothlist={};
smooth_1d_list={};
ymat_smooth_all=[];
ymat_1d_all=[];
temp_smooth=[];
temp_1d=[];
for ridi=1:size(ymat_all_trans,2)
  ridi
  measy=ymat_all_trans(:,ridi);
  res=smooth_derivative(measy,xvec,loglambda_vec,1,false);
  esty=eval_fd(xvec,res.spfd_selec);
  esty_1d=eval_fd(xvec,res.spfd_selec,1);
  res.variance=var((measy-esty)/mean(esty));
  listres{ridi}=res;
  temp_smooth=[temp_smooth esty];
  temp_1d=[temp_1d esty_1d];
  if any(samp_ind==(ridi+1))
    prelen=length(smoothlist);
    smoothlist{prelen+1}=temp_smooth;
    smooth_1d_list{prelen+1}=temp_1d;
    ymat_smooth_all=[ymat_smooth_all temp_smooth];
    ymat_1d_all=[ymat_1d_all temp_1d];
    temp_smooth=[];
    temp_1d=[];
  end
end
prelen=length(smoothlist);
smoothlist{prelen+1}=temp_smooth;
smooth_1d_list{prelen+1}=temp_1d;
ymat_smooth_all=[ymat_smooth_all temp_smooth];
ymat_1d_all=[ymat_1d_all temp_1d];
save([workdir 'tracing.smoothed.mat'],'listres','samp_ind','ymat_all','ymat_smooth_all','ymat_1d_all','matlist','smoothlist','smooth_1d_list','ppmlist','namelist','ridgerangelist','samp_loc','listloc','xvec','ppmmat_all','samp_loc_list','timelist','timemat_all')%
% listres: all the fda output for each ridge
% samp_ind: the start index for each replicate on ymat_all
% ymat_all: the whole data matrix. ntime*(nreplicate*nridge)
% ymat_smooth_all: the whole data smoothed matrix.  ntime*(nreplicate*nridge)
% ymat_1d_all: the whole data 1st derivative matrix.  ntime*(nreplicate*nridge)
% matlist: the whole data matrix for each replicate. ntime*nridge
% smoothlist: the whole data smoothed matrix for each replicate. ntime*nridge
% smooth_1d_list: the whole data 1st derivative matrix for each replicate. ntime*nridge
% ppmlist: the whole data ppm matrix for each replicate. ntime*nridge
% namelist: the whole data name vector for each replicate.
% ridgerangelist: time range for each ridge in each replicate. nridge*2
% samp_loc: the sample id and ridge id for each ridge in each replicate. (nreplicate*nridge)*2
% listloc: the sample id and ppm range for each ridge in each replicate. (nreplicate*nridge)*3
% xvec: the time index vector
% ppmmat_all: the whole ppm matrix. ntime*(nreplicate*nridge)
% samp_loc_list: the sample id and ridge id for ridge. nridge*2
% timelist: the whole data time matrix for each replicate. ntime*nridge
% timemat_all: the whole time matrix. ntime*(nreplicate*nridge)
