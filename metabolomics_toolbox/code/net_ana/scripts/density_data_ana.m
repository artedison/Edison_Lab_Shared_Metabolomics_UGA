%%this is part of the new dataset
% WT high vs low density
set(0,'DefaultFigureVisible','on');
close all;
clear all;
comp='/Users/yuewu/';%the computer user location
pardir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge.net/result_addition/fpca_cn_wt_density/'];
% OR
% pardir=[pwd() '/'];
workdir=[pardir];%%the working folder
datadir=[pardir];%the data folder
cd(workdir);
load([datadir 'wt_density_2exp.mat']);%the ridge tracking result as input
% densvec: density of each experiment
% expnames: name of each experiment
% matlist: intensity of tracked ridges of each experiment (time*n_ridges)
% ppmlist: ppm of tracked ridges of each experiment (time*n_ridges)
% namesall: name of matched ridges between the two experiment
% ppmmatch_ind_all: ridge match index (n_experiment*n_matched_ridges)
nsample=length(expnames);
rng(1);
sampseq_all=1:nsample;
plot_ylim=[-2.2,2.2];
replicates=[];
names={};
ppm_all={};
ymat_raw=[];
ind_shift_vec=[0];
% whole data matrix
for isample=sampseq_all
  ymat_raw=[ymat_raw matlist{isample}];
  locname=namelist{isample};
  ch_ind=ppmmatch_ind_all(isample,:);
  no_nanind= ~isnan(ch_ind);
  locname(ch_ind(no_nanind))=namesall(no_nanind);
  locname_front={};
  for namesin=locname
    nameparts=split(namesin,' ');
    locname_front=[locname_front nameparts(1)];
  end
  names=[names locname_front];
  %
  ppmvec=num2cell(mean(ppmlist{isample},1));
  ppmadd=cellfun(@(x) num2str(x),ppmvec,'UniformOutput',false);
  ppm_all=[ppm_all ppmadd];
  replicates=[replicates repmat(isample,[1,length(namelist{isample})])];
  ind_shift_vec=[ind_shift_vec size(ymat_raw,2)];
end
%
ymat=ymat_raw;
nridges=size(ymat,2);
ntime=30;
nDer=1;%maximal derivative to consider
xvec=1:ntime;
%% center and scale for each feature
for j=1:nridges
  shiftv=ymat(:,j)-mean(ymat(:,j));
  ymat(:,j)=shiftv./std(shiftv);
end
% select matched part
accumu_ind=[];
for irid=1:size(ppmmatch_ind_all,2)
  ch_ind=ppmmatch_ind_all(:,irid);
  no_nanind= ~isnan(ch_ind);
  real_ind=ind_shift_vec(no_nanind)+ch_ind(no_nanind)';
  accumu_ind=[accumu_ind real_ind];
end
ymat_match=ymat(:,accumu_ind);
replicates_match=replicates(accumu_ind);
names_match=names(accumu_ind);
ppm_all_match=ppm_all(accumu_ind);
% Funcitonal/smoothing representation for each feature
loglambda_vec= -4:0.25:4;
res=smooth_derivative(ymat_match,xvec,loglambda_vec,nDer);
% Funcitonal/smoothing representation for each PC
npc=5;%maximum number of pc to calculat
lambdapc=0.0;%lambda for pc
fdres=res.spfd_selec;
nsmooth=nDer+2;
norder=nsmooth+2;
nbasis=length(xvec)+norder-2;
bbasis=create_bspline_basis([min(xvec) max(xvec)],nbasis,norder,xvec);
fdParpc=fdPar(bbasis,nsmooth,lambdapc);
res=fda_pca(fdres,fdParpc,npc);
close all;
%
harmscr=res.fdapcastr.harmscr;
names_unique=unique(names_match);
ncompd=length(names_unique);
%
annotnames_ind=find(~cellfun(@(x) any(regexp(x,'[uU]nknown')),names_match,'UniformOutput',true));
markerlist={'o' '^'};%for each conditions high vs low density
names_unique_pres={'QA', 'ethanol', 'adenosine', 'uracil', 'alanine'};%instead of unique(names) just select a few to show
unsele_ind=cellfun(@(x) ~ismember(x,names_unique_pres),names_match,'UniformOutput',true);
names_pres=names_match;
names_pres(unsele_ind)={'unknown'};
names_unique_all=['unknown' names_unique_pres];
ncompd=length(names_unique_all);
%
colorsmap=hsv(ncompd);
colorsmap=colorsmap(randperm(size(colorsmap,1)),:);
color_noanno_ind=cellfun(@(x) any(regexp(x,'[uU]nknown')),names_unique_all,'UniformOutput',true);
colorsmap(color_noanno_ind,:)=repmat([0.5 0.5 0.5],[1,1]);
legends={};
h2=figure();
hold on;
for namei=1:length(names_unique_all)
  name=names_unique_all{namei};
  for replicate=unique(replicates_match)
    ind=find(replicates_match==replicate&strcmp(names_pres,name));
    if strcmp(name,'unknown')
      scatter(harmscr(ind,1),harmscr(ind,2),[],colorsmap(namei,:),markerlist{replicate});%
    else
      scatter(harmscr(ind,1),harmscr(ind,2),200,colorsmap(namei,:),markerlist{replicate},'filled');
    end
    legends=[legends ['density ' num2str(densvec(replicate)) ' ' name]];
  end
end
xlabel('pc1');
ylabel('pc2');
title(['pca epl']);
legend(legends);
saveas(h2,[workdir,'fdapca_score_plot.density.fig']);
close all;

% loading plot and variance plot
plot_pca_fd(res.fdapcastr,1,[1 2]);
figh=gcf;
allaxes=findall(figh,'type','axes')
set(allaxes(1),'ylim',plot_ylim);
set(allaxes(2),'ylim',plot_ylim);
saveas(figh,[workdir,'fdapca_loading_plot.fig']);
close all;
fig_prop=figure();
plot(1:npc,res.fdapcastr.varprop(1:npc),'-o');
xlabel('Eigenvalue Number');
ylabel('Eigenvalue');
saveas(fig_prop,[workdir,'fdapca_prop_plot.fig']);
close all;

% plot for selected compounds
markerlist={'o' '^'};%for each conditions high vs low density
colorline={'r','g'};
for compd=names_unique_pres
  legendsrec={};
  h2=figure();
  hold on;
  for replicate=unique(replicates_match)
    expname=expnames{replicate};
    ind=find(replicates_match==replicate&strcmp(names_match,compd));
    yplot=ymat_match(:,ind);
    yplot=yplot(:);
    xplot=repmat(timelist{replicate}(:,1),[length(ind),1]);
    scatter(xplot,yplot,200,colorline{replicate},markerlist{replicate});
    for linei=1:length(ind)
      lineind=(((linei-1)*ntime)+1):((linei*ntime));
      line(xplot(lineind),yplot(lineind),'LineWidth',2,'LineStyle','--','Color',colorline{replicate});
    end
    legendsrec=[legendsrec {[expname]}];
  end
  xlabel('time');
  ylabel('quantification');
  title([compd]);
  % legend(legendsrec);
  saveas(h2,[workdir,'time_series_plot.',compd{1},'_in_allsample.fig']);
  close all;
end
