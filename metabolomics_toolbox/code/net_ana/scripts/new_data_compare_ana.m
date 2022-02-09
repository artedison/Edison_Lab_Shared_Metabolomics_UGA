%%this is part of the new dataset
% WT glucose vs QA as carbon sources
set(0,'DefaultFigureVisible','on');
close all;
clear all;
comp='/Users/yuewu/';%the computer user location
pardir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge.net/result_addition/fpca_cn_wt_qa/'];
% OR
% pardir=[pwd() '/'];
workdir=[pardir];%%the working folder
datadir=[pardir];%the data folder
continue_dir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge.net/result_reprod/result_data/'];
cd(workdir);
load([datadir 'wt_comp_exp.mat']);%the ridge tracking result as input
load([continue_dir 'tracing.smoothed.mat']);
%
aerobic_ind=[1 2 3];
matlist_comb=[matlist(aerobic_ind) matlist_add(1)];
namelist_comb=[namelist(aerobic_ind) namelist_add(1)];
ppmlist_comb=[ppmlist(aerobic_ind) ppmlist_add(1)];
timelist_comb=[timelist(aerobic_ind) timelist_add(1)];
expnames={'glucose 1','glucose 2','glucose 3','QA'};
%
nsample=length(expnames);
rng(1);
sampseq_all=1:nsample;
%
plot_ylim=[-2.2,2.2];
replicates=[];
names={};
ppm_all={};
ymat_raw=[];
ind_shift_vec=[0];
% whole data matrix
for isample=sampseq_all
  ymat_raw=[ymat_raw matlist_comb{isample}];
  names=[names namelist_comb{isample}];
  replicates=[replicates repmat(isample,[1,length(namelist_comb{isample})])];
  ppmvec=num2cell(mean(ppmlist_comb{isample},1));
  ppmadd=cellfun(@(x) num2str(x),ppmvec,'UniformOutput',false);
  ppm_all=[ppm_all ppmadd];
end
%
names(strcmp(names,'glucose-1-phosphate'))={'g1p'};
names(strcmp(names,'Valine'))={'valine'};
%
ymat=ymat_raw;
nridges=size(ymat,2);
ntime=52;
nDer=1;%maximal derivative to consider
xvec=1:ntime;
%% center and scale for each feature
for j=1:nridges
  shiftv=ymat(:,j)-mean(ymat(:,j));
  stdval=std(shiftv);
  if stdval==0.0
    stdval=1000;%practially ignore ridges with no changes during the selected time range
  end
  ymat(:,j)=shiftv./stdval;
end
% Funcitonal/smoothing representation for each feature
loglambda_vec= -4:0.25:4;
res=smooth_derivative(ymat,xvec,loglambda_vec,nDer);
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
names_unique=unique(names);
ncompd=length(names_unique);
%
annotnames_ind=find(~cellfun(@(x) any(regexp(x,'[uU]nknown')),names,'UniformOutput',true));
markerlist={'o' 'o' 'o' '^'};%for each conditions high vs low density
names_unique_pres={'ethanol', 'adenosine', 'alanine', 'g1p', 'trehalose', 'choline'};%instead of unique(names) just select a few to show
unsele_ind=cellfun(@(x) ~ismember(x,names_unique_pres),names,'UniformOutput',true);
names_pres=names;
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
  for replicate=unique(replicates)
    ind=find(replicates==replicate&strcmp(names_pres,name));
    if strcmp(name,'unknown')
      scatter(harmscr(ind,1),harmscr(ind,2),[],colorsmap(namei,:),markerlist{replicate});%
    else
      scatter(harmscr(ind,1),harmscr(ind,2),200,colorsmap(namei,:),markerlist{replicate},'filled');
    end
    legends=[legends ['experiment ' expnames{replicate} ' ' name]];
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
selecompds=['glucose', 'QA', names_unique_pres];
markerlist={'o' 'o' 'o' '^'};%for each conditions high vs low density
colorline={'r' 'r' 'r' 'g'};
for compd=selecompds
  legendsrec={};
  h2=figure();
  hold on;
  for replicate=unique(replicates)
    expname=expnames{replicate};
    ind=find(replicates==replicate&strcmp(names,compd));
    yplot=ymat(:,ind);
    yplot=yplot(:);
    xplot=repmat(timelist_comb{replicate}(:,1),[length(ind),1]);
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
