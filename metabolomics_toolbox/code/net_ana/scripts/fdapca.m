%%fda pca to show distribution of different time series curves
%% the PCA is done on each sample and different ridges is treated as different samples
close all;
clear all;
% ADD PATH
% Metabolic toolbox toolbox found @  https://github.com/artedison/Edison_Lab_Shared_Metabolomics_UGA
localPaths.public_toolbox='/Users/yuewu/Documents/GitHub/Edison_Lab_Shared_Metabolomics_UGA/';
% some wrappers for fdaM @ https://github.com/mikeaalv/fda_learn
localPaths.fdalearn='/Users/yuewu/Documents/GitHub/fda_learn/';
% functiona data analysis matlab package fdaM @ https://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/
localPaths.fdam='/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/matalb.lib/fdaM/';
%
addpath(genpath(localPaths.public_toolbox));
addpath(genpath(localPaths.fdalearn));
addpath(genpath(localPaths.fdam));
%
comp='/Users/yuewu/';%the computer user location
pardir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge.net/result_reprod/'];
% OR
% pardir=[pwd() '/'];
workdir=[pardir 'result/smoothing/'];%%the working folder
datadir=[pardir 'data/'];%the data folder
continue_dir=[pardir 'result_data/'];
cd(workdir);
load([continue_dir 'tracing.smoothed.mat']);
rng(1);
sampeseq=[1 4];
npc=5;%maximum number of pc to calculate
lambdapc=0.0;%lambda for pc
plot_ylim=[-2.2,2.2];
for isample=sampeseq
  ymat=matlist{isample};
  names=namelist{isample};
  timevec=timelist{isample}(:,1);
  nridges=size(ymat,2);
  ntime=52;
  nDer=1;%maximal derivative to consider
  xvec=1:ntime;
  %% center and scale for each feature (before pca)
  for j=1:nridges
    shiftv=ymat(:,j)-mean(ymat(:,j));
    ymat(:,j)=shiftv./std(shiftv);
  end
  loglambda_vec= -4:0.25:4;
  res=smooth_derivative(ymat,xvec,loglambda_vec,nDer);
  fdres=res.spfd_selec;
  nsmooth=nDer+2;
  norder=nsmooth+2;
  nbasis=length(xvec)+norder-2;
  bbasis=create_bspline_basis([min(xvec) max(xvec)],nbasis,norder,xvec);
  fdParpc=fdPar(bbasis,nsmooth,lambdapc);
  close all;
  res=fda_pca(fdres,fdParpc,npc);
  harmscr=res.fdapcastr.harmscr;
  annotnames_ind=find(~cellfun(@(x) any(regexp(x,'[uU]nknown')),names,'UniformOutput',true));
  colors=repmat([0 0 0],[nridges,1]);
  colors(annotnames_ind,1)=1;
  h2=figure();
    scatter(harmscr(:,1),harmscr(:,2),[],colors,'filled');
    xlabel('pc1');
    ylabel('pc2');
  %% add names
  text(harmscr(annotnames_ind,1),harmscr(annotnames_ind,2),names(annotnames_ind),'Color','red')
  saveas(h2,[workdir,'fdapca_score_plot.',num2str(isample),'.fig']);
  close all;
  plot_pca_fd(res.fdapcastr,1,[1 2]);
  figh=gcf;
  allaxes=findall(figh,'type','axes')
  set(allaxes(1),'ylim',plot_ylim);
  set(allaxes(2),'ylim',plot_ylim);
  saveas(figh,[workdir,'fdapca_loading_plot.',num2str(isample),'.fig']);
  close all;
  fig_prop=figure();
  plot(1:npc,res.fdapcastr.varprop(1:npc),'-o');
  xlabel('Eigenvalue Number');
  ylabel('Eigenvalue');
  saveas(fig_prop,[workdir,'fdapca_prop_plot.',num2str(isample),'.fig']);
  close all;
end

%%fda pca on 1sd derivative
%% the PCA is done on each sample and different ridges is treated as different samples
lambdapc=0;%lambda for pc
for isample=sampeseq
  ymat=smooth_1d_list{isample};
  names=namelist{isample};
  timevec=timelist{isample}(:,1);
  nridges=size(ymat,2);
  ntime=52;
  nDer=0;%maximal derivative to consider
  xvec=1:ntime;
  %% center and scale for each feature (before pca)
  for j=1:nridges
    shiftv=ymat(:,j)-mean(ymat(:,j));
    ymat(:,j)=shiftv./std(shiftv);
  end
  loglambda_vec= -4:0.25:4;
  res=smooth_derivative(ymat,xvec,loglambda_vec,nDer);
  fdres=res.spfd_selec;
  nsmooth=nDer+2;
  norder=nsmooth+2;
  nbasis=length(xvec)+norder-2;
  bbasis=create_bspline_basis([min(xvec) max(xvec)],nbasis,norder,xvec);
  fdParpc=fdPar(bbasis,nsmooth,lambdapc);
  res=fda_pca(fdres,fdParpc,npc);
  close all;
  harmscr=res.fdapcastr.harmscr;
  annotnames_ind=find(~cellfun(@(x) any(regexp(x,'[uU]nknown')),names,'UniformOutput',true));
  colors=repmat([0 0 0],[nridges,1]);
  colors(annotnames_ind,1)=1;
  h2=figure();
    scatter(harmscr(:,1),harmscr(:,2),[],colors,'filled');
    xlabel('pc1');
    ylabel('pc2');
  %% add names
  text(harmscr(annotnames_ind,1),harmscr(annotnames_ind,2),names(annotnames_ind),'Color','red')
  saveas(h2,[workdir,'fdapca_score_plot_1d.',num2str(isample),'.fig']);
  close all;
  plot_pca_fd(res.fdapcastr,1,[1 2]);
  figh=gcf;
  allaxes=findall(figh,'type','axes')
  set(allaxes(1),'ylim',plot_ylim);
  set(allaxes(2),'ylim',plot_ylim);
  saveas(figh,[workdir,'fdapca_loading_plot_1d.',num2str(isample),'.fig']);
  close all;
  fig_prop=figure();
  plot(1:npc,res.fdapcastr.varprop(1:npc),'-o');
  xlabel('Eigenvalue Number');
  ylabel('Eigenvalue');
  saveas(fig_prop,[workdir,'fdapca_prop_plot_1d.',num2str(isample),'.fig']);
  close all;
end

% FPCA aerobic vs anaerobic
% merged fpca
new_workdir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge.net/result_addition/fpca_multi_aero_vs_anerobic/'];
cd(new_workdir);
rng(1);
npc=5;%maximum number of pc to calculate
lambdapc=0.0;%lambda for pc
plot_ylim=[-2.2,2.2];
sampeseq=[1 2 3 4 5 6];
condseq=[1 1 1 2 2 2];
ymat=[];
names={};
replicates=[];
for isample=sampeseq
  ymat=[ymat matlist{isample}];
  names=[names namelist{isample}];
  replicates=[replicates repmat(condseq(isample),[1,length(namelist{isample})])];
end
%
nridges=size(ymat,2);
ntime=52;
nDer=1;%maximal derivative to consider
xvec=1:ntime;
%% center and scale for each feature (before pca)
for j=1:nridges
  shiftv=ymat(:,j)-mean(ymat(:,j));
  ymat(:,j)=shiftv./std(shiftv);
end
loglambda_vec= -4:0.25:4;
res=smooth_derivative(ymat,xvec,loglambda_vec,nDer);
fdres=res.spfd_selec;
nsmooth=nDer+2;
norder=nsmooth+2;
nbasis=length(xvec)+norder-2;
bbasis=create_bspline_basis([min(xvec) max(xvec)],nbasis,norder,xvec);
fdParpc=fdPar(bbasis,nsmooth,lambdapc);
close all;
res=fda_pca(fdres,fdParpc,npc);
harmscr=res.fdapcastr.harmscr;
annotnames_ind=find(~cellfun(@(x) any(regexp(x,'[uU]nknown')),names,'UniformOutput',true));
markerlist={'o' '^'};%for each conditions aerobic vs anaerobic
names_unique_pres={'glucose', 'ethanol', 'glucose-1-phosphate', 'arginine', 'alanine'};%instead of unique(names) just select a few to show
unsele_ind=cellfun(@(x) ~ismember(x,names_unique_pres),names,'UniformOutput',true);
names_pres=names;
names_pres(unsele_ind)={'unknown'};
names_unique_all=['unknown' names_unique_pres];
ncompd=length(names_unique_all);
colorsmap=hsv(ncompd);
colorsmap=colorsmap(randperm(size(colorsmap,1)),:);
color_noanno_ind=cellfun(@(x) any(regexp(x,'[uU]nknown')),names_unique_all,'UniformOutput',true);
colorsmap(color_noanno_ind,:)=repmat([0.5 0.5 0.5],[1,1]);
% comparing compounds
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
    legends=[legends ['replicate ' num2str(replicate) ' ' name]];
  end
end
xlabel('pc1');
ylabel('pc2');
title(['pca epl']);
legend(legends);
saveas(h2,[new_workdir,'fdapca_score_plot.aero_vs_anaero.fig']);
close all;

% plot for selected compounds
markerlist={'o' '^'};
colorline={'r','g'};
for compd=names_unique_pres
  legendsrec={};
  h2=figure();
  hold on;
  for replicate=unique(replicates)
    % expname=expnames{replicate};
    ind=find(replicates==replicate&strcmp(names_pres,compd));
    yplot=ymat(:,ind);
    yplot=yplot(:);
    xplot=repmat(timelist{replicate}(:,1),[length(ind),1]);
    scatter(xplot,yplot,200,colorline{replicate},markerlist{replicate});
    for linei=1:length(ind)
      lineind=(((linei-1)*ntime)+1):((linei*ntime));
      line(xplot(lineind),yplot(lineind),'LineWidth',2,'LineStyle','--','Color',colorline{replicate});
    end
    % legendsrec=[legendsrec {[expname]}];
  end
  xlabel('time');
  ylabel('quantification');
  title([compd]);
  % legend(legendsrec);
  saveas(h2,[workdir,'time_series_plot.',compd{1},'_in_allsample.fig']);
  close all;
end
