%%fda pca to show distribution of different time series curves
%% the PCA is done on each sample and different ridges is treated as different samples
close all;
clear all;
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
