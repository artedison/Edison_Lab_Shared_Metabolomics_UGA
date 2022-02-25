%% analysis with the pyruvate dataset
set(0,'DefaultFigureVisible','on');
close all;
clear all;
comp='/Users/yuewu/';%the computer user location
pardir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge.net/result_addition/pyruvate_data/res/analysis/'];
workdir=[pardir];%%the working folder
datadir=[pardir];%the data folder
cd(workdir);
load(['pyruvate_exp.mat']);%the ridge tracking result as input
load(['tracing.smoothed.mat']);
%
expnames={'glucose','glucose','glucose','pyruvate (6.02 mg)','pyruvate (10 mg)','pyruvate (6.02 mg)','pyruvate (10 mg)'};
exptypes={'H','H','H','H','H','C13','C13'};
aerobic_ind=[1 2 3];
matlist_comb_new=[matlist(aerobic_ind) matlist_comb{1}(1) matlist_comb{1}(2) matlist_comb{2}(1) matlist_comb{2}(2)];
namelist_comb_new=[namelist(aerobic_ind) namelist_comb{1}(1) namelist_comb{1}(2) namelist_comb{2}(1) namelist_comb{2}(2)];
ppmlist_comb_new=[ppmlist(aerobic_ind) ppmlist_comb{1}(1) ppmlist_comb{1}(2) ppmlist_comb{2}(1) ppmlist_comb{2}(2)];
timelist_comb_new=[timelist(aerobic_ind) timelist_comb{1}(1) timelist_comb{1}(2) timelist_comb{2}(1) timelist_comb{2}(2)];
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
ppmmat=[];
ind_shift_vec=[0];
% whole data matrix
for isample=sampseq_all
  ymat_raw=[ymat_raw matlist_comb_new{isample}];
  ppmmat=[ppmmat ppmlist_comb_new{isample}];
  names=[names namelist_comb_new{isample}];
  replicates=[replicates repmat(isample,[1,length(namelist_comb_new{isample})])];
  ppmvec=num2cell(mean(ppmlist_comb_new{isample},1));
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
markerlist={'o' 'o' 'o' '^' '^' '^' '^'};%
filledlist=[true true true true true false false];
names_unique_pres={'uracil','alanine','choline','glucose','pyruvate','ethanol'}; %instead of unique(names) just select a few to show
unsele_ind=cellfun(@(x) ~ismember(x,names_unique_pres),names,'UniformOutput',true);
names_pres=names;
names_pres(unsele_ind)={'unknown'};
names_unique_all=['unknown' names_unique_pres];
ncompd=length(names_unique_all);
%
colorsmap=hsv(ncompd);
colorsmap=colorsmap(randperm(size(colorsmap,1)),:);
color_noanno_ind=cellfun(@(x) any(regexp(x,'[uU]nknown')),names_unique_all,'UniformOutput',true);
colorsmap(color_noanno_ind,:)=repmat([0.5 0.5 0.5],[length(find(color_noanno_ind)),1]);
legends={};
h2=figure();
hold on;
for namei=1:length(names_unique_all)
  name=names_unique_all{namei};
  for replicate=unique(replicates)
    ind=find(replicates==replicate&strcmp(names_pres,name));
    if length(ind)==0
      continue;
    end
    if strcmp(name,'unknown')
      scatter(harmscr(ind,1),harmscr(ind,2),[],colorsmap(namei,:),markerlist{replicate});%
    else
      if filledlist(replicate)
        scatter(harmscr(ind,1),harmscr(ind,2),200,colorsmap(namei,:),markerlist{replicate},'filled');
      else
        scatter(harmscr(ind,1),harmscr(ind,2),200,colorsmap(namei,:),markerlist{replicate});
      end
    end
    legends=[legends [expnames{replicate} ' ' exptypes{replicate} ' ' name]];
  end
end
xlabel('pc1');
ylabel('pc2');
title(['pca epl']);
legend(legends);
saveas(h2,[workdir,'fdapca_score_plot.fig']);
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
selecompds=[names_unique_pres];
markerlist={'o' 'o' 'o' '^' '^' '^' '^'};%glucose vs pyruvate
% colorline={'r' 'r' 'r' 'g' 'g' 'g' 'g'};
colorline={'r' 'r' 'r' 'r' 'r' 'g' 'g'};%C12 vs C13
% filledlist=[true true true true true false false];
for compd=selecompds
  legendsrec={};
  h2=figure();
  hold on;
  for replicate=unique(replicates)
    expname=expnames{replicate};
    ind=find(replicates==replicate&strcmp(names,compd));
    yplot=ymat(:,ind);
    yplot=yplot(:);
    xplot=repmat(timelist_comb_new{replicate}(:,1),[length(ind),1]);
    scatter(xplot,yplot,200,colorline{replicate},markerlist{replicate});
    if strcmp(colorline{replicate},'g')
      legendsrec=[legendsrec {[expname ' C13']}];
    else
      legendsrec=[legendsrec {[expname ' C12']}];
    end
    for linei=1:length(ind)
      lineind=(((linei-1)*ntime)+1):((linei*ntime));
      line(xplot(lineind),yplot(lineind),'LineWidth',2,'LineStyle','--','Color',colorline{replicate});
      legendsrec=[legendsrec {['']}];
    end
  end
  xlabel('time');
  ylabel('quantification');
  title([compd]);
  legend(legendsrec);
  saveas(h2,[workdir,'time_series_plot.',compd{1},'_in_allsample.fig']);
  close all;
end

% plot C12 vs C13 quantification of the same compound
% normalize by pyruvate quantification
quant_list_c13={'alanine','ethanol'}
ref_compd='pyruvate';
ref_ratio_low=0.1;
% calculate the scale factor
combexps={'pyruvate (6.02 mg)','pyruvate (10 mg)'};
scal_factor=[];
for combexp=combexps
  replicate_comb=find(strcmp(expnames,combexp));
  % C12 parts
  ind=find(replicates==replicate_comb(1) & strcmp(names,ref_compd));
  py_sum_vec_c12=sum(ymat_raw(:,ind),2);
  % C13 parts
  ind=find(replicates==replicate_comb(2) & strcmp(names,ref_compd));
  py_sum_vec_c13=ymat_raw(:,ind);
  % quantile range
  c12cutoff_ind=max(find(py_sum_vec_c12>(max(py_sum_vec_c12)*ref_ratio_low)));
  c13cutoff_ind=max(find(py_sum_vec_c13>(max(py_sum_vec_c13)*ref_ratio_low)));
  seleind=min(c12cutoff_ind,c13cutoff_ind)
  seleseq=1:seleind;
  scal_factor=[scal_factor mean(py_sum_vec_c12(seleseq)./py_sum_vec_c13(seleseq))];
end
% plot for compounds that have both C12 and C13 parts
% 1. convert C13 quantification of the same level as C12 spectra, 2. scale by maximum, 3. scale each peak by maximum and to the whole maximum
colorline={'r' 'r' 'r' 'r' 'r' 'g' 'g'};%C12 vs C13
markerlist={'o' 'o' 'o' 'o' '^' 'o' '^'};%the two experiments
pyexp_rep=[4 5 6 7];
ethanol_quan=[1.1 1.3];
for compd=quant_list_c13
  h2=figure();
  hold on;
  glob_fac=[];
  legends={};
  for replicate=pyexp_rep
    expname=expnames{replicate};
    ind=find(replicates==replicate&strcmp(names,compd));
    if strcmp(compd,'ethanol')
      ppmmean=mean(ppmmat(:,ind),1);
      ind=ind(ppmmean>ethanol_quan(1)&ppmmean<ethanol_quan(2));
    end
    yplot=ymat_raw(:,ind);
    dataty='c13';
    if replicate==6
      yplot=yplot*scal_factor(1)/glob_fac(1);
    elseif replicate==7
      yplot=yplot*scal_factor(2)/glob_fac(2);
    else
      elefac=max(yplot(:));
      yplot=yplot/elefac;
      glob_fac=[glob_fac elefac];
      dataty='c12';
    end
    yplot=yplot./max(yplot,[],1)*max(yplot(:));
    yplot=yplot(:);
    xplot=repmat(timelist_comb_new{replicate}(:,1),[length(ind),1]);
    scatter(xplot,yplot,200,colorline{replicate},markerlist{replicate});
    legends=[legends {[expname dataty]}];
    for linei=1:length(ind)
      lineind=(((linei-1)*ntime)+1):((linei*ntime));
      line(xplot(lineind),yplot(lineind),'LineWidth',2,'LineStyle','--','Color',colorline{replicate});
      legends=[legends {['']}];
    end
  end
  xlabel('time');
  ylabel('quantification');
  title([compd]);
  legend(legends);
  saveas(h2,[workdir,'time_series_plot.',compd{1},'_in_allsample_c12_vs_c13.fig']);
  close all;
end
