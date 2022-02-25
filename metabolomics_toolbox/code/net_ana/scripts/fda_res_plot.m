% Some example explotary plot can be produced from the FDA smoothed data
% example spectra and correlation heatmap
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
%% the user will need to modify the path here for local run
comp='/Users/yuewu/';%the computer user location
pardir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge.net/result_reprod/'];
% OR
% pardir=[pwd() '/'];
workdir=[pardir 'result/smoothing/'];%%the working folder
datadir=[pardir 'data/'];%the data folder
continue_dir=[pardir 'result_data/'];
cd(workdir);
load([continue_dir 'tracing.smoothed.mat']);
%%plot smoothing result (randome selected examples)
nsampsize=16;
layout=[4 4];
for randi=[1 2 3 4 5 6]
  rng(randi);
  rangcheck=randsample(size(ymat_all,2),nsampsize);%%
  ypartmat=ymat_all(:,rangcheck);
  % a scale factor can be included or not
  meany=abs(mean(ypartmat,1));
  % stdy=std(ypartmat,1);
  meany=ones([1 length(meany)]);
  ypartmat_scal=ypartmat./meany;%./stdy
  %%plot original data
  for sampi=1:nsampsize
    sampind=rangcheck(sampi);
    xvectime=timemat_all(:,sampind);
    fdres=listres{sampind}.spfd_selec;
    xfinetime=linspace(min(xvectime),max(xvectime),201)';
    xfine=linspace(min(xvec),max(xvec),201)';
    intmatfine=eval_fd(xfine,fdres);%1
    intmatfine_scal=intmatfine/meany(sampi);%./stdy(sampi)
    % figure();
    subplot(layout(1),layout(2),sampi);
    phdl=plot(xfinetime,intmatfine_scal,'-');
    set(phdl,'LineWidth',2);
    hold on;
    plot(xvectime,ypartmat_scal(:,sampi),'o')
    hold off;
    xlabel('\fontsize{19} time')
    ylabel('\fontsize{19} intensity')
    locele=listloc{sampind};
    title([num2str(locele(1)) ' ' num2str(locele(2),'%4.3f') ' ' num2str(locele(3),'%4.3f')]);
  end
  fig=gcf;
  saveas(fig,[workdir,'examp_rand_',num2str(randi),'_original_2','.fig']);
  close(fig);
  %%plot first order derivative data
  for sampi=1:nsampsize
    sampind=rangcheck(sampi);
    xvectime=timemat_all(:,sampind);
    fdres=listres{sampind}.spfd_selec;
    xfinetime=linspace(min(xvectime),max(xvectime),201)';
    xfine=linspace(min(xvec),max(xvec),201)';
    intmatfine_d1=eval_fd(xfine,fdres,1);%1
    % figure();
    subplot(layout(1),layout(2),sampi);
    phdl=plot(xfinetime,intmatfine_d1,'-');
    set(phdl,'LineWidth',2);
    xlabel('\fontsize{19} time')
    ylabel('\fontsize{19} intensity')
    locele=listloc{sampind};
    title([num2str(locele(1)) ' ' num2str(locele(2),'%4.3f') ' ' num2str(locele(3),'%4.3f')]);
  end
  fig=gcf;
  saveas(fig,[workdir,'examp_rand_',num2str(randi),'_d1_2','.fig']);
  close(fig);
end

% plot a few targeted compounds
compdlist={'glucose-1-phosphate' 'choline' 'glucose' 'uridine' 'alanine' 'lactate' 'ethanol' 'phenylalanine' 'trehalose' 'tyrosine' 'arginine'};
sampi=1;
checked_ind=[];
for compdelei=1:length(compdlist)
  compdele=compdlist{compdelei};
  matchind=sort(find(strcmp(namelist{1},compdele)));
  checked_ind=[checked_ind matchind(1)];
end
checked_ind=[checked_ind 257];%addon plotting for selected unknown peaks
ypartmat=ymat_all(:,checked_ind);
% a scale factor can be included or not
meany=abs(mean(ypartmat,1));
% stdy=std(ypartmat,1);
meany=ones([1 length(meany)]);
ypartmat_scal=ypartmat./meany;%./stdy
%%plot original data
layout=[4 4];
for sampi=1:length(checked_ind)
  sampind=checked_ind(sampi);
  if sampi<=length(compdlist)
    compdname=compdlist{sampi};
  else
    compdname=['unknown'];
  end
  xvectime=timemat_all(:,sampind);
  fdres=listres{sampind}.spfd_selec;
  xfinetime=linspace(min(xvectime),max(xvectime),201)';
  xfine=linspace(min(xvec),max(xvec),201)';
  intmatfine=eval_fd(xfine,fdres);%1
  intmatfine_scal=intmatfine/meany(sampi);%./stdy(sampi)
  % figure();
  subplot(layout(1),layout(2),sampi);
  phdl=plot(xfinetime,intmatfine_scal,'-');
  set(phdl,'LineWidth',2);
  hold on;
  plot(xvectime,ypartmat_scal(:,sampi),'o')
  hold off;
  xlabel('\fontsize{19} time')
  ylabel('\fontsize{19} intensity')
  locele=listloc{sampind};
  title([num2str(locele(1)) ' ' num2str(locele(2),'%4.3f') ' ' num2str(locele(3),'%4.3f') '' compdname]);
end
fig=gcf;
saveas(fig,[workdir,'examp_sele_original','.fig']);
close(fig);
%%plot first order derivative data
for sampi=1:length(checked_ind)
  sampind=checked_ind(sampi);
  if sampi<=length(compdlist)
    compdname=compdlist{sampi};
  else
    compdname=['unknown'];
  end
  xvectime=timemat_all(:,sampind);
  fdres=listres{sampind}.spfd_selec;
  xfinetime=linspace(min(xvectime),max(xvectime),201)';
  xfine=linspace(min(xvec),max(xvec),201)';
  intmatfine_d1=eval_fd(xfine,fdres,1);%1
  % figure();
  subplot(layout(1),layout(2),sampi);
  phdl=plot(xfinetime,intmatfine_d1,'-');
  set(phdl,'LineWidth',2);
  xlabel('\fontsize{19} time')
  ylabel('\fontsize{19} intensity')
  locele=listloc{sampind};
  title([num2str(locele(1)) ' ' num2str(locele(2),'%4.3f') ' ' num2str(locele(3),'%4.3f') '' compdname]);
end
fig=gcf;
saveas(fig,[workdir,'examp_sele_d1_2','.fig']);
close(fig);

% check ppm range in ridge tracking
load([continue_dir 'tracing.newmeth.experiment.manual.completerid.mat'])
randi=1;
nridges=size(Sample_complete_rid(1).ridges,2);
rng(randi);
nsampsize=16;
sampinds=randsample(nridges,nsampsize);
fprintf('sample_id\tstart_range_ppm\tend_range_ppm\tmean_ppm\tpeak_in_range\n');
for sampind=sampinds'
  region=Sample_complete_rid(1).ridges(sampind+1).parameters.region;
  meanppm=mean(Sample_complete_rid(1).ridges(sampind+1).result.ppm);
  fprintf('%d:\t%f\t%f\t%f\t',sampind,region(1),region(2),meanppm);
  if meanppm>region(1)&meanppm<region(2)
    fprintf('true');
  else
    fprintf('false');
  end
  fprintf('\n');
end
