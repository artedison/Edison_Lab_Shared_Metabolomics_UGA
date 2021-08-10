%% Formualte the input for CausalKinetiX (and for all related analysis that apply to features existing in all replicates)
%% Combine ridges from different replicates and annotate compound peaks with annotation
%% The combination is an automatic process and visual check is provided at the end.

close all;
clear all;

%% the user will need to modify the path here for local run
comp='/Users/yuewu/';%the computer user location
pardir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/ridge.net/result_reprod/'];
% OR
% pardir=[pwd() '/'];
workdir=[pardir 'result/quality_check/'];%%the working folder
datadir=[pardir 'data/'];%the data folder
continue_dir=[pardir 'result_data/'];
cd(workdir);
load([continue_dir 'tracing.newmeth.experiment.manual.completerid.mat']);
load([continue_dir 'tracing.smoothed.mat']);
load([datadir 'unshared/sampleData.mat']);%please refer to the CIVM-NMR paper for this data set: Continuous in vivo Metabolism by NMR

%%%formulate the data for input of CausalKinetiX
%%% mat_all matched peaks existing in all 6 replicates
%%%% the matching is based on mean ppm and iterated for all replicate from 1-6 [1 2 3 7 8 9]. no manual criteria is involved.

ntimes=52;
timevec=sampleData(1).timesCollapsed_1h1d;
timeindvec=1:ntimes;
deltapm_threshold=0.01;%distance threshold for peak matching
wind_size_view=0.01;%%in term of ppm
sample_ind_match=[1 2 3 7 8 9];
compoundnames=unique([namelist{:}]);
mat_all=[];
mat_reshape_all=[];
ppmmatch_ind_all=[];
namesall={};
nsample=length(namelist);
% visual parameters
fig_layout=[4,7];
fig_plot_i=1;
figind=0;
region_sep_dist=0.05;
fig=figure();
for compdname=compoundnames
  compdname=compdname{1};
  n1names=length(find(cellfun(@(x) strcmp(x,compdname),namelist{1},'UniformOutput',true)));
  %% find out ppms matches use order of tracked peaks in ppm for known peaks and ppm distance for unknown peaks
  if strcmp(compdname,'unknown')
    [ppmmatch_ind ppmvec]=ppm_list_match(ppmlist,namelist,'^unknown$',deltapm_threshold);
  else
    %% for compound peak match it is only by ppm order and this is possible because those sample without corresponding compound peaks are totally missing all peaks for the compound
    %% it is checked by ploting here
    ppmmatch_ind=[];
    ppmmatch_ppm_array=[];
    ppmvec=[];
    for isample=1:nsample
      namevec=namelist{isample};
      ppmmat=ppmlist{isample};
      matchind=cellfun(@(x) strcmp(x,compdname),namevec,'UniformOutput',true);
      unanno_ind=find(matchind);
      ppmmean=mean(ppmmat,1);
      ppmmean_unannot=ppmmean(unanno_ind);
      [ppmmatch_ppm,sortind]=sort(ppmmean_unannot);
      ppmmatch_ind=[ppmmatch_ind; unanno_ind(sortind)];
      ppmmatch_ppm_array=[ppmmatch_ppm_array; ppmmatch_ppm];
    end
    ppmvec=mean(ppmmatch_ppm_array,1);
  end
  %%the intensity matrix
  mat_temp=[];
  mat_reshape_temp=[];
  namevec_temp=strcat({compdname},cellstr(num2str(ppmvec'))');
  for isample=1:size(ppmmatch_ind,1)%%some ridges will not be trcked in all samples
    intenmat=matlist{isample};
    ppmmat=ppmlist{isample};
    ppmind=ppmmatch_ind(isample,:);
    intenarray=intenmat(:,ppmind);
    mat_temp=[mat_temp; intenarray(:)'];
    mat_reshape_temp=[mat_reshape_temp; intenarray];
  end
  %%plotting check
  display([compdname,' ',num2str(length(namevec_temp)) '/' num2str(n1names)]);
  % separate ridge regions
  endind_vec=[0];
  ppmi=2;
  startppm=ppmvec(1);
  while ppmi<=length(ppmvec)
    if abs(ppmvec(ppmi)-startppm) >= region_sep_dist
      endind_vec=[endind_vec ppmi-1];
      startppm=ppmvec(ppmi);
    end
    ppmi=ppmi+1;
  end
  endind_vec=[endind_vec length(ppmvec)];
  for endind_i=2:length(endind_vec)
    range_region=[endind_vec(endind_i-1)+1 endind_vec(endind_i)];
    ppmmatch_ind_plot=ppmmatch_ind(:,range_region(1):range_region(2));
    % plot surface for each sample
    intensi_comb=[];
    for samplei=1:size(ppmmatch_ind_plot,1)%%some ridges will not be tracked in all samples
      surf_ind=samplei+7*(fig_plot_i-1);
      data=sampleData(sample_ind_match(samplei));
      mathere=data.Xcollapsed_1h1d;
      ppmhere=data.ppm_1h1d;
      %%%decide window size based on ridge ppm range
      ppmridind=ppmmatch_ind_plot(samplei,:);
      nridges=length(ppmridind);
      ppmmat=ppmlist{samplei};
      intenmat=matlist{samplei};
      ridppm=ppmmat(:,ppmridind);
      ppmrange=[min(ridppm(:))-wind_size_view max(ridppm(:))+wind_size_view];
      ppmrangeind=matchPPMs(ppmrange,ppmhere);
      ppmrangeseq=ppmrangeind(1):ppmrangeind(2);
      ppmvecloc=ppmhere(ppmrangeseq);
      intensityarray=mathere(:,ppmrangeseq);
      %%
      ridlen=ntimes;
      clusters=1:nridges;
      clustshere=repmat(clusters,[ridlen,1]);
      clustshere=clustshere(:);
      cindallhere=matchPPMs(ridppm(:),ppmhere)-ppmrangeind(1)+1;
      rindallhere=repmat(timeindvec,[1,nridges]);
      ridvalallhere=intenmat(timeindvec,ppmridind);
      ridvalallhere=ridvalallhere(:);
      intensi_comb=[intensi_comb; ridvalallhere'];
      titlehere=[compdname ' ' num2str(ppmrange(1)) '-' num2str(ppmrange(2))];
      sub=subplot_tight(fig_layout(1),fig_layout(2),surf_ind);
      plotRidgesherenew(intensityarray,ppmvecloc,timevec,clustshere,cindallhere,rindallhere,ridvalallhere,clusters,titlehere);
      subfig=gca;
      subfig_hand=gcf;
      copyobj(get(subfig,'children'),sub);
      close(subfig_hand);
      view([0 45]);
      set(gca,'XDir','reverse');
    end
    % plot time trajectory
    % concatenate multiple ridges together for easy visuzlization
    timevec_comb=repmat(1:size(intensi_comb,2),[size(intensi_comb,1),1]);
    subplot_tight(fig_layout(1),fig_layout(2),7*(fig_plot_i-1)+7);
    plot(timevec_comb',intensi_comb');
    title(titlehere);
    if fig_plot_i==fig_layout(1) | (strcmp(compdname,compoundnames{end}) & endind_i==length(endind_vec))
      saveas(fig,['comb_quality_check ' num2str(figind) '.fig']);
      close(fig);
      fig=figure();
      fig_plot_i=1;
      figind=figind+1;
    else
      fig_plot_i=fig_plot_i+1;
    end
  end
  %%combine data
  ucrlen=size(ppmmatch_ind,1);
  if ucrlen<6
    mat_temp=[mat_temp; zeros([6-ucrlen,size(mat_temp,2)])];
    mat_reshape_temp=[mat_reshape_temp; zeros([(6-ucrlen)*ntimes,size(mat_reshape_temp,2)])];
    ppmmatch_ind=[ppmmatch_ind; NaN([6-ucrlen,size(ppmmatch_ind,2)])];
  end
  mat_all=[mat_all mat_temp];
  mat_reshape_all=[mat_reshape_all mat_reshape_temp];
  ppmmatch_ind_all=[ppmmatch_ind_all ppmmatch_ind];
  namesall=[namesall namevec_temp];
end
envvec=[1 1 1 2 2 2];
save('CausalKinetiX_input_peakwise.mat','mat_all','timevec','envvec','namesall','ppmmatch_ind_all');
save('network_data.mat','namesall','mat_reshape_all','ppmmatch_ind_all');
% 'mat_all': feature intenstiy matrix (nreplicate*(nfeature*ntime))
% 'timevec': time vector (ntime)
% 'envvec': environment vector (nreplicate)
% 'namesall': compound name vector (nfeature)
% 'ppmmatch_ind_all': ppm match matrix (nreplicate*(nfeature))
% 'mat_reshape_all': feature intenstiy matrix ((nreplicate*ntime)*nfeature)
