function [flag_run]=spec_vis_gissmo(id,sampi,matchnode,paralist)
% visualize the spectra of tracked ridges with the gissmo reference spectra. This function is often used with the R script network_clust_annotation.R. Annotated peaks will be highlighted by stars while unkown peaks will be highlighted by circles.
%
% Arugment:
%         id: string. the BMSE id for the matched compound. Must be provided
%         sampi: integer. index of the sample to visualize with. Must be provided
%         matchnode: array. string. the matched nodes for the ridges. Must be provided.
%         paralist: structure. add on parameters. Default nan. If changes is needed for data folder please provide them.
%                   datadir: the path to spectra data
%                   gissmodir: the path to gissmo data
%                   match_dir: the path to local working directory
% Return:
%         flag_run: the function has been run properly 0, otherwise 1(the compound wasn't matched in the gissmo library).
%
% Examples:
%
% matchnode={'unknown     2.1218','unknown     2.3687','glutamate2.3822','unknown     2.1113','unknown     2.0986','unknown     4.1845','unknown     2.1338','unknown     2.3621','unknown      2.144'};
% id='bmse000913';%glutamate spectra
% sampi=1;
% paralist=struct();
% paralist.datadir='/PATH/TO/SPECTRA/DATA/';%% the user can download from links to our previous publication as stated in the manuscript
% paralist.gissmodir='/PATH/TO/GISSMO/LIB/';
% paralist.match_dir='/PATH/TO/CURRENT/WORKING/DIR/WITH/DATA/';
% spec_vis_gissmo(id,sampi,matchnode,paralist);
%
% Yue Wu 07/02/2020
% Tested with MATLAB R2018b

if ~exist('id','var')
  error('please provide input BMSE id');
end
if ~exist('sampi','var')
  error('please provide input sample id');
end
if ~exist('matchnode','var')
  error('please provide input matched ppm');
end
if ~exist('paralist','var')
  paralist=struct();
end
if ~isfield(paralist,'datadir')%location of spectra data
  datadir='.';
else
  datadir=paralist.datadir;
end
if ~isfield(paralist,'gissmodir')%location of gissmo library
  gissmo_dir='.';
else
  gissmo_dir=paralist.gissmodir;
end
if ~isfield(paralist,'match_dir')%location of ridge result
  match_dir='.';
else
  match_dir=paralist.match_dir;
end

% match compound name
load([gissmo_dir 'gissmo_tab.mat']);
load([gissmo_dir 'wholegissmo_spectral.complex.more.mat']);
indtab=find(strcmp(id,tabinfor{:,'EntryID'}));
matchcompd=tabinfor{indtab,'CompoundName'};%compound name
match_spec_id=[id '_simulation_1'];%use the information from the first simulation in gissmo for each compound
if any(cellfun(@(x) strcmp(match_spec_id,x),fieldnames(strdata),'UniformOutput',true))
  matchspec=strdata.(match_spec_id);%compound spectra
else
  warning('the id not matched in current local gissmo database');
  flag_run=1;
  return
end
matchppm=ppm;%compound ppm

% the experimental nmr spectra
load([match_dir 'network_data.mat']);%network related data
load([match_dir 'tracing.smoothed.mat']);%ridge tracking result
load([datadir 'sampleData.mat']);%spectra data

matchind=cellfun(@(x) find(strcmp(namesall,x)),matchnode,'UniformOutput',true);
anno_peak_mask= ~contains(matchnode,'unknown','IgnoreCase',true);
if size(matchind,2)~=size(matchnode,2)
  stop('unmatched nodes. probably indicating problem in name format');
end
locind=ppmmatch_ind_all(sampi,matchind);
% decide the time range
timeind=ridgerangelist{sampi};
timerange_arra=timeind(locind,:);
showrange=[max(timerange_arra(:,1)),min(timerange_arra(:,2))];
showind=[showrange(1) floor(mean(showrange)) showrange(2)];
ppmmat=ppmlist{sampi};
ppmmat_loc=ppmmat(showind,locind);
% matched spectra
exp_ppm=sampleData(sampi).ppm_1h1d;
exp_mat=sampleData(sampi).Xcollapsed_1h1d;
exp_mat=exp_mat(showind,:);

% plotting
fig=figure();
peak_inten_array=[];
peak_ind_array=[];
colorvec={'r','b','g'};
ploti=1;
for showi=1:length(showind)
  peakind=matchPPMs(ppmmat_loc(showi,:),exp_ppm);
  peak_ind_array=[peak_ind_array; peakind];
  peakintensity=exp_mat(showi,peakind);
  peak_inten_array=[peak_inten_array; peakintensity];
  [~, plines(showi)]=plotr(exp_ppm,exp_mat(showi,:),'Color',colorvec{showi});
  hold on;
  if any(~anno_peak_mask)
    pdots(ploti)=plot(exp_ppm(peakind(~anno_peak_mask)),peakintensity(~anno_peak_mask),'.','Color',colorvec{showi});
    ploti=ploti+1;
  end
  if any(anno_peak_mask)
    pdots(ploti)=plot(exp_ppm(peakind(anno_peak_mask)),peakintensity(anno_peak_mask),'s','Color',colorvec{showi});
    ploti=ploti+1;
  end
end
% rescale matched gissmo spectra
convfactor=max(peak_inten_array(:))/max(matchspec);
matchspec=matchspec*abs(convfactor);
[~, plines(4)]=plotr(matchppm,matchspec,'Color','k');
legend([plines(1) plines(2) plines(3) plines(4)],{'start','middle','end','reference'});
title(matchcompd);
set(gca,'fontsize',20);
flag_run=0;
