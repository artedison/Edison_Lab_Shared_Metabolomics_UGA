% Use correlation information to construct network for annotation
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
workdir=[pardir 'result/annotation/'];%%the working folder
datadir=[pardir 'data/'];%the data folder
continue_dir=[pardir 'result_data/'];
cd(workdir);

%% correlation based network
threhold_corr=0.9;
load([continue_dir 'network_data.mat'])
cor_intensity=corr(mat_reshape_all,mat_reshape_all,'Type','Spearman');
indmat=zeros(size(cor_intensity));
lenmat=size(cor_intensity,1);
for i=1:lenmat
  indmat(i,(i+1):lenmat)=1;
end
triaind=find(indmat);
cor_intensity_triag_vec=cor_intensity(triaind);
cor_intensity_triag=cor_intensity;
cor_intensity_triag(~indmat)=0;
%%edge threhold
thres_display_vec_inten=[quantile(cor_intensity_triag_vec,threhold_corr)];

%%%positive correlation of intensity
[rowind colind]=find(cor_intensity_triag>=thres_display_vec_inten);
node1={};
node2={};
correlation=[];
for elei=1:length(rowind)
  roweleind=rowind(elei);
  coleleind=colind(elei);
  node1{elei}=namesall{roweleind};
  node2{elei}=namesall{coleleind};
  correlation(elei)=cor_intensity_triag(roweleind,coleleind);
end
cornet_table_corr=table(node1',node2',correlation','VariableNames',{'source' 'target' 'association'});
writetable(cornet_table_corr,'causalkinetix_corr.txt','Delimiter','\t');
