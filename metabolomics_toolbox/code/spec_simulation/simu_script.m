%% the script to simulate multiple sample of spectra with different concentration of different compounds
%% the simualtion will be composed of GISSMO spectra,
%% concentration differences in chemicals between samples can be introduced to formulate the conditions in experiments
%%%% YUE WU 12112019
close all;
clear all;
comp='/Users/yuewu/';%the computer user location
workdir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/spec_simu_gissmo/'];%%the working folder
libdir=[comp 'Dropbox (Edison_Lab@UGA)/Resources/gissmo_lib/cleandata/'];
cd(workdir);
rng(1);%%random seed
args=struct();%% all the arguments

%%%%%%%%%%%%%%USER ARGUMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%sample control
args.nsample=100;%%total number of samples
ngroup=2;
subsampnum=[50 50];%%please make sure sum(subsampnum)=args.nsample
args.sampleindvec=[repmat(1,1,subsampnum(1)) repmat(2,1,subsampnum(2))];%%use this to control samples of each group
%%concentration control
args.conc_range=[0 1000];%%the whole dynamic range of all compounds
args.totalnoiseregion=[-0.5 -0.2];%% the region of no signal to shift spectra
args.factor_std=0.3*[repmat(1,1,subsampnum(1)) 1.5*repmat(2,1,subsampnum(2))];%0.3;%%the relative variance for compound concentration in samples can also be a vector for indicate all the compound specifically
args.sigma_noise=15;%white noise on the whole spectra
%%type of simulation: 'random_chemical_sample' OR 'biological_sample'
args.simutype='biological_sample' %'biological_sample';
args.group_vec={'Light' 'Dark'}; %different biological groups
%%phenotype control
args.pheno_vec={'growth_rate' 'weight' 'offsprings'};
args.phenorang=[0 10];% the range of phenotype value
%% compound control
args.compd_vec={'ATP' 'Acetic-acid' 'Adenosine' 'Betaine' 'Choline' 'Citrate' 'Cytosine' 'DSS' 'Ethanol' 'Formate' 'Fumaric-acid' 'Glycerol' 'Glycine' 'L-Serine' 'L-Threonine' 'L-Valine' 'L-alanine' 'L-arginine' 'L-glutamic-acid' 'L-malic-acid'};
%%concentration that are changed between groups (can have no this field)
%%%{compound, relative mean concentration in different samples}
args.modif_conc_list=cell2table({'L-arginine' 1 3; 'Choline' 2 1});
%%correlation between compound and phenotypes (can have no this field)
%%%group index, compound, phenotype index, simulated correlation
args.corr_conc_list=cell2table({2 'Glycine' 1 0.8});
%%technical control
args.lambda=-0.00035;% the exponential window function for fid, need to be modified according to inspection of result spectra
args.exitprob=1;% only used in 'random_chemical_sample' the compound can be totally missing from the spectra if 1 (by p=0.5). the listed compounds must be in the spetra if 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%sample group name vector
args.group_ident_vec=args.group_vec(args.sampleindvec);

%%% phenotypes simulation
%%% more complex ways of simulation can be provided or the user can provide the table themselves
%%%% phenotype value are uniform random generated. the user can provide their matrix as well
pheno_mat=zeros(length(args.sampleindvec),length(args.pheno_vec));
pheno_mat=rand(length(args.sampleindvec),length(args.pheno_vec))*(args.phenorang(2)-args.phenorang(1))+args.phenorang(1);
args.pheno_mat=pheno_mat;
%%concentration simulation
%% OR the user can provide the concentration matrix themselves
simures=spec_conc_simu(args);
conc_mat=simures.conc_mat;

%%% spectra simulation
load([libdir 'wholegissmo_spectral.complex.more.mat']);
load([libdir 'wholegissmo_spectral.real.list.mat']);
load([libdir 'gissmo_tab.mat']);
spec_args=struct();
spec_args.ppm=ppm;
spec_args.sampleindvec=args.sampleindvec;
spec_args.compd_vec=args.compd_vec;
spec_args.tabinfor=tabinfor;
spec_args.conc_mat=conc_mat;
spec_args.strdata=strdata;
spec_args.sigma_noise=args.sigma_noise;
spec_args.lambda=args.lambda;
spec_args.totalnoiseregion=args.totalnoiseregion;
spec_mat=nmr_spec_simu(spec_args);

%%saving data
group_vec=args.group_vec;
sampleindvec=args.sampleindvec;
group_ident_vec=args.group_ident_vec;
group_ident_vec=args.group_ident_vec;
pheno_vec=args.pheno_vec;
compd_vec=args.compd_vec;
pheno_tab=array2table(pheno_mat,'VariableNames',args.pheno_vec);
save('simu_spect.mat','conc_mat','spec_mat','ppm','pheno_mat','pheno_tab','group_ident_vec','compd_vec','pheno_vec','group_vec','sampleindvec');

%%checks
% surf(conc_mat);
% heatmap(corrcoef(conc_mat));
% corrcoef(conc_mat(:,5),conc_mat(:,18))
% mean(conc_mat(51:100,5))/mean(conc_mat(1:50,5))%0.5
% mean(conc_mat(51:100,18))/mean(conc_mat(1:50,18))%3.0
% corrcoef(conc_mat(51:100,13),pheno_mat(51:100,1))%0.8
% plot(conc_mat(51:100,13),pheno_mat(51:100,1));

% X=spec_mat;
% Yvec=args.sampleindvec;
% % save('simu_spect_student.mat','X','ppm','Yvec');%%spec_mat: X sampleindvec: Yvec
% %plot test
% plotr(ppm,spec_mat(2,:));
% stackSpectra(spec_mat,ppm,0.0,10,'simulated spectra');
