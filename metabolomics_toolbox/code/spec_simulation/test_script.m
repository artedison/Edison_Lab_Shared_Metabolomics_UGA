%% this script will
%%1. load the data and test statistical analysis on the data set
%% PCA and STOCKSY will be used for analyzing distinguishing compounds
%% the pair with high correlation between phenotypes and compound (peaks) will be found
%% No alignement, region removal, scaling, normalization should be needed
%%2. some condition test of the function spec_conc_simu()
%%%% yue wu
close all;
clear all;
comp='/Users/yuewu/';%the computer user location
% workdir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/Documents/ta/simu_exam/'];
workdir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/spec_simu_gissmo/'];
load([workdir 'simu_spect.mat']);
%differences
%% 'L-arginine' 1 3; 'Choline' 2 1
% spec_mat=X;
% sampleindvec=Yvec;
stackSpectra(spec_mat,ppm,0.0,10,'simulated spectra');
testregion=[3.7 3.8]%[3.1 3.3] [1.6 1.8]%% [1.8 2.0] [3.7 3.8]
ppmrang=matchPPMs(testregion,ppm);
ppmind=ppmrang(1):ppmrang(2);
mattest=spec_mat(:,ppmind);
ppmtest=ppm(ppmind);
samples=1:length(group_ident_vec);
fig=figure(), hold on
    surf(ppmtest,samples',mattest,'FaceColor','Interp');
    % surf(mattest);
    shading interp;
    set(gca,'xdir','reverse');
    ylabel('sample');
    zlabel('intensity');
    title(['shiftregion']);
    xlabel('ppm');

%correlation
%% {2 'Glycine' 1 0.8}
testregion=[3.54 3.55];
ppmrang=matchPPMs(testregion,ppm);
ppmind=ppmrang(1):ppmrang(2);
groupind=find(strcmp(group_ident_vec,'Dark'));
mattest=spec_mat(groupind,ppmind);
ppmtest=ppm(ppmind);
fig=figure(), hold on
    surf(ppmtest,groupind',mattest,'FaceColor','Interp');
    % surf(mattest);
    shading interp;
    set(gca,'xdir','reverse');
    ylabel('sample');
    zlabel('intensity');
    title(['shiftregion']);
    xlabel('ppm');
fig=figure(), hold on
    surf(ppmtest,samples',spec_mat(:,ppmind),'FaceColor','Interp');
    % surf(mattest);
    shading interp;
    set(gca,'xdir','reverse');
    ylabel('sample');
    zlabel('intensity');
    title(['shiftregion']);
    xlabel('ppm');
corr(max(mattest,[],2),table2array(pheno_tab(groupind,1)))
corr(table2array(pheno_tab(groupind,1)),conc_mat(groupind,13))
%PCA
% spec_mat=X;
% sampleindvec=Yvec;
normcheck(spec_mat)
% spec_mat_N=normalize(spec_mat,ppmR,'PQN');
% normcheck(spec_mat_N)
varcheck(spec_mat)
spec_mat_S=scale(spec_mat,'logoff');
varcheck(spec_mat_S)
% PCA=nipalsPCA(spec_mat,5);
PCA=nipalsPCA(spec_mat_S,5);
figure
hold
VisScores(spec_mat_S,PCA,[1 2],'Y',sampleindvec,'conf_ellipse',true,'showlegend',sampleindvec);
%% Loadings plots: PC-1 seems related to the glycine correlation
VisLoadings1D(spec_mat,PCA.loadings(1,:),ppm)
%% PC-2
VisLoadings1D(spec_mat,PCA.loadings(2,:),ppm)
%STOCKSY
STOCSY(3.19498,spec_mat,ppm);%3.2363


%%%%test on the function spec_conc_simu()
close all;
clear all;
comp='/Users/yuewu/';%the computer user location
workdir=[comp 'Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/spec_simu_gissmo/'];%%the working folder
libdir=[comp 'Dropbox (Edison_Lab@UGA)/Resources/gissmo_lib/'];
cd(workdir);
rng(1);%%random seed
args=struct();%% all the arguments
%%%%%%%%%%%%%%USER ARGUMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%concentration control
args.conc_range=[0 1000];%%the whole dynamic range of all compounds
args.totalnoiseregion=[-0.5 -0.2];%% the region of no signal to shift spectra
args.factor_std=0.3;%%the relative variance for compound concentration in samples
args.sigma_noise=15;%white noise on the whole spectra
%%sample control
args.nsample=50;%%number of samples for each group
args.sampleindvec=[repmat(1,1,args.nsample) repmat(2,1,args.nsample)];%%use this to control samples
%%type of simulation: 'random_chemical_sample' OR 'biological_sample'
args.simutype='biological_sample';
args.group_vec={'Light' 'Dark'}; %different biological groups
%%phenotype control
args.pheno_vec={'growth_rate' 'weight' 'offsprings'};
args.phenorang=[0 10];% the range of phenotype value
%% compound control
args.compd_vec={'ATP' 'Acetic-acid' 'Adenosine' 'Betaine' 'Choline' 'Citrate' 'Cytosine' 'DSS' 'Ethanol' 'Formate' 'Fumaric-acid' 'Glycerol' 'Glycine' 'L-Serine' 'L-Threonine' 'L-Valine' 'L-alanine' 'L-arginine' 'L-glutamic-acid' 'L-malic-acid'};
%%concentration that are changed between groups
%%%{compound, relative mean concentration in different samples}
args.modif_conc_list=cell2table({'L-arginine' 1 3; 'Choline' 2 1});
%%correlation between compound and phenotypes
%%%group index, compound, phenotype index, simulated correlation
args.corr_conc_list=cell2table({2 'Glycine' 1 0.8});
%%technical control
args.lambda=-0.00035;% the exponential window function for fid, need to be modified according to inspection of result spectra
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

%% not enough arguemnts: ERROR message
argstemp=rmfield(args,'simutype');
simures=spec_conc_simu(argstemp);
conc_mat=simures.conc_mat;

%%random_chemical_sample
args.simutype="random_chemical_sample";
simures=spec_conc_simu(args);
conc_mat=simures.conc_mat;
heatmap(conc_mat);
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
stackSpectra(spec_mat,ppm,0.001,10,'simulated spectra');
testregion=[1.8 2.0]%[1.6 1.8]%[3.1 3.3]% [1.8 2.0] [3.7 3.8]
ppmrang=matchPPMs(testregion,ppm);
ppmind=ppmrang(1):ppmrang(2);
mattest=spec_mat(:,ppmind);
ppmtest=ppm(ppmind);
samples=1:length(args.group_ident_vec);
%NO group pattern
fig=figure(), hold on
    surf(ppmtest,samples',mattest,'FaceColor','Interp');
    % surf(mattest);
    shading interp;
    set(gca,'xdir','reverse');
    ylabel('sample');
    zlabel('intensity');
    title(['shiftregion']);
    xlabel('ppm');


%no concentration diff between groups: no pattern
argstemp=rmfield(args,'modif_conc_list');
simures=spec_conc_simu(argstemp);
conc_mat=simures.conc_mat;
heatmap(conc_mat);
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
testregion=[1.8 2.0]%[1.6 1.8]%[3.1 3.3]% [1.8 2.0] [3.7 3.8]
ppmrang=matchPPMs(testregion,ppm);
ppmind=ppmrang(1):ppmrang(2);
mattest=spec_mat(:,ppmind);
ppmtest=ppm(ppmind);
samples=1:length(args.group_ident_vec);
%NO group pattern
fig=figure(), hold on
    surf(ppmtest,samples',mattest,'FaceColor','Interp');
    % surf(mattest);
    shading interp;
    set(gca,'xdir','reverse');
    ylabel('sample');
    zlabel('intensity');
    title(['shiftregion']);
    xlabel('ppm');

%%correlation between compound and phenotypes: not clear correlation
argstemp=rmfield(args,'corr_conc_list');
simures=spec_conc_simu(argstemp);
conc_mat=simures.conc_mat;
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
testregion=[3.54 3.55];
ppmrang=matchPPMs(testregion,ppm);
ppmind=ppmrang(1):ppmrang(2);
groupind=find(strcmp(args.group_ident_vec,'Dark'));
mattest=spec_mat(groupind,ppmind);
ppmtest=ppm(ppmind);
corr(max(mattest,[],2),args.pheno_mat(groupind,1))
corr(args.pheno_mat(groupind,1),conc_mat(groupind,13))
