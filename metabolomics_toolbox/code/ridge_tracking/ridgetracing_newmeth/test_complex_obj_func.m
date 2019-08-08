close all;
clear all;
global dirp;
dirp='/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/';
workdir=[dirp 'Bioinformatics_modeling/spectral.related/result/simulate.pH/'];
% addpath([dirp 'Bioinformatics_modeling/spectral.related/code/simulatespectra/'])
cd(workdir);
%%multiple dataset comparison
%%%%HRMAS experiment data 3 aerobic + 3 anaerobic
load(['./experiment.data.nc/sampleData.mat']);
load([workdir 'experiment.data.nc/tracing.newmeth.experiment.manual.5.mat']);
samp_i=6;%1:6
data=sampleData(samp_i);
mat=data.Xcollapsed_1h1d;
ppm=data.ppm_1h1d;
time=data.timesCollapsed_1h1d;
objstr=complex_obj_func(mat,ppm,time',Sample(samp_i),[-0.4 -0.2],10,[-0.4 10]);
%%%%simulated data
load('simulated.timeseries.ph3.mat');
load([workdir 'manual.res/tracing.newmeth.simulated.manual.mat'])
samp_i=3;%1:3
mat=matrixstr{samp_i};
objstr=complex_obj_func(mat,ppm,times',Sample(samp_i),[-0.4 -0.2],10,[-0.4 10])

%%%%HMRAS experiment long time
workdir2=[dirp 'Bioinformatics_modeling/spectral.related/result/newlongdata/'];
load([workdir2 '/HRMAS_ncrassa_conidia_1c.mat']);
load([workdir2 'tracing.newmeth.experiment.manual.newlong.full.mat'])
samp_i=1;
mat=X_1h1d;
ppm=ppm_1h1d;
time=times_1h1d;
objstr=complex_obj_func(mat,ppm,time',Sample(samp_i),[-0.4 -0.2],10,[-0.4 10])
