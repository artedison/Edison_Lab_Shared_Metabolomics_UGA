%add this to startup.m
addpath(genpath('C:\Users\Camera\Documents\edisongit'))

%% Preprocessing 1Ds
loadallft

%% class membership 
Y = [0 0 0 1 1 1];
colors = YtoRGB(Y);
%% make 2D matrix of 1D spectra
[H,ppm,HTitles]=Setup1D(spectra);
figure, [hax,h]=plotr(ppm,HR);
set(h, {'Color'},num2cell(colors,2)) % color by group
%% remove unnecessary regions, check alignment
HR=remove_region(H,ppm,3.2, 3.38);
HR=remove_region(HR,ppm,-Inf,0);
HR=remove_region(HR,ppm,9.0,Inf);
figure, [h1,h]=plotr(ppm,HR);
title('removed regions')
%% alignment
HRA=guide_align1D(HR,ppm,'spearman','PAFFT');
figure, [h2,h]=plotr(ppm,HRA);
title('guided alignment')
linkaxes([h1,h2])
%%
HRA=star_align1D(X,ppm,'median','PAFFT',0.08,.5);
figure, plotr(ppm,HRA)
title('star alignment')
XAL=star_align1D(XAL,ppm,'mean','PAFFT',0.08,.5);
figure, plotr(ppm,XAL)
title('star alignment')

%% peak picking
HRAB=HRA;
HRAB(HRAB<1e4)=0;
[peakmatrix,shifts]=Peakpick1D(HRAB,ppm,'median',.83,'Simple');
%%
[H,~,I_b,~]=opt_bucket(ppm,HRA,.05,.5);
%% or
peaks=peakpick(XAL(1,:),30,0.009*max(XAL(1,:))); %0.001 is the threshold
peakmatrix=XAL(:,peaks);
shifts=ppm(peaks);


%% Check Scaling and Normalization

%check distribution of fold-changes before and after normalization for full
%resolution data
%[~,c]=min(abs(ppm-0));
normcheck(HRA)
HALN=normalize(HRA,ppm,'pqn');
figure, plotr(ppm,HALN)
normcheck(HALN)

%check distribution of variance before and after stabilization for full
%resolution data
varcheck(HALN)
HALSN=scale(HALN,'logoff'); %pwr pareto logoff for H %pareto  for X
HALSN(HALSN==0 | isnan(HALSN) | isinf(HALSN))=0.001;
varcheck(HALN)

%check distribution of fold-changes before and after normalization for
%peaks
normcheck(peakmatrix)
peakmatN=normalize(peakmatrix,'PQN');
normcheck(peakmatN)

%check distribution of variance before and after stabilization for peaks
varcheck(peakmatN)
peakmatSN=scale(peakmatN,'logoff');
varcheck(peakmatSN)

%% Analysis 

%HCA analysis of peaks
[sample_order,variable_order]=two_way_cluster(peakmatSN,'weighted','spearman',Y2);
[sample_order,variable_order]=two_way_cluster(peakmatSN,'weighted','euclidean',Y2);
displaypeak1D(HALNd,ppm,shifts([133]),Y2)
displaypeak1D(HALNd,ppm,shifts(variable_order(308:317)),Y2)

%PCA analysis of peaks
PCA_peaks=nipalsPCA(peakmatSN,5);
VisScores(peakmatSN,PCA_peaks,[1 2],Y);
VisLoadings1Dpeak(peakmatSN,PCA_peaks.loadings(1,:),shifts)
figure, imagesc(PCA_peaks.loadings(:,variable_order))

%PLS regression of peaks
PLS_peaks=plsPV(peakmatN,Y,4,'da',10,'logoff');
VisScores(peakmatN,PLS_peaks,[1 2],Y);
predictors=PLS_peaks.betas(1:end);
figure, imagesc(predictors(variable_order))
VisLoadings1Dpeak(peakmatSN,predictors,shifts)

%t-tcest of peaks
[p,sigs]=MWAS(peakmatSN,Y,'bonferroni');
manhattan(peakmatSN,Y,shifts,p,sigs,'bonferroni')
displaypeak1D(HALNd,ppm,shifts(sigs==1),Y2)

%PCA analysis of full-resolution data
PCA_fullH=nipalsPCA(HALN,5);
VisScores(HALN,PCA_fullH,[1 2],'y',Y);
VisLoadings1D(HALN,PCA_fullH.loadings(1,:),ppm)

%PLS regression of full-resolution data
PLS_full=plsPV(HALSN,Y,4,'da',10,'logoff');
VisScores(HALSN,PLS_full,[1 2],Y);
VisLoadings1D(HALSN,PLS_full.loadings(1,:),ppm)
