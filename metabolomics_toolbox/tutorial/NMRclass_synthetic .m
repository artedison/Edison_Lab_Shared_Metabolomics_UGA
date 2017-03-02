%NMR class: Synthetic Data

%% Load spectra into workspace

loadallft

clear d ftlist i 

figure, plot(spectra(1).ppm,spectra(1).real)
set(gca,'XDir','rev')

%% Reference to Internal Standard
[spectra]=refspec(spectra,0.01);

figure, plot(spectra(1).ppm,spectra(1).real)
set(gca,'XDir','rev')
%% Set up your X matrix and visualize the spectra to make sure they look
%% correct



[X,ppm,XTitles]=Setup1D(spectra);


figure, plot(ppm,X)
set(gca,'XDir','rev')

%Notice that although there was water supression but the residual water
%should still be taken out sot hat it is not taken into account for our
%anaysis. 

%% Make your identity vector

y=[0 0 0 0 0 1 1 1 1 1];

%% Spectra colored by group


figure
hold
plot(ppm,X(y==0,:),'b')
plot(ppm,X(y==1,:),'r')
set(gca,'XDir','reverse')


%% Take out the residual water


XR=remove_region(X,ppm,4.85,5.1);

% and plot so that you can be sure that the correct region was removed

figure, plot(ppm,XR)
set(gca,'XDir','rev')

% Be sure to take a look at the spectra and see how the peaks are
% overlayed. 

%% Alignment of peaks

XAL=star_align1D(XR,ppm,'mean','PAFFT');

% and plot
figure, plot(ppm,XAL)
set(gca,'XDir','reverse')

%% Check normalization 

normcheck(XAL)

%reference:
% Craig, A., et al., Scaling and normalization effects in NMR spectroscopic metabonomic 
%data sets. Anal Chem, 2006. 78(7): p. 2262-2267.


%% Normalize - there are many different ways to normalize your data. Need
%% to try them out to see which one works best for your data

%PQN is a method is based on the calculation of a most probable dilution factor 
%by looking at the distribution of the quotients of the amplitudes of a test spectrum by 
%those of a reference spectrum

% Dieterle, F., et al., Probabilistic quotient normalization as robust method to account for 
%dilution of complex biological mixtures. Application in H-1 NMR metabonomics.
%Anal Chem, 2006. 78(13): p. 4281-4290.


XALN=normalize(XAL,'PQN');
normcheck(XALN)

% and visualize

figure, plot(ppm,XALN)
set(gca,'XDir','rev')

%% If you have an internal standard that you would like to normalize to.. for
%% example TSP is a reference standard that should be at 0 ppm. 

[~,index]=min(abs(ppm-0.0));
XALN=normalize(XAL,'intensity',index);

figure, plot(ppm,XALN)
set(gca,'XDir','rev')
normcheck(XALN)


%% Check the scaling

varcheck(XALN)

%% Try different scaling methods 

%logoff
XALSN=scale(XALN,'logoff');
varcheck(XALSN)


%% pareto 

XALSN=scale(XALN,'pareto');
varcheck(XALSN)

%% Principle Component Analysis- lets look at 5 components


PCA=nipalsPCA(XALSN,5);



%% Scaled


VisScores(XALSN,PCA,[1 2],y);
VisScores(XALSN,PCA,[1 3],y);
VisScores(XALSN,PCA,[1 4],y);
VisScores(XALSN,PCA,[1 5],y);
VisScores(XALSN,PCA,[2 3],y);
VisScores(XALSN,PCA,[2 4],y);
VisScores(XALSN,PCA,[2 5],y);
VisScores(XALSN,PCA,[3 4],y);
VisScores(XALSN,PCA,[3 5],y);
VisScores(XALSN,PCA,[4 5],y);

%% Visualize the loadings plot associated with the scores plot that
%% separates best

VisLoadings1D(XALSN,PCA.loadings(1,:),ppm)

%% STOCSY with the peak that separates the two groups from your loadings
%% plot

STOCSY(30.95,XALN,ppm);
