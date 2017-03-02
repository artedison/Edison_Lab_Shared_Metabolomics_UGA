%%
loadStudyFTdata('ft_all.csv')

%%
loadStudyFTdata('ft_edit.csv')

%%
[X,ppm,Xtitles]=Setup1D(spectra);

%%
spectra = ref_spectra_manual(spectra, [3.02,3.14], 3.0204)

%%
[X,ppm,Xtitles]=Setup1D(spectra);

%% 706=1,729=1,734=2,770=2

Y=[1 1 1 1 1 1 2 2 2 2 2 2];

%% RV=1,LV=2,S=3

Y2=[1 2 3 3 1 2 1 3 2 2 1 3];

%% 706=1,729=1,734=2,770=2

Y=[1 1 1 1 2 2 2 2 2 2];

%% RV=1,LV=2,S=3

Y2=[1 2 3 2 1 3 2 2 1 3];


    
%%
plotr(ppm,X)

%%
XR=remove_region(X,ppm,4.793,5.144);

%%
plotr(ppm,XR)

%%
[XR1,ppmR]=remove_ends(XR,ppm,0.5,9);

%%
plotr(ppmR,XR1)

%% baseline correction
[A,S]=showBaseline(XR1, ppmR, 9e6, 100);

%% baseline correction
XNew=CorrectBl(XR1, A, S);           

%%
plotr(ppmR,XNew)

%%
XAL=star_align1D(XNew,ppmR,'mean','CCOW');

%%
XALg=guide_align1D(XNew,ppmR,'correlation','CCOW');


%%
plotr(ppmR,XAL)

%%
figure, plotr(ppmR,XALg)

%%
normcheck(XAL)

%%
W=[17.2 24.7 24.4 22.4 25.6 27.6 23.1 25.4 25.7 27.3 24.4 27.7];
W=W'

%% normalize to weight
XN1=W(1,:).*XAL(1,:);
XN2=W(2,:).*XAL(2,:);
XN3=W(3,:).*XAL(3,:);
XN4=W(4,:).*XAL(4,:);
XN5=W(5,:).*XAL(5,:);
XN6=W(6,:).*XAL(6,:);
XN7=W(7,:).*XAL(7,:);
XN8=W(8,:).*XAL(8,:);
XN9=W(9,:).*XAL(9,:);
XN10=W(10,:).*XAL(10,:);
XN11=W(11,:).*XAL(11,:);
XN12=W(12,:).*XAL(12,:);
XALN=vertcat(XN1,XN2,XN3,XN4,XN5,XN6,XN7,XN8,XN9,XN10,XN11,XN12);
normcheck(XALN)

%% normalize script 

for i=1:length(Y)
    XALN(i,:)=XAL(i,:)/W(i);
end

%%

normcheck(XAL)

%% don't use
XALN=normalize(XAL,ppmR,'PQN');

%% 
figure,plotr(ppmR,XALN(Y==1,:),'b')
hold
plotr(ppmR,XALN(Y==2,:),'r')

%% 
plotr(ppmR,XALN(Y2==1))

%%
XALNS=scale(XALN,'pareto');

%%
PCA=nipalsPCA(XALNS,5);

%%
VisScores(XALNS,PCA,[1 2],'Y',Y,'showlabels',true,'mahalanobis',[1 2]); 
VisScores(XALNS,PCA,[1 3],'Y',Y,'showlabels',true,'mahalanobis',[1 2]);  
VisScores(XALNS,PCA,[1 4],'Y',Y,'showlabels',true,'mahalanobis',[1 2]); 
figure, VisScores(XALNS,PCA,[1 5],'Y',Y,'showlabels',true,'mahalanobis',[1 2]);  
figure, VisScores(XALNS,PCA,[2 4],'Y',Y,'showlabels',true,'mahalanobis',[1 2]); 
VisScores(XALNS,PCA,[2 5],'Y',Y,'showlabels',true,'mahalanobis',[1 2]); 
VisScores(XALNS,PCA,[3 4],'Y',Y,'showlabels',true,'mahalanobis',[1 2]); 
VisScores(XALNS,PCA,[3 5],'Y',Y,'showlabels',true,'mahalanobis',[1 2]); 
VisScores(XALNS,PCA,[4 5],'Y',Y,'showlabels',true,'mahalanobis',[1 2]); 

%%

VisScores(XALNS,PCA,[1 2],'Y',Y2,'showlabels',true); 
VisScores(XALNS,PCA,[1 3],'Y',Y2,'showlabels',true);  
VisScores(XALNS,PCA,[1 4],'Y',Y2,'showlabels',true); 
VisScores(XALNS,PCA,[1 5],'Y',Y2,'showlabels',true);  
VisScores(XALNS,PCA,[2 3],'Y',Y2,'showlabels',true); 
VisScores(XALNS,PCA,[2 4],'Y',Y2,'showlabels',true); 
VisScores(XALNS,PCA,[2 5],'Y',Y2,'showlabels',true); 
VisScores(XALNS,PCA,[3 4],'Y',Y2,'showlabels',true); 
VisScores(XALNS,PCA,[3 5],'Y',Y2,'showlabels',true); 
VisScores(XALNS,PCA,[4 5],'Y',Y2,'showlabels',true); 

%% pls
PLS=plsPV(XALN,Y,7,'da',10,'pareto');

%% pls
PLS=plsCV(Y,XALNS,5,'da',10)

%%
VisScores(XALNS,PLS,[1 2],'Y',Y,'showlabels',true,'mahalanobis',[1 2]); 
VisScores(XALNS,PLS,[1 3],'Y',Y,'showlabels',true,'mahalanobis',[1 2]);  
VisScores(XALNS,PLS,[1 4],'Y',Y,'showlabels',true,'mahalanobis',[1 2]); 
VisScores(XALNS,PLS,[1 5],'Y',Y,'showlabels',true,'mahalanobis',[1 2]);  
VisScores(XALNS,PLS,[2 3],'Y',Y,'showlabels',true,'mahalanobis',[1 2]); 
VisScores(XALNS,PLS,[2 4],'Y',Y,'showlabels',true,'mahalanobis',[1 2]); 
VisScores(XALNS,PLS,[2 5],'Y',Y,'showlabels',true,'mahalanobis',[1 2]); 
VisScores(XALNS,PLS,[3 4],'Y',Y,'showlabels',true,'mahalanobis',[1 2]); 
VisScores(XALNS,PLS,[3 5],'Y',Y,'showlabels',true,'mahalanobis',[1 2]); 
VisScores(XALNS,PLS,[4 5],'Y',Y,'showlabels',true,'mahalanobis',[1 2]); 

%% pls
PLS2=plsPV(XALNS,Y2,5,'da',10,'pareto');

%%
VisScores(XALNS,PLS2,[1 2],'Y',Y2);
VisScores(XALNS,PLS2,[1 3],'Y',Y2);
VisScores(XALNS,PLS2,[1 4],'Y',Y2);
VisScores(XALNS,PLS2,[1 5],'Y',Y2);
VisScores(XALNS,PLS2,[2 3],'Y',Y2);
VisScores(XALNS,PLS2,[2 4],'Y',Y2);
VisScores(XALNS,PLS2,[2 5],'Y',Y2);
VisScores(XALNS,PLS2,[3 4],'Y',Y2);
VisScores(XALNS,PLS2,[3 5],'Y',Y2);
VisScores(XALNS,PLS2,[4 5],'Y',Y2);

%%
VisLoadings1D(XALN,PCA.loadings(1,:),ppmR)

%%
STOCSY
%%
[hValine,pValine] = ttest_roi(XALN,ppmR,Y,1.149,1.101,[1 2],'bonferroni')
[halanine,palanine] = ttest_roi(XALN,ppmR,Y,1.488,1.44,[1 2],'bonferroni')
[hCreatine,pCreatine] = ttest_roi(XALN,ppmR,Y,3.049,2.985,[1 2],'bonferroni')
[hcholines,pcholines] = ttest_roi(XALN,ppmR,Y,3.297,3.159,[1 2],'bonferroni')
[hPEA,pPEA] = ttest_roi(XALN,ppmR,Y,3.572,3.495,[1 2],'bonferroni')
[hTaurine1,pTaurine2] = ttest_roi(XALN,ppmR,Y,3.437,3.376,[1 2],'bonferroni')
[hTaurine2,pTaurine2] = ttest_roi(XALN,ppmR,Y,3.635,3.574,[1 2],'bonferroni')
[hCreatine2,pCreatine2] = ttest_roi(XALN,ppmR,Y,3.933,3.911,[1 2],'bonferroni')
[hSugar1,pSugar1] = ttest_roi(XALN,ppmR,Y,4.071,4.002,[1 2],'bonferroni')
[hSugar2,pSugar2] = ttest_roi(XALN,ppmR,Y,4.139,4.073,[1 2],'bonferroni')
[hLactate,pLactate] = ttest_roi(XALN,ppmR,Y,4.136,4.081,[1 2],'bonferroni')
[hAromatic1,pAromatic1] = ttest_roi(XALN,ppmR,Y,8.254,8.203,[1 2],'bonferroni')
[hAromatic2,pAromatic2] = ttest_roi(XALN,ppmR,Y,8.355,8.324,[1 2],'bonferroni')


%%

[hlipidA,plipidA] = ttest_roi(XALN,ppmR,Y,0.9708,0.7638,[1 2],'bonferroni')
[hlipidBC,plipidBC] = ttest_roi(XALN,ppmR,Y,1.44,1.203,[1 2],'bonferroni')
[hlipidD,plipidD] = ttest_roi(XALN,ppmR,Y,1.706,1.505,[1 2],'bonferroni')
[hlipidE,plipidE] = ttest_roi(XALN,ppmR,Y,2.091,1.923,[1 2],'bonferroni')
[hlipidF,plipidF] = ttest_roi(XALN,ppmR,Y,2.308,2.174,[1 2],'bonferroni')
[hlipidG,plipidG] = ttest_roi(XALN,ppmR,Y,2.705,2.687,[1 2],'bonferroni')
[hlipidK,plipidK] = ttest_roi(XALN,ppmR,Y,5.257,5.183,[1 2],'bonferroni')
[hlipidLM,plipidLM] = ttest_roi(XALN,ppmR,Y,5.504,5.257,[1 2],'bonferroni')


%% unsaturation (p=0.2220)
LM=IntegralPeak_roi(XALN,ppmR,5.504,5.257);
A=IntegralPeak_roi(XALN,ppmR,0.9708,0.7638);
unsaturation=LM./((2/3).*A)
[h,p]=ttest2(unsaturation(Y==1),unsaturation(Y==2))

%% polysaturation (p=0.0018)

G=IntegralPeak_roi(XALN,ppmR,2.705,2.687);
A=IntegralPeak_roi(XALN,ppmR,0.9708,0.7638);
polysaturation=G./((2/3).*A)
[h,p]=ttest2(polysaturation(Y==1),polysaturation(Y==2))
PS1=mean(polysaturation(Y==1));
PS2=mean(polysaturation(Y==2));
PS=horzcat(PS1,PS2);
bar(PS)

%% unsaturated FA (p=0.2534)

E=IntegralPeak_roi(XALN,ppmR,2.091,1.923);
F=IntegralPeak_roi(XALN,ppmR,2.308,2.174);
FU=E./(2.*F)
[h,p]=ttest2(FU(Y==1),FU(Y==2))
FU1=mean(FU(Y==1));
FU2=mean(FU(Y==2));
%FUm=horzcat(FU1,FU2);
%bar(FUm)

%% diunsaturated FA (p=0.0082)

G=IntegralPeak_roi(XALN,ppmR,2.705,2.687)
F=IntegralPeak_roi(XALN,ppmR,2.308,2.174)
FD=G./F
[h,p]=ttest2(FD(Y==1),FD(Y==2))
FD1=mean(FD(Y==1));
FD2=mean(FD(Y==2));


%% mcl
A=IntegralPeak_roi(XALN,ppmR,0.9708,0.7638);
BC=IntegralPeak_roi(XALN,ppmR,1.44,1.203);
D=IntegralPeak_roi(XALN,ppmR,1.706,1.505);
E=IntegralPeak_roi(XALN,ppmR,2.091,1.923);
F=IntegralPeak_roi(XALN,ppmR,2.308,2.174);
G=IntegralPeak_roi(XALN,ppmR,2.705,2.687);
LM=IntegralPeak_roi(XALN,ppmR,5.504,5.257);
mcl = [((2/3).*A)+BC+D+E+F+G+(2.*LM)./((2/3).*A)]


