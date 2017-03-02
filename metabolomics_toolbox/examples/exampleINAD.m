%% Preprocessing 1Ds

% % class membership 
% Y=[0 0 0 0 1 1 1 1 1 1 1 1]';
% Y2=[1 1 1 1 2 2 2 2 3 3 3 3]';

d=dir;
for i=3:length(dir)
    if strfind(d(i).name,'.ft')
        ftlist{i-2}=[cd '/' d(i).name];
    end
end
ftlist(cellfun(@isempty,ftlist))=[];

for i=1:length(ftlist)
    spectra(i)=pipe2matlab(ftlist{i});
end

%% plot first spectrum in 'fast' mode and 'contour' mode.  
vis2D(spectra(1),'f')
%vis2D(spectra(1),'c',20)

%make 3D matrix of 2D spectra
Setup2D(spectra);

%test basic segmentation with high S/N threshold
% Xg=X(1:8:end,1:8:end,:);
% Xng=XNoise(1:8:end,1:8:end,:);
% ppm1g=ppm1(1:8:end);
% ppm2g=ppm2(1:8:end);

%% Here I reduce the dimensions (chemical shift). 

% % WWater
% X=X(:,1:8:end,:);
% XNoise=XNoise(:,1:8:end,:);
% ppm1=ppm1(1:8:end);
% ppm2=ppm2(1:8:end);
% 
% % WPellet
% X=X(1:6:end,1:3:end,:);
% XNoise=XNoise(1:6:end,1:3:end,:);
% ppm1=ppm1(1:3:end);
% ppm2=ppm2(1:6:end);
%% segment the INADEQUATE
%labels=segment2D(Xg,Xng,ppm1g,ppm2g,1);
%label=segmentINAD(X(:,:,[1:4 6:8]),XNoise(:,:,[1:4 6:8]),ppm1,ppm2);
label=segmentINAD(X,XNoise,ppm1,ppm2);
figure, imagesc(ppm1,ppm2,label)
set(gca,'XDir','reverse')

%%

%load label matrix calculated with default settings
% load('../data/label.mat')
% figure, imagesc(ppm1,ppm2,label)
% set(gca,'XDir','reverse')

%Check if alignment is necessary
%stackplot(X,XNoise,ppm1,ppm2,Y)
% 
%if so, use either guided alignment (HATS) or star alignment. For this
%data, alignment is not necessary.  
%XAL=star_align2D(X,label,ppm1,ppm2);
%XALg=HATS(X(:,1:8:end,[1:4 6:8]),XNoise,ppm1(1:8:end),ppm2,label);
XAL=HATS(X,XNoise,ppm1,ppm2,label);

%Bin 2D spectra using segmentation in label matrix
%binmat=bin2D(XAL(:,:,[1:2 4:5 7:8]),label);
binmat=bin2D(XAL,label);

%% Check Scaling and Normalization
%save WP_INAD_referenced_noreduced.mat -v7.3

%check distribution of fold-changes before and after normalization for full
%resolution data
% normcheck(binmat)
binmatN=normalize(binmat,'PQN');
% normcheck(binmatN)

%check distribution of variance before and after stabilization for full
%resolution data
% varcheck(binmatN)
binmatSN=scale(binmatN,'logoff');
% varcheck(binmatSN)


%% Analysis 

%HCA analysis of INAD features
%[sample_order,variable_order]=two_way_cluster(binmatSN,'weighted','spearman',Y);

%% PCA analysis of INAD features
PCA=nipalsPCA(binmatSN,5);
VisScores(binmatSN,PCA,[1 2],Y);
%figure, imagesc(PCA.loadings(:,variable_order))
% caxis([-1*max(max(abs(PCA.loadings))) 1*max(max(abs(PCA.loadings)))])
VisLoadings2D(XAL,PCA.loadings(1,:),ppm1,ppm2,'c',label);

%% PLS regression of INAD features
PLS=plsPV(binmatN,Y,4,'da',10,'logoff');
VisScores(binmatN,PLS,[1 2],Y);
predictors=PLS.betas(2:end);
figure, imagesc(predictors(variable_order))
caxis([-1*max(abs(predictors)) 1*max(abs(predictors))])
VisLoadings2D(X,predictors,ppm1,ppm2,'c',label);

%% t-test of INAD features
[p,sigs]=MWAS(binmatSN,Y,'bonferroni');
manhattan(binmatSN,Y,1:size(binmat,2),p,sigs,'bonferroni')
displaypeak2D(X,XNoise,label,find(sigs==1),ppm1,ppm2)

