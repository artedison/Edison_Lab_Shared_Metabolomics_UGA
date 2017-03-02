% set working directory to 'code' 
importXCMS('../data/Elegans_report_data_XCMS.xls')
reorder=[9 10 11 12 4 1 2 3 5 6 7 8];
msmatrix=msmatrix(reorder,:);

% class membership 
Y=[0 0 0 0 1 1 1 1 1 1 1 1]';
Y2=[1 1 1 1 2 2 2 2 3 3 3 3]';

%check distribution of fold-changes before and after normalization 
normcheck(intensity)
MNmatrix=normalize(intensity,'PQN');
normcheck(MNmatrix)

%check distribution of variance before and after stabilization
varcheck(MNmatrix)
MSNmatrix=scale(MNmatrix,'logoff');
varcheck(MSNmatrix)

%HCA analysis of MS features
[sample_order,variable_order]=two_way_cluster(MSNmatrix,'weighted','euclidean',Y);

%PCA analysis of MS features
PCA=nipalsPCA(MSNmatrix,10);
figure, VisScores(MSNmatrix,PCA,[1 2],Y);
figure, VisLoadingsMass(MSNmatrix,PCA.loadings(1,:),CMZ);
figure, imagesc(PCA.loadings(:,variable_order))
caxis([-1*max(max(abs(PCA.loadings))) 1*max(max(abs(PCA.loadings)))])

%PLS regression of MS features
PLS=plsPV(MNmatrix,Y,4,'da',10,'logoff');
VisScores(MSNmatrix,PLS,[1 2],Y);
predictors=PLS.betas(2:end);
figure, imagesc(predictors(variable_order))
caxis([-1*max(abs(predictors)) 1*max(abs(predictors))])

%t-test of MS features
[p,sigs]=MWAS(MSNmatrix,Y,'bonferroni');
manhattan(MSNmatrix,Y,CMZ,p,sigs,'bonferroni')
