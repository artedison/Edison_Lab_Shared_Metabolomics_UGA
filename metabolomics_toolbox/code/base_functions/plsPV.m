function PLS=plsPV(X,Y,nfold,type,permutations,method)

% PLS=plsPV(X,Y,nfold,type,permutations,method)
%
% Calculates a partial least squares regression of Y on X using PLS components of X by SIMPLS algorithm and
% outputs an array with scores, and loadings for each PLS component and the
% PLS coefficients vector.  This function uses nfold permutation-validation to
% avoid overfitting, and the user will be asked to select the optimal
% number of PLS components to avoid over- and under-fitting the data.  This
% function can perform PLS regression or PLS discriminant analysis.
% Scaling is performed during validation in order to avoid the introduction
% of bias.
%
% Arguments:
%
% X            Data matrix from spectra: binmat or full X matrix
% Y            Response vector taking real values in the case of regression
%              or dummy variables (ie, ones and zeros) in the case of 2 class
%              discriminant analysis
% nfold        N-fold cross validation: ie, 5 specifies five-fold cross
%              validation which splits X up into 5 sets of 80% training and 20%
%              validation data.  The model constructed from the training data is
%              evaluated on the validation data and R2 and Q2 values are returned.
% type         're' for regression, or 'da' for discriminant analysis
% permutation  Number of permutations to test cross-validation.
% method       (Optional) - 'log' for log2 fold-change vs. median scaling, 'logoff'
%              for offset log2 fold-change vs. median,
%              'mc' for mean-centering (only), 'auto' for autoscaling
%              (mean-center and univariance scaling), 'pareto' for
%              pareto scaling. 'integral' for normalization to sum of
%              set of features
% Return Values:
% PLS.scores      Scores for each spectrum in X along each PLS component
% PLS.loadings    Coefficients of each variable contributing to generation of the scores in each PLS component.
% PLS.beta        PLS coefficients constructed from user selected optimal
%                 number of PLS components
% PLS.r2          Fit of model to training data
% PLS.q2          Fit of model to test data
% PLS.variance    Variance in X for each component

% Edits:
%       MJ 14NOV2017 legend call had extra parameter "1". Removed. 

if isempty(ver('stats'))==1
    error('This function requires the Matlab Statistics Toolbox')
end

if size(X,3)>1
    for k=1:size(matrix,3)
        matrix(k,:)=reshape(X(:,:,k),1,[]);
    end
    maxsize=20;
else
    matrix=X;
    maxsize=permutations;
end

if size(Y,2)>size(Y,1)
    Y=Y';
end

nlength=round(size(matrix,1)/nfold);

NFoldCombinations=nchoosek(size(matrix,1),nlength);
if NFoldCombinations>maxsize
    permutation_set=zeros(maxsize,size(matrix,1));
    for k=1:maxsize
        permutation_set(k,:)=randperm(size(matrix,1));
    end
    sets=permutation_set(:,1:nlength);
else
    sets=nchoosek(1:size(matrix,1),nlength);
end
crossval=size(sets,1);

if (size(matrix,1)-(nlength+1))>10
    numcomp=10;
else
    numcomp=size(matrix,1)-(nlength+1);
end

for z=1:numcomp
    for k=1:crossval;
        
        test=sets(k,:);
        train=setdiff(1:size(matrix,1),test);
        
        if exist('method')~=1
            plsmat=matrix;
        else
            switch method
                case 'log'
                    matrix=abs(matrix);
                    if exist('offset')~=1
                        offset=0.0001;
                    end
                    median_X=median(matrix(train,:)+offset);
                    FC=zeros(size(matrix));
                    for index=1:size(matrix,1)
                        FC(index,:)=[matrix(index,:)+offset]./median_X;
                    end
                    plsmat=log2(FC);
                case 'logoff'
                    matrix=abs(matrix);
                    if exist('offset')~=1
                        trainset=matrix(train,:);
                        offset=median(trainset(trainset>0));
                    end
                    median_X=median(matrix(train,:)+offset);
                    FC=zeros(size(matrix));
                    for index=1:size(matrix,1)
                        FC(index,:)=[matrix(index,:)+offset]./median_X;
                    end
                    plsmat=log2(FC);
                case 'mc'
                    plsmat=matrix-repmat(mean(matrix(train,:)),[size(matrix,1),1]);
                case 'auto'
                    XMC=matrix-repmat(mean(matrix(train,:)),[size(matrix,1),1]);
                    plsmat=XMC./repmat(std(matrix(train,:)),[size(matrix,1),1]);
                case 'pareto'
                    XMC=matrix-repmat(mean(matrix(train,:)),[size(matrix,1),1]);
                    plsmat=XMC./repmat(sqrt(std(matrix(train,:))),[size(matrix,1),1]);
                    
            end
        end
        
        [XL{z}(k,:,:),yl,XS{z}(k,:,:),YS,beta{z}(k,:,:),PCTVAR] = plsregress(plsmat(train,:),Y(train),z);
        fit1=[ones(size(plsmat(test,:),1),1) plsmat(test,:)]*beta{z}(k,:)';
        fit2 = [ones(size(plsmat(train,:),1),1) plsmat(train,:)]*beta{z}(k,:)';
        
        
        if type(1)=='r'
            sy=norm(Y-mean(Y))./sqrt(length(Y-mean(Y)));
            yfit_test = fit1;
            yfit_train = fit2;
            q2(k,z)=((sy.^2)-sum((yfit_test-Y(test)).^2)./length(yfit_test))/(sy.^2);
            r2(k,z)=((sy.^2)-sum((yfit_train-Y(train)).^2)./length(yfit_train))/(sy.^2);
        elseif type(1)=='d'
            yfit_test=fit1;
            yfit_train=fit2;
            
            %from Westerhuis et al Metabolomics 2008
            q2(k,z)=1-(sum((Y(test)-yfit_test).^2)/sum((Y(test)-mean(Y)).^2));
            r2(k,z)=1-(sum((Y(train)-yfit_train).^2)/sum((Y(train)-mean(Y)).^2));
            
            %percentage right
            percent_correct_test(k,z)=length(Y(test))-sum([Y(test)-round(fit1)].^2)/length(Y(test));
            percent_correct_train(k,z)=length(Y(train))-sum([Y(train)-round(fit2)].^2)/length(Y(train));
        else
            error('Type must be "re" or "da"')
        end
        
        pervarexpX(k,z)=PCTVAR(1,z);
        pervarexpY(k,z)=PCTVAR(2,z);
    end
end

figure, plot(mean(q2),'r')
hold on; plot(mean(r2),'g')
hold on; plot(cumsum(mean(pervarexpX)),'m')
hold on; plot(cumsum(mean(pervarexpY)),'k')
legend('Q2','R2','Percent Variance in X','Percent Variance in Y');
%legend('Q2','R2','Percent Variance in X','Percent Variance in Y',1);
xlabel('Number of PLS components');

numcomponents=input('Select best number of components');

PLS.loadings=(squeeze(mean(XL{numcomponents}))');
PLS.scores=plsmat*PLS.loadings';
PLS.betas=mean(beta{numcomponents});
PLS.q2=mean(q2(:,1:numcomponents));
PLS.r2=mean(r2(:,1:numcomponents));
PLS.variance=mean(pervarexpX(:,1:numcomponents));




