function PLS=plsCV(Y,matrix,nfold,type,permutations)

% PLS=plsCV(Y,matrix,nfold,type,permutations)
%
% Calculates a partial least squares regression of Y on X using PLS components of X by SIMPLS algorithm and
% outputs an array with scores, and loadings for each PLS component and the
% PLS coefficients vector.  This function uses nfold cross-validation to
% avoid overfitting, and the user will be asked to select the optimal
% number of PLS components to avoid over- and under-fitting the data.  This
% function can perform PLS regression or PLS discriminant analysis.
%
% Arguments:
% 
% Y            Response vector taking real values in the case of regression
%              or dummy variables (ie, ones and zeros) in the case of 2 class
%              discriminant analysis 
% matrix       Data matrix from spectra: binmat or full X matrix
% nfold        N-fold cross validation: ie, 5 specifies five-fold cross
%              validation which splits X up into 5 sets of 80% training and 20%
%              validation data.  The model constructed from the training data is
%              evaluated on the validation data and R2 and Q2 values are returned.
% type         're' for regression, or 'da' for discriminant analysis
% permutation  Number of permutations to test cross-validation.  
%
% Return Values:
% PLS.scores      Scores for each spectrum in X along each PLS component
% PLS.loadings    Coefficients of each variable contributing to generation of the scores in each PLS component.
% PLS.beta        PLS coefficients constructed from user selected optimal
%                 number of PLS components
%
% Log: MTJ added parameter reporting 30NOV2020 (see lines 33, 124-126).

p = reportParams('exclude',{'X','ppm'},'sizeLimit',3);

if isempty(ver('stats'))==1
    error('This function requires the Matlab Statistics Toolbox')
end

if size(matrix,3)>1
    for k=1:size(matrix,3)
        plsmat(k,:)=reshape(matrix(:,:,k),1,[]);
    end
    maxsize=20;
else
    plsmat=matrix;
    maxsize=permutations;
end

if size(Y,2)>size(Y,1)
    Y=Y';
end

nlength=round(size(plsmat,1)/nfold);

NFoldCombinations=nchoosek(size(plsmat,1),nlength);
if NFoldCombinations>maxsize
    permutation_set=zeros(maxsize,size(plsmat,1));
    for k=1:maxsize
        permutation_set(k,:)=randperm(size(plsmat,1));
    end
    sets=permutation_set(:,1:nlength);
else
    sets=nchoosek(1:size(plsmat,1),nlength);
end
crossval=size(sets,1);

if (size(plsmat,1)-(nlength+1))>10
    numcomp=10;
else
    numcomp=size(plsmat,1)-(nlength+1);
end

for z=1:numcomp
    for k=1:crossval;
        
        validate=sets(k,:);
        train=setdiff(1:size(plsmat,1),validate);

        [XL,yl,XS,YS,beta,PCTVAR] = plsregress(plsmat(train,:),Y(train),z);
        fit1=[ones(size(plsmat(validate,:),1),1) plsmat(validate,:)]*beta;
        fit2 = [ones(size(plsmat(train,:),1),1) plsmat(train,:)]*beta;
        if type(1)=='r'
            sy=norm(Y-mean(Y))./sqrt(length(Y-mean(Y)));
            yfit_test = fit1;
            yfit_train = fit2;
            q2(k,z)=((sy.^2)-sum((yfit_test-Y(validate)).^2)./length(yfit_test))/(sy.^2);
            r2(k,z)=((sy.^2)-sum((yfit_train-Y(train)).^2)./length(yfit_train))/(sy.^2);
        elseif type(1)=='d'
            yfit_test=fit1;
            yfit_train=fit2;
            
            %from Westerhuis et al Metabolomics 2008
            q2(k,z)=1-(sum((Y(validate)-yfit_test).^2)/sum((Y(validate)-mean(Y)).^2));
            r2(k,z)=1-(sum((Y(train)-yfit_train).^2)/sum((Y(train)-mean(Y)).^2));
             
            %percentage right
            percent_correct_test(k,z)=length(Y(validate))-sum([Y(validate)-round(fit1)].^2)/length(Y(validate));
            percent_correct_train(k,z)=length(Y(train))-sum([Y(train)-round(fit2)].^2)/length(Y(train));
        else
            error('Type must be "re" or "da"')
        end
        
        pervarexp(k,z)=sum(PCTVAR(1,:));
    end
end

figure, plot(mean(q2),'r')
hold on; plot(mean(r2),'g')
hold on; plot(mean(pervarexp))
legend('Q2','R2','Percent Variance in X');
xlabel('Number of PLS components');

numcomponents=input('Select best number of components: ');

[XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(plsmat,Y,numcomponents);

PLS.scores=XS;
PLS.loadings=XL';
PLS.betas=beta';
PLS.q2=mean(q2(:,1:numcomponents));
PLS.r2=mean(r2(:,1:numcomponents));
PLS.variance=PCTVAR(1,1:numcomponents);

%% MTJ added parameter reporting 30NOV2020
    PLS.params = p;
    PLS.params.userInput.components = numcomponents;

end


