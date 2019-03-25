function normcheck(X)

% normcheck(X)
%
% Displays histogram and box plots of log-fold change vs. median for all
% features in each spectrum.  Dilution / normalization effects are often
% visible as distributions not centered at 0.  
%
% Arguments:
% 
% X        N x P matrix of spectral features for each sample
%

X=abs(X);

F=X./repmat(median(X),[size(X,1),1]);
x=-4:.16:4;
n=histc(log(F)',x);
figure, 
hold

subplot(2,1,1)
plot(x,n)
title('Histogram of log-fold-change vs. median values')
xlabel('log-fold-change vs. median  feature value')
ylabel('Features (count)')



    if isempty(ver('stats'))==1
        disp('Boxplots require the Matlab Statistics Toolbox')
    else
subplot(2,1,2)    
    boxplot(log(F)','plotstyle','compact','symbol',' ')
    ylim([-4 4])
    title('Boxplots of log-fold-change vs. median values')
    xlabel('Observation')
    ylabel('log-fold-change vs. median  feature value')
end