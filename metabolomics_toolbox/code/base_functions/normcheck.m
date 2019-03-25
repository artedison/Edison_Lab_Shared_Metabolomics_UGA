function normcheck(X)

% Author: Edison Lab
    % Ver 0.1
    % Tested on Matlab Version R2016b, R2018b
    % Date: 25MAR2019
    %
    %
    % Description:
    %   Displays histogram and box plots of log-fold change vs. median for all
    %   features in each spectrum.  Dilution / normalization effects are often
    %   visible as distributions not centered at 0.  
    % Input:
    %   X: N x P matrix of spectral features for each sample 
    %
    % Output:
    % Histogram and box plots of log-fold change vs. median for all
    % features in each spectrum.
    %
    % Log:
    %   Edited by : SZ, OS
    %   Date      : 25MAR2019
    %   Ver       : 0.1
    %      

X=abs(X);

F=X./repmat(median(X),[size(X,1),1]);
x=-4:.16:4;
n=histc(log(F)',x);
figure, plot(x,n)
title('Histogram of log-fold-change vs. median values')
xlabel('log-fold-change vs. median  feature value')
ylabel('Features (count)')

if isempty(ver('stats'))==1
    disp('Boxplots require the Matlab Statistics Toolbox')
else
    figure, boxplot(log(F)','plotstyle','compact','symbol',' ')
    ylim([-4 4])
    title('Boxplots of log-fold-change and median values')
    xlabel('Observation')
    ylabel('log-fold-change and median  feature value')
end