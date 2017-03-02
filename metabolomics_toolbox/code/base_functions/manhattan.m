function manhattan(X,Y,features,p,sigs,correction)

% manhattan(X,Y,features,p,sigs,correction)
% 
% Manhattan plot of p values calculated by MWAS.m for all features, and
% boxplots of all significant features by class.  
%
% Arguments:
% 
% X                Data matrix of spectra 
% Y                Vector of classes
% features         Vector of features (chemical shifts or m/z)
% p                P values
% sigs             Optional: logical vector of significant features
% correction       Optional: Either Bonferroni, 'bonferroni', or Sidak,
%                  'sidak'.  Default is 'bonferroni'

if exist('correction')~=1
    beta=0.05/size(X,2);
elseif strcmp(correction,'bonferroni')==1
    beta=0.05/size(X,2);
elseif strcmp(correction,'sidak')==1
    beta=1-((1-0.05)^(1/size(X,2)));
else
    Error('Correction must be either bonferroni or sidak')
end

if exist('sigs')~=1
    sigs=zeros(1,size(X,2));
    sigs(p<beta)=1;
end
    
[p_sort,p_order]=sort(log10(p));

%Manhattan plot
significant=find(sigs==1);
figure, n=stem(features,abs(log10(p)),'k');
hold on;plot(features,-log10(beta)*ones(1,length(sigs)),'--b','LineWidth',1.5)
hold on;m=stem(features(significant),abs(log10(p(significant))),'fill','k');
set(m,'MarkerFaceColor','red')
ylim([.1*min(log10(p)),-1.3*min(log10(p))]);
ylabel('-Log(p)')

%Box plots
if isempty(ver('stats'))==1
    disp('Boxplots require the Matlab Statistics Toolbox')
else
    if length(significant>0) && length(unique(Y))<=10
        figure
        subplot_size=ceil(sqrt(length(significant)));
        for k=1:length(significant)
            subplot(subplot_size,subplot_size,k)
            boxplot(X(:,significant(k)),Y);
            ylabel(['Feature at ',num2str(features(significant(k)))])
            title(['P value = ',num2str(p(significant(k)))])
        end
    end
end