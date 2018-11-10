
function [h,p] = ttest_roi(X,ppm,y,ppmLeft,ppmRight,groups,correction)

%% demo
% X=rand(10,1000);
% ppm=linspace(0,10,1000);
% y=[0,0,0,0,0,1,1,1,1,1];
% ppmLeft=[2.31, 4.32];
% ppmRight=[2.50, 4.54];
% groups=[0,1];
% correction Either Bonferroni, 'bonferroni', or Sidak, 'sidak'
%%
if exist('correction')~=1
    beta=0.05/length(ppmLeft);
elseif strcmp(correction,'bonferroni')==1
    beta=0.05/length(ppmLeft);
elseif strcmp(correction,'sidak')==1
    beta=1-((1-0.05)^(1/length(ppmLeft)));
else
    Error('Correction must be either bonferroni or sidak')
end

for k=1:length(ppmLeft)
    [~,idxL]=min(abs(ppmLeft(k)-ppm));
    [~,idxR]=min(abs(ppmRight(k)-ppm));

    if idxL>idxR
        Xseg=X(:,idxR:idxL);
    else
        Xseg=X(:,idxL:idxR);
    end
    Xroi(:,k)=sum(Xseg,2);
end

%%
group1=min(groups);
group2=max(groups);
pos=find(y==group1);
neg=find(y==group2);
for k=1:length(ppmLeft)
    [h(k),p(k)]=ttest2(Xroi(pos,k),Xroi(neg,k),beta);
end
%%
%Manhattan plot
sig=find(h==1);
figure, n=stem(ppmLeft,abs(log10(p)),'k');
hold on; plot(ppmLeft,-log10(beta)*ones(1,length(h)),'--b','LineWidth',1.5)
hold on; m=stem(ppmLeft(sig),abs(log10(p(sig))),'fill','k');
set(m,'MarkerFaceColor','red')
ylim([.1*min(log10(p)),-1.3*min(log10(p))]);
ylabel('-Log(p)')

%Box plots
if isempty(ver('stats'))==1
	disp('Boxplot require the MATLAB Statistics Toolbox')
else
	if length(sig>0) && length(unique(y)) <=10
		figure
		subplot_size=ceil(sqrt(length(sig)));
		for k=1:length(sig)
			subplot(subplot_size, subplot_size,k)
			boxplot(X(:,sig(k)),y);
			ylabel(['Feature at ' num2str(ppmLeft(sig(k)))])
			title(['P value = ', num2str(p(sig(k)))])
		end
	end
end

