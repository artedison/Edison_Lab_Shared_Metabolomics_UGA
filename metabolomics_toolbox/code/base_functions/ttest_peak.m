
function [h,p] = ttest_peak(X,ppm,y,ppmLeft,ppmRight,groups, varargin)

%% demo
% X=rand(10,1000);
% ppm=linspace(0,10,1000);
% y=[0,0,0,0,0,1,1,1,1,1];
% ppmLeft=2.31;
% ppmRight=2.50;
% plot(ppm,X);
% groups=[0,1];
%%

[~,idxL]=min(abs(ppmLeft-ppm));
[~,idxR]=min(abs(ppmRight-ppm));

if idxL>idxR
    X=X(:,idxR:idxL);
else
    X=X(:,idxL:idxR);
end
X=sum(X,2);

%%
group1=min(groups);
group2=max(groups);
X1=X(y==group1);
X2=X(y==group2);

[h,p]=ttest2(X1,X2,varargin{:});
