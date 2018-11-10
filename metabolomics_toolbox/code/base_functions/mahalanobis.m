function M=mahalanobis(groups,y,scores,alpha)

% groups       Groups to compare (a subset or all of the values in 'Y') ex: [1, 3]
% Y            Your Y vector or group labels (size: 1xN). ex: [0,0,0,1,1,1,2,2,2,3,3,3]
% scores       Scores values from your model. Each column is one PC. (size: Nx2)

if length(groups)~=2
    error('must have only 2 groups')
end
if ~exist('alpha','var')==1
    alpha=0.05;
end

group1=[scores(y==groups(1),1),scores(y==groups(1),2)];
group2=[scores(y==groups(2),1),scores(y==groups(2),2)];

n1=size(group1,1);
n2=size(group2,1);
mean1=mean(group1);
mean2=mean(group2);
d=mean2-mean1;
c=(n1*cov(group1)+n2*cov(group2))/(n1+n2);
M.mdis=sqrt(d*inv(c)*d');

T2=((n1*n2)/(n1+n2))*d*inv(c)*d';
M.F=(n1+n2-2-1)/(2*(n1+n2-2))*T2;
M.G=finv((1-alpha),1,n1+n2-1);


M.X=[mean1(1),mean2(1)];
M.Y=[mean1(2),mean2(2)];
% text(mean(X),mean(Y)+0.05,['mdis= ',num2str(mdis),' '])
% line(X,Y)

%%
if M.F<M.G
    M.H=0; %not signficant difference between groups
else
    M.H=1; %sig diff
end
    

