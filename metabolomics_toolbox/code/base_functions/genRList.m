function [SRList,yVec]=genRList(grpSizes, blankSpacing, tocsy)
%
% Arguments:
% grpSizes = [10,12,15,8]
% nsamp = number of samples in each group
%

if ~exist('tocsy','var')==1
    tocsy=0;
end

ngrp=length(grpSizes);
grps=('A':char(64+ngrp));
ns=[1:sum(grpSizes)];

SListNumbers=[];
grpsList=[];
for i=1:length(grpSizes)
    SListNumbers=[SListNumbers, [1:grpSizes(i)]];
    grpsList=[grpsList, repmat(grps(i),1,grpSizes(i))];
end

SList={};
for i = 1:sum(grpSizes)
    SList=[SList,{[grpsList(i), int2str(SListNumbers(i))]}];
end

SRList=SList(randperm(length(SList)));

rgrpsList=[];
for i = 1:length(SRList)
    rgrpsList=[rgrpsList,SRList{i}(1)];
end
yVec=rgrpsList-'A';

if tocsy
    tocsyIdx=randi(sum(grpSizes),[1,tocsy])
    disp(SRList(tocsyIdx))
end

nBlank=length(SRList)/blankSpacing;
for i = 1:blankSpacing:length(SRList)+nBlank
    SRList=[SRList(1:i-1) {'Z'} SRList(i:end)];
end


