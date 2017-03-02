function [groupidx,all,idx]=groupPeaks(all, range)
%% example

% all{1}=[10 21 50 88];
% all{2}=[15 16 18 91 150];
% all{3}=[19 22 24 92 151];

%%
clear idx
for i=1:length(all)
     idx{i}=repmat(i,1,length(all{i}));
end
all=cat(1,all{:});
idx=cat(2,idx{:});
%% sort, keep indices
%range=10; %10 seconds!
[all,sortidx]=sort(all);
idx=idx(sortidx);
%%
grp=1;
groupidx=zeros(size(all));
distance=[];
for i=1:length(all)
    nearest=find(all>=all(i)-range & all<=all(i)+range); %find all within range
    nearest=nearest(idx(nearest)~=idx(i)); %only keep the ones that were in a different list
    % need to check if there are others with duplicate samples, then
    % pick the closest one of a sample
    newnearest=[];
    if length(unique(idx(nearest)))~=length(nearest)
        j=1;
        while j<=length(idx(nearest))
            if sum(idx(nearest(j))==idx(nearest))>=1
                qidx=find(idx(nearest(j))==idx(nearest));
                [~,lidx]=min(abs(nearest(qidx)-i));
                newnearest=[newnearest nearest(qidx(lidx))];
            end
            j=j+1;
        end
    nearest=unique(newnearest);    
    end
    
    if nearest
        if size(nearest,1)>size(nearest,2)
            nearest=nearest';
        end
        nearest=[i nearest];
        p=nchoosek(all(nearest),2);
        distance(grp)=mean(abs((p(:,2)-p(:,1))));
        %check to see if we might be stealing a member of a group, and only
        %do it if the distance of this group is less
        remidx=[];
        for j=1:length(nearest)
            gidx=groupidx(nearest);
            if gidx(j) & distance(gidx(j))<distance(grp) %don't steal it
                remidx=[remidx j];
            end
        end
        nearest(remidx)=[];
        groupidx(nearest)=grp;
    else
        groupidx(i)=grp;
    end
    grp=grp+1;
end
groupidx(sortidx)=groupidx;
idx(sortidx)=idx;
all(sortidx)=all;
groupidx=remgaps(groupidx);
end

function d_out=remgaps(diffcomp)
%remove gaps in the numbering in diffcomp
a=1; 
b=1;
df2=zeros(size(diffcomp));
dfl=length(diffcomp);
while a<=dfl
    df2(diffcomp==diffcomp(a))=b;
    diffcomp(diffcomp==diffcomp(a))=0;
    b=b+1;
    a=find(diffcomp,1,'first');
end
d_out=df2;
end