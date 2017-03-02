function idx=qtclusteuclid(G, d, sample, D )
% QT clustering algorithm as described in:
%
% Heyer, L. J., Kruglyak, S., Yooseph, S. (1999). Exploring expression
% data: Identification and analysis of coexpressed genes. Genome Research
% 9, 1106–1115.
%
% http://genome.cshlp.org/content/9/11/1106.full
% http://genome.cshlp.org/content/9/11/1106/F5.large.jpg
%
% Modified by GSS, so it takes another argument, sample. Peaks will not be grouped if
% they are in the sample sample.
% 
% if two sets A{i} have same cardinality, we pick first one
% our distance measure is Euclidean distance
%
% input:
% G-nxp data to cluster
% d-diameter threshold
% D-Euclidean distance for all 0<i<j<=n
%
% output:
% idx-nx1 vector containing cluster indices
%
% Misha Koshelev
% January 20th, 2009
% Montague Laboratory

n=size(G,1);
if n<=1
    idx=ones(n,1);
    return;
end
if nargin<4
    D=Inf*ones(n,n);
    for i=1:n
        D(i,i)=0;
        for j=i+1:n
            D(i,j)=sqrt(sum((G(i,:)-G(j,:)).^2));D(j,i)=D(i,j);
        end
    end
    for i=1:size(D,1) %make distance large if the same samples
        idx=sample(i)==sample(1:size(D,2));
        D(idx,idx)=d*2;
    end
end
%length(unique(sample([A,j])))~=length(sample([A,j]))

C=[];Ccard=0;Cdiam=0;
for i=1:n
    flag=true;
    A=[i];Acard=1;Adiam=0;
    while flag&&length(A)<n
        pts=1:n;
        pts(A)=[];
        jdiam=zeros(length(pts),1);
        for pidx=1:length(pts)
            % We only need to compute maximal distance from new point j to all
            % existing points in cluster
            jdiam(pidx)=max(D(pts(pidx),A));
        end
        [minjdiam,pidx]=min(jdiam);
        j=pts(pidx);
        
        if max(Adiam,minjdiam)>d
            flag=false;
        else
            A=[A,j];
            Acard=Acard+1;
            Adiam=max(Adiam,minjdiam);
        end
    end
    
    A=sort(A);
    if Acard>Ccard
        C=A;
        Ccard=Acard;
        Cdiam=Adiam;
    end
end

idx=ones(n,1);
GmC=1:n;GmC(C)=[];
idx(GmC)=qtclusteuclid(G(GmC,:), d, sample, D(GmC,GmC))+1;

function d=diam(G,clust,D)
% http://thesaurus.maths.org/mmkb/entry.html?action=entryByConcept&id=3279
% The largest distance between any two points in a set is called the  set's diameter.
%
% input:
% G-nxp data for our cluster
% clust-vector of row indices into G
% D-Euclidean distance for all 0<i<j<=n
%

d=0;
for k=1:length(clust)-1
    d=max(d,max(D(clust(k),clust(k+1:end))));
end
