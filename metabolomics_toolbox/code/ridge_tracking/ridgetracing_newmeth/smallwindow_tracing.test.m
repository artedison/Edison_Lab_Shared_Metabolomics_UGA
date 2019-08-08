%% the test script for the small window function
global dirp;
dirp='/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/';
load '/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/clock/In_vivo_circadian_metabolomics/meta_analysis_1/workflow_combineSamples_1_1.mat'
%% visual parameter
horzshift= -0.0025;
vertshift=1E-4;
regionsele=[1.46 1.49; 1.02 1.05; 0.97 0.99; 0.992 1.01; 2.4 2.6; 2.5 2.9; 5.46 5.5; 6.51 6.55; 8.25 8.36; 0.8 0.9; 1.6 1.7; 0.925 0.955];
regionsele=regionsele([10 11 12],:);
sample=1:6;
% preparation of data structure
Sample=[];
Sample(1).ridges(1).parameters=[];
Sample(1).ridges(1).result=[];
for i = 2:length(sample)
    Sample(i).ridges(1)=Sample(1).ridges(1);
end
path=strcat(dirp,'Bioinformatics_modeling/spectral.related/result/');
%%%%%%%funciton internal codes
thredseg=4; % the main tune parameter
samp_i=3;
i=3;
maxaddon=6;
regionhere=regionsele(i,:);
data=sampleData(samp_i);
mat=data.Xcollapsed_1h1d;
ppm=data.ppm_1h1d;
time=data.timesCollapsed_1h1d;
regionhere=regionsele(i,:);

region=[2.33 2.44];%[2.4 2.56]%[2.5 2.9]%[2.5 2.9]%[1.46 1.49]%[0.8 0.9];%[1.6 1.7];%[2.4 2.6];%[0.925 0.955];2.7 2.8 2 2.1
% data=sampleData(1);
% mat=data.Xcollapsed_1h1d;
% timesCollapsed_1h1d=data.timesCollapsed_1h1d;
% time=timesCollapsed_1h1d;
% ppm=data.ppm_1h1d;
windsize=7;
threhod_K=1;
threhold_H=0;
thredseg=5;
maxaddon=7;
threddisprop=1.0;% the threhold in ratio to connect segments
thedtheta=0.5;% the cos theta that the limit the largest theta variation between segments
%% argument parameters
if ~exist('windsize', 'var')
  windsize=7;
end
if ~exist('threhold_K', 'var')
  threhold_K=1;
end
if ~exist('threhold_H', 'var')
  threhold_H=0;
end
if ~exist('thredseg', 'var')
  thredseg=5;%20
end
if ~exist('lengthrid', 'var')
  lengthrid=5;%20
end
if ~exist('maxaddon', 'var')
  maxaddon=1;%20
end
if ~exist('peakflag', 'var')
  peakflag=true;
end
if ~exist('remflag', 'var')
  remflag=true;
end
%% storing parameters
parameters=struct();
parameters.region=region;
parameters.threddisprop=threddisprop;
parameters.thedtheta=thedtheta;
parameters.windsize=windsize;
parameters.threhold_K=threhold_K;
parameters.threhold_H=threhold_H;
parameters.thredseg=thredseg;
parameters.lengthrid=lengthrid;
parameters.maxaddon=maxaddon;
returndata=struct();
%%initialtion
indreg=matchPPMs(region,ppm);
ind=indreg(1):indreg(2);
wholsize=size(mat);
mat=mat(:,ind);
ppm=ppm(ind);
halfwidplus=(windsize+1)/2;
halfwidminus=(windsize-1)/2;
sizes=size(mat);
rown=sizes(1);
coln=sizes(2);
matuse=[repmat(mat(1,:),(windsize-1)/2,1); mat; repmat(mat(end,:),(windsize-1)/2,1)];
%% construct the othogonal basis
x=1:windsize;
x=x-halfwidminus-1;
miu=[];
miu(1)=windsize;%miu0
miu(3)=sum(((-halfwidminus):halfwidminus).^2);%miu2
phi=[];
phi(1,:)=repmat(1,1,windsize);
phi(2,:)=x;
phi(3,:)=x.^2-repmat(miu(3)/miu(1),1,windsize);
b=[];
for i=1:3
  b(i,:)=phi(i,:)./sqrt(sum(phi(i,:).^2));
end
%% H-K matrix calculation
ridpoint=zeros(size(mat));
Hmat=zeros(size(mat));
Kmat=zeros(size(mat));
for i=1:rown
  for j=halfwidplus:(coln-halfwidminus)
    tempmat=matuse(i:(i+windsize-1),(j-halfwidminus):(j+halfwidminus));
    a=[];
    for pi=1:3
      for pj=1:3
        summat=(b(pi,:)'*b(pj,:)).*tempmat;
        a(pi,pj)=sum(summat(:));
      end
    end
    %f fu fv fuu fvv fuv fvu
    yfuncvec=[tempmat(halfwidplus,halfwidplus) a(2,1) a(1,2) a(3,1) a(1,3) a(2,2) a(2,2)];
    %% vec: x xu xv xuu xvv xuv xvu
    surfvecmat=[0 0 yfuncvec(1); 1 0 yfuncvec(2); 0 1 yfuncvec(3); 0 0 yfuncvec(4); 0 0 yfuncvec(5); 0 0 yfuncvec(6); 0 0 yfuncvec(7)];
    nv=cross(surfvecmat(2,:),surfvecmat(3,:));
    nval=nv/sqrt(sum(nv.^2));
    G=[surfvecmat(2,:)*surfvecmat(2,:)'  surfvecmat(2,:)*surfvecmat(3,:)'; surfvecmat(2,:)*surfvecmat(3,:)'  surfvecmat(3,:)*surfvecmat(3,:)'];
    B=[surfvecmat(4,:)*nval' surfvecmat(6,:)*nval'; surfvecmat(6,:)*nval' surfvecmat(5,:)*nval'];
    detG=det(G);
    invG=[G(2,2) -G(1,2); -G(2,1) G(1,1)]/detG;
    H=trace(invG*B)/2;
    K=det(B)/detG;
    Ridflag=(H<threhold_H)&(abs(K)<threhold_K);
    ridpoint(i,j)=Ridflag;
    Hmat(i,j)=H;
    Kmat(i,j)=K;
  end
end
% filter on local H value adaptively
% ridpointh=zeros(size(mat));
% for i=1:rown
%   for j=1:(coln-windsize+1)
%     tempmatH=Hmat(i,j:(j+windsize-1));
%     indmin=find(tempmatH==min(tempmatH(:)));
%     [rind cind]=ind2sub(size(tempmatH),indmin);
%     ridpointh(i,j+cind-1)=1;
%   end
% end

% plot soly ridges
indln=find(ridpoint==1);
[rindall cindall]=ind2sub(size(ridpoint),indln);
ridvalall=mat(indln);
% plotRidgesherenew(mat,ppm,time,repmat(1,size(cindall)),cindall,rindall,ridvalall,[1]);
% filter on local minimum
%% ridge & local minimum
ridmin=zeros(size(mat));
for i = 1:size(mat,1)
  tempmat=mat(i,:);
  ridmin(i,:)=islocalmax(tempmat);
end
% plot soly local minimum
% indln=find(ridmin==1);
% [rindall cindall]=ind2sub(size(ridmin),indln);
% ridvalall=mat(indln);
% plotRidgesherennew(mat,ppm,time,repmat(1,size(cindall)),cindall,rindall,ridvalall,[1]);
% add maximum as ridge for each row by maximum
for i =1:size(ridmin,1)
  indopt=find(ridmin(i,:));
  tempmat=mat(i,indopt);
  [temp indtemp]=sort(tempmat,'descend');
  ridpoint(i,indopt(indtemp(1:maxaddon)))=1;
end
ridpointfinal=(ridpoint+ridmin==2);
% ridpointfinal=ridmin;
indln=find(ridpointfinal==1);
[rindall cindall]=ind2sub(size(ridpointfinal),indln);
ridvalall=mat(indln);
% fig=figure(), hold on
%     surf(mat,'FaceColor','Interp');
%     ylabel('y')
%     zlabel('z')
%     title(['example'])
%     xlabel('x')
%     scatter3(cindall,rindall,ridvalall,'r','linewidth',3);
%plot 2
% fig=figure(), hold on
%   surf(ppm,time',mat,'FaceColor','Interp');
%   shading interp;
%   set(gca,'xdir','reverse');
%   ylabel('time(h)');
%   zlabel('intensity');
%   title(['example']);
%   xlabel('ppm');
%   set(gcf,'InvertHardCopy','off');
%   set(gca,'fontsize',20);
%   set(gca,'box','off');
%   scatter3(ppm(cindall),time(rindall),ridvalall,'r','linewidth',3);
% saveas(gcf,'1.3-1.5pointrid.pdf')
mattab=[cindall rindall ridvalall];% ppm,time,intensity
%% connection the point into segments
%% choose the closest one in the window size if the nearby window has point classified as ridge
%% window size is defined by thredseg
segments=struct();
clusts=[];
connected=[];% the used ridge point
sizeinfortab=size(mattab);
sizewind=size(ridpointfinal);
i=1;
sizeinfortabseq=1:sizeinfortab(1);
%%% connect nearby ridge point into segments
while 1
  connected=[connected i];
  tab=mattab(i,:);
  locus=tab(1:2);
  groupstore=[i];
  for j=(tab(2)+1):sizewind(1)
    seleind=find(abs(mattab(:,1)-locus(1))<=thredseg&(mattab(:,2)-locus(2)==1));
    seleind=setdiff(seleind,connected);
    [temp,minind]=min(abs(mattab(seleind,1)-locus(1)));
    locind=seleind(minind);
    if(length(locind)==0)
      break;
    end
    locind=locind(1);
    locus=mattab(locind,1:2);
    % locus(2)
    connected=[connected locind];
    groupstore=[groupstore locind];
  end
  locus=tab(1:2);
  for j=(tab(2)-1):-1:1
    seleind=find(abs(mattab(:,1)-locus(1))<=thredseg&(mattab(:,2)-locus(2)== -1));
    seleind=setdiff(seleind,connected);
    [temp,minind]=min(abs(mattab(seleind,1)-locus(1)));
    locind=seleind(minind);
    if(length(locind)==0)
      break;
    end
    locind=locind(1);
    locus=mattab(locind,1:2);
    % locus(2)
    connected=[connected locind];
    groupstore=[groupstore locind];
  end
  name=['group' num2str(i)];
  segments.(name)=groupstore;
  clusts(groupstore)=i;
  temp=setdiff(sizeinfortabseq,connected);
  if length(temp)==0
    break;
  end
  i=temp(1);
end
%scale
norma=mean(mattab);
mattabscal=mattab./norma;
%% connect the segment into traces
%% calculate segment mean height and formulate segments
hei=[];
pointtab=[];% the segment information table
segmentscell=struct2cell(segments);
%%% sort the segments by heights
for i =1:length(fieldnames(segments))
  hei(i)=mean(mattab(segmentscell{i},3));
  tabseleind=mattab(segmentscell{i},2);
  tabseleunscal=mattab(segmentscell{i},:);
  tabsele=mattabscal(segmentscell{i},:);
  timeseq=tabsele(:,2);
  [minval minind]=min(timeseq);
  [maxval maxind]=max(timeseq);
  [minraw minrawind]=min(tabseleind);
  [maxraw maxrawind]=max(tabseleind);
  pointtab(i,:)=[minval maxval tabsele(minind,1) tabsele(maxind,1) maxval-minval+1/sizewind(1) minraw maxraw tabseleunscal(minrawind,1) tabseleunscal(maxrawind,1)]; %min_time max_time min_time_ppm max_time_ppm length_time min_time_unscaled max_time_unscaled min_time_ppm_unscaled max_time_ppm_unscaled
end
[temp indsort]=sort(hei,'descend');
i=indsort(1);
rids=struct();
ridrefinetab=[];% ridsrefine=[];%cind(ppm) rind(time) intensity ridge_group
clustsrid=[];
connectedseg=[]; % used segments
sizeseqseq=1:length(segmentscell);
while 1
  connectedseg=[connectedseg i];
  timerag=pointtab(i,6:7);
  locus=pointtab(i,[2 4]);
  groupstore=[segmentscell{i}];
  ridrefinetab=[ridrefinetab; [mattab(segmentscell{i},:) repmat(i,length(segmentscell{i}),1)]];
  tempi=i;
  for j=(timerag(2)+1):sizewind(1)
    dist=sqrt((locus(1)-pointtab(:,1)).^2+(locus(2)-pointtab(:,3)).^2);
    % dist=abs(locus(2)-pointtab(:,3));
    dist(find(pointtab(:,1)-locus(1)<0))=Inf;%max(dist);
    dist(tempi)=max(dist);
    vec1=pointtab(tempi,[2 4])-pointtab(tempi,[1 3]);
    vec2=pointtab(:,[2 4])-pointtab(:,[1 3]);
    % vec2(i,:)=[];
    if norm(vec1)~=0
      cal1=vec1*vec2'/norm(vec1);
      C=num2cell(vec2,2);
      cos=cal1./cellfun(@norm,C)';
      cos(find(isnan(cos)))=0;
      dist(find(abs(cos)<(1-thedtheta)))=max(dist);
    end
    [minval minind]=min(dist);
    minind=minind(find(minval<pointtab(minind,5)*threddisprop));
    minind=setdiff(minind,connectedseg);
    if length(minind)==0
      break;
    end
    minind=minind(1);
    locus=pointtab(minind,[2 4]);
    connectedseg=[connectedseg minind];
    groupstore=[groupstore segmentscell{minind}];
    ppmrag=[pointtab(tempi,9) pointtab(minind,8)];
    recregion=mat((pointtab(tempi,7)+1):(pointtab(minind,6)-1),min(ppmrag):max(ppmrag));
    [addonmax addonmaxind]=max(recregion,[],2);
    temppart=recregion(sub2ind(size(recregion),1:length(addonmaxind),addonmaxind'))';
    if size(temppart,1)==1
      temppart=temppart';
    end
    addonpart=[addonmaxind+min(ppmrag)-1 ((pointtab(tempi,7)+1):(pointtab(minind,6)-1))' temppart repmat(i,length(addonmaxind),1)];
    ridrefinetab=[ridrefinetab; addonpart; [mattab(segmentscell{minind},:) repmat(i,length(segmentscell{minind}),1)]];
    tempi=minind;
  end
  locus=pointtab(i,[1 3]);
  tempi=i;
  for j=(timerag(1)-1):-1:1
    dist=sqrt((locus(1)-pointtab(:,2)).^2+(locus(2)-pointtab(:,4)).^2);
    % dist=abs(locus(2)-pointtab(:,4));
    dist(find(pointtab(:,2)-locus(1)>0))=Inf;%max(dist);
    dist(tempi)=max(dist);
    vec1=pointtab(tempi,[2 4])-pointtab(tempi,[1 3]);
    vec2=pointtab(:,[2 4])-pointtab(:,[1 3]);
    % vec2(i,:)=[];
    if norm(vec1)~=0
      cal1=vec1*vec2'/norm(vec1);
      C=num2cell(vec2,2);
      cos=cal1./cellfun(@norm,C)';
      cos(find(isnan(cos)))=0;
      dist(find(abs(cos)<(1-thedtheta)))=max(dist);
    end
    [minval minind]=min(dist);
    minind=minind(find(minval<pointtab(minind,5)*threddisprop));
    minind=setdiff(minind,connectedseg);
    if length(minind)==0
      break;
    end
    minind=minind(1);
    locus=pointtab(minind,[1 3]);
    connectedseg=[connectedseg minind];
    groupstore=[groupstore segmentscell{minind}];
    ppmrag=[pointtab(tempi,8) pointtab(minind,9)];
    recregion=mat((pointtab(tempi,6)-1):-1:(pointtab(minind,7)+1),min(ppmrag):max(ppmrag));
    [addonmax addonmaxind]=max(recregion,[],2);
    temppart=recregion(sub2ind(size(recregion),1:length(addonmaxind),addonmaxind'))';
    if size(temppart,1)==1
      temppart=temppart';
    end
    addonpart=[addonmaxind+min(ppmrag)-1 ((pointtab(tempi,6)-1):-1:(pointtab(minind,7)+1))' temppart repmat(i,length(addonmaxind),1)];
    ridrefinetab=[ridrefinetab; addonpart; [mattab(segmentscell{minind},:) repmat(i,length(segmentscell{minind}),1)]];
    tempi=minind;
  end
  name=['group' num2str(i)];
  rids.(name)=groupstore;
  % ridsrefine.(name)=ridrefinetab;
  clustsrid(groupstore)=i;
  temp=setdiff(sizeseqseq,connectedseg);
  if length(temp)==0
    break;
  end
  i=temp(1);
end
%% refine of each ridge
newridrefinetab=[];
for class=unique(ridrefinetab(:,4),'stable')'
  indclass=find(ridrefinetab(:,4)==class);
  ridrefinetabref=ridrefinetab(indclass,:);
  timesrefine=unique(ridrefinetabref(:,2),'stable')';
  for timerefine=timesrefine
    indtime=find(ridrefinetabref(:,2)==timerefine);
    intensityrefine=ridrefinetabref(indtime,3)';
    indintensity=find(intensityrefine==max(intensityrefine));
    newridrefinetab=[newridrefinetab; ridrefinetabref(indtime(indintensity),:)];
  end
end
ridrefinetab=newridrefinetab;
[]=smallwindow_tracing(mat,time,ppm,ridrefinetab,lentrain);
