function [returndata]=ridgetrace_power2_ext(mat,ppm,time,region,path,thredseg,maxaddon,varargin)
%% this script use orthogonal polynomial, H-K sign graph to trce peaks in time series spectral
%%  the ridge will be traced mostly on the time direction
%%  the ridge will be refined in small window and by removing end.
%%  argument:
%%%  mat: the input spectral matrix size(mat)=[time ppm]
%%%  ppm: the ppm of the spectral, expected to be a size(ppm)=[1 length(spectra)]
%%%  time: the time of each sample, expected size(time)=[length(time) 1]
%%%  region: the region to do ridge tracing, expected to be a size(time)=[sample_number 1]
%%%  thredseg: the maximum distance to connect a segment
%%%%%% default 5 can be changed for ridge with different curvy
%%%  maxaddon: the addon nth maximum as ridge points.
%%%%%% this is for the condition that many peak points are not classified as ridges default 1. can be changed when the high points are not classified as ridge
%%% the following parameters usually need no changes. They can be varied by varargin by ('name' value)
%%%  windsize: the window size for H-K sign classification
%%%%%% default 7 most time no variation is needed
%%%  threhold_K: the threhold for K in the H-K sign classification,
%%%%%% default 1 most time no variation is needed
%%%  threhold_H: the threhold for H in the H-K sign classification,
%%%%%% default 0 most time no variation is needed
%%%  lengthrid: the minimum length of ridge.  default 5. This parameter affect all manual steps.
%%%%%% Ridges smaller than the threhold will not be shown and so never be selected and returned.
%%%  peakflag:  whether use manual peak picking for  ridges. default, true
%%%  remflag: whether use manual boundary removing for ridges. default, true
%%% flagsmallwid: whether do a small window refine. true need a refine false not do a refine. default false
%%% smalwid_lentrain: the prior data length for small window refine. default 5
%%% smalwid_thredseg: the thredseg for the small window ridge tracing. default 2
%%% totalautoflag: whether the ridge tracing process is totally automatic. default false
%%% defaultinput: the default input for compound name and quantifiablity
%% RETURN: returndata, a structure contains the parameters, result(cind(ppm) rind(time) intensity ridge_group time_vec ppm_vec), and ridge ©
%%%% The parameters include both input parameters in the function call and the record for the manual procedure.
%%%% The manual procedure:
%%%%% ridge_remove: the table storing information for deleted ridges. Colume names: cind(ppm) rind(time) intensity ridge_group
%%%%% refined_region: the region that do the peak refinement. an array [xmin ymin; xmax ymax]
%%%%% ridge_picking: the table include picked ridges. Colume names: cind(ppm) rind(time) intensity ridge_group
%%%%% end_remove: the table include removed ridge part. Colume names: cind(ppm) rind(time) intensity ridge_group
%%%% yue wu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%test%%%%%%
% i=1;
% region=regionhere;%[0.9250 0.9550]%[1.3 1.5];%[2.4 2.56]%[2.5 2.9]%[2.5 2.9]%[1.46 1.49]%[0.8 0.9];%[1.6 1.7];%[2.4 2.6];%[0.925 0.955];2.7 2.8 2 2.1
% % data=sampleData(1);
% % time=times;
% % mat=data.Xcollapsed_1h1d;
% % timesCollapsed_1h1d=data.timesCollapsed_1h1d;
% time=times';%timesCollapsed_1h1d;
% % ppm=data.ppm_1h1d;
% windsize=7;
% threhold_K=1;
% threhold_H=0;
% thredseg=1;
% lengthrid=5;
% maxaddon=1;
% flagsmallwid=false;
% smalwid_lentrain=5;
% smalwid_thredseg=2;
% totalautoflag=true;
% %%%%%%%%%%%%%%%%%%%%%%%%%
if size(ppm,1)>1||size(time,2)>1
  error('wrong dimension for ppm and time vector, please refer to the documents, ppm should be a row vector and time should be a column vector');
end

disp('*****************************');
% fixed parameters
threddisprop=0.1;% the threhold in ratio to connect segments
thedtheta=0.5;% the cos theta that the limit the largest theta variation between segments

%multiple argmument seting
opt=struct('peakflag',[],'remflag',[],'lengthrid',[],'threhold_K',[],'threhold_H',[],'windsize',[],'flagsmallwid',[],'smalwid_lentrain',[],'smalwid_thredseg',[],'totalautoflag',[],'defaultinput',[]);
optNames=fieldnames(opt);
nArgs=length(varargin);
if round(nArgs/2)~=nArgs/2
  error('Input Argument-Value in pairs')
end
for pairs=reshape(varargin,2,[])
  inpName=pairs{1};
  if any(strcmp(inpName,optNames))
    opt.(inpName)=pairs{2};
  else
    error('%s not a valid parameter',inpName)
  end
end
%% default argument parameters
if length(opt.windsize)==0%~exist('windsize', 'var')
  windsize=7;
else
  windsize=opt.windsize;
end
if length(opt.threhold_K)==0%~exist('threhold_K', 'var')
  threhold_K=1;
else
  threhold_K=opt.threhold_K;
end
if length(opt.threhold_H)==0%~exist('threhold_H', 'var')
  threhold_H=0;
else
  threhold_H=opt.threhold_H;
end
if ~exist('thredseg', 'var')
  thredseg=5;%20
end
if length(opt.lengthrid)==0%~exist('lengthrid', 'var')
  lengthrid=5;%20
else
  lengthrid=opt.lengthrid;
end
if ~exist('maxaddon', 'var')
  maxaddon=1;%20
end
if length(opt.peakflag)==0%~exist('peakflag', 'var')
  peakflag=true;
else
  peakflag=opt.peakflag;
end
if length(opt.remflag)==0%~exist('remflag', 'var')
  remflag=true;
else
  remflag=opt.remflag;
end
if length(opt.flagsmallwid)==0%~exist('flagsmallwid', 'var')
  flagsmallwid=true;
else
  flagsmallwid=opt.flagsmallwid;
end
if length(opt.smalwid_lentrain)==0%~exist('smalwid_lentrain', 'var')
  smalwid_lentrain=5;
else
  smalwid_lentrain=opt.smalwid_lentrain;
end
if length(opt.smalwid_thredseg)==0%~exist('smalwid_thredseg', 'var')
  smalwid_thredseg=2;
else
  smalwid_thredseg=opt.smalwid_thredseg;
end
if length(opt.totalautoflag)==0%~exist('totalautoflag', 'var')
  totalautoflag=false;
else
  totalautoflag=opt.totalautoflag;
end
if length(opt.defaultinput)==0%~exist('windsize', 'var')
  defaultinput=struct();
  defaultinput.compd='unknown';
  defaultinput.quan='N';
else
  defaultinput=opt.defaultinput;
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
parameters.peakflag=peakflag;
parameters.remflag=remflag;
parameters.flagsmallwid=flagsmallwid;
parameters.smalwid_lentrain=smalwid_lentrain;
parameters.smalwid_thredseg=smalwid_thredseg;
parameters.totalautoflag=totalautoflag;
returndata=struct();
%%initialtion
indreg=matchPPMs(region,ppm);
ind=indreg(1):indreg(2);
wholsize=size(mat);
matraw=mat;
ppmraw=ppm;
timeraw=time;
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
miu(2)=0;
miu(3)=sum(x.^2);%miu2
phi=[];
phi(1,:)=repmat(1,1,windsize);
phi(2,:)=x;
phi(3,:)=x.^2-repmat(miu(3)/miu(1),1,windsize);
b=[];
for i=1:3
  b(i,:)=phi(i,:)./sum(phi(i,:).^2);
end
%% H-K matrix calculation
ridpoint=zeros(size(mat));
Hmat=zeros(size(mat));
Kmat=zeros(size(mat));
for i=1:rown
  for j=halfwidplus:(coln-halfwidminus)
    tempmat=matuse(i:(i+windsize-1),(j-halfwidminus):(j+halfwidminus));
    a=zeros(3,3);
    for pi=1:3
      for pj=1:3
        summat=(b(pi,:)'*b(pj,:))'.*tempmat;
        a(pi,pj)=sum(summat(:));
      end
    end
    %f fu fv fuu fvv fuv fvu
    % yfuncvec=[tempmat(halfwidplus,halfwidplus) a(2,1) a(1,2) a(3,1) a(1,3) a(2,2) a(2,2)];
    yfuncvec=[tempmat(halfwidplus,halfwidplus) a(2,1) a(1,2) 2*a(3,1) 2*a(1,3) a(2,2) a(2,2)];
    %% vec: x xu xv xuu xvv xuv xvu the first row not used.
    surfvecmat=[0 0 yfuncvec(1); 1 0 yfuncvec(2); 0 1 yfuncvec(3); 0 0 yfuncvec(4); 0 0 yfuncvec(5); 0 0 yfuncvec(6); 0 0 yfuncvec(7)]; %% the first row is never used in calculation
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
%%plot 1
% plotRidgesherenew(mat,ppm,time,repmat(1,size(cindall)),cindall,rindall,ridvalall,[1]);
% %plot 2
% fig=figure(), hold on
%     surf(mat,'FaceColor','Interp');
%     ylabel('y')
%     zlabel('z')
%     title(['example'])
%     xlabel('x')
%     scatter3(cindall,rindall,ridvalall,'r','linewidth',3);

% filter on local maximum
%% ridge & local maximum
ridmin=zeros(size(mat));
for i = 1:size(mat,1)
  tempmat=mat(i,:);
  ridmin(i,:)=islocalmax(tempmat);
end
% plot soly local maximum
% indln=find(ridmin==1);
% [rindall cindall]=ind2sub(size(ridmin),indln);
% ridvalall=mat(indln);
% %plot 1
% plotRidgesherenew(mat,ppm,time,repmat(1,size(cindall)),cindall,rindall,ridvalall,[1]);
% %plot2
% fig=figure(), hold on
%     surf(mat,'FaceColor','Interp');
%     ylabel('y')
%     zlabel('z')
%     title(['example'])
%     xlabel('x')
%     scatter3(cindall,rindall,ridvalall,'r','linewidth',3);

% add maximum as ridge for each row by maximum
for i =1:size(ridmin,1)
  indopt=find(ridmin(i,:));
  tempmat=mat(i,indopt);
  [temp indtemp]=sort(tempmat,'descend');
  if length(indtemp)<maxaddon
    warning('there is no maxaddon number of maximum for some spectra, use the real maximal number of local maximum instead');
    maxaddon=length(indtemp);
  end
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
%%%%removed ridges that are really small/short.
newridrefinetab=[];
for class=unique(ridrefinetab(:,4),'stable')'
  indclass=find(ridrefinetab(:,4)==class);
  if length(indclass)>=lengthrid
    ridrefinetabref=ridrefinetab(indclass,:);
    newridrefinetab=[newridrefinetab; ridrefinetabref];
  end
end
ridrefinetab=newridrefinetab;
ridgenamesvec={};
quantifyvec={};
%% whether not do interactive process and make the ridge tracing process totally automatic
if ~totalautoflag
  %%% refine for selected region
  if flagsmallwid
    [resstr]=smallwindow_tracing(mat,time,ppm,region,ridrefinetab,smalwid_lentrain,smalwid_thredseg,lengthrid);
    if strcmp(resstr.refinereturndata,'C')
      warning("user stop the program");
      return;
    elseif strcmp(resstr.refinereturndata,'S')
      flagstage='ssw';
    else
      if length(fieldnames(resstr.para))~=0
        ridrefinetab=resstr.refinereturndata;
        parameters.ridge_remove=resstr.para.ridge_remove;
        parameters.refined_region=resstr.para.refined_region;
      end
      flagstage='ft';
    end
  end
  %%% manual pick
  clustsrid=ridrefinetab(:,4);
  groups=unique(clustsrid,'stable');
  % newRidges=groups(find(histc(clustsrid,groups)>lengthrid));
  newRidges=groups(cell2mat(arrayfun(@(x)length(find(clustsrid==x)),groups,'UniformOutput',false))>lengthrid);
  if peakflag
    disp('3: final pick peaks');
    [strres]=interridpickinnernew(ridrefinetab,mat,time,ppm,'final pick peaks',lengthrid,defaultinput,flagstage);
    newRidges=strres.clusterreturn;
    if strcmp(newRidges,'C')
      warning("user stop the program");
      return
    end
    ridgenames=strres.namevec;
    quantifys=strres.quantifyvec;
    newRidges=unique(newRidges,'stable');
    ridrefinetab=ridrefinetab(find(ismember(ridrefinetab(:,4),newRidges)),:);
    parameters.ridge_picking=ridrefinetab;
  end
  if remflag
    disp('4: peak end refinement');
    resstr=interridremnew(ridrefinetab,mat,time,ppm,'peak end refinement');
    ridrefinetab=resstr.tabtochange;
    parameters.end_remove=resstr.para;
  end
  %%update ridge table according to new refinement
  ridrefinetab=ridrefinetab(ismember(ridrefinetab(:,4),newRidges),:);
  ridgenamesvec=cell([1 size(ridrefinetab,1)]);
  quantifyvec=ridgenamesvec;
  for newridgesi=1:length(newRidges)
    idx=find(ridrefinetab(:,4)==newRidges(newridgesi));
    ridgenamesvec(idx)={ridgenames{newridgesi}};
    quantifyvec(idx)={quantifys{newridgesi}};
  end
  % ridrefinetab=[ridrefinetab ridgenamesvec];%add name colume
else
  ridgenamesvec=cell([1 size(ridrefinetab,1)]);
  quantifyvec=ridgenamesvec;
  ridgenamesvec(:)=num2cell(ridrefinetab(:,4));
  quantifyvec(:)={'N'}; %repmat('N',[1 length(ridgenamesvec)]);
  parameters.ridge_remove=[];
  parameters.refined_region=[];
  parameters.ridge_picking=[];
  parameters.end_remove=[];
end
%%plot the surface
close all;
clustsrid2=ridrefinetab(:,4);
groups=unique(clustsrid2,'stable');
% groupsele=groups(find(histc(clustsrid,groups)>lengthrid));
% groupsele=groups(cell2mat(arrayfun(@(x)length(find(clustsrid==x)),groups,'UniformOutput',false))>lengthrid);
% ind2=find(ismember(clustsrid2,groupsele));
% mattabhere2=ridrefinetab(ind2,:);
ind2=1:size(ridrefinetab,1);
clustshere2=clustsrid2(ind2);
cindallhere2=ridrefinetab(ind2,1);
rindallhere2=ridrefinetab(ind2,2);
ridvalallhere2=ridrefinetab(ind2,3);
clusters2=unique(clustshere2,'stable')';
plotRidgesherenew(mat,ppm,time,clustshere2,cindallhere2,rindallhere2,ridvalallhere2,clusters2,'plot ridges');
%%%%%%%%%%%%%%%%%%%%
%%%output data format
ppmind=ridrefinetab(:,1)+indreg(1)-1;
timeind=ridrefinetab(:,2);
linearind=sub2ind(wholsize,timeind,ppmind);
ridrefinetab(:,1)=ppmind;
ridrefinetab=[linearind ridrefinetab];
ridrefinetab=[ridrefinetab timeraw(timeind) ppmraw(ppmind)'];
result=ridrefinetab;%linear_ind column_ind row_ind intensity group time_vec ppm_vec
returndata.para=parameters;
returndata.result=result;
returndata.names=ridgenamesvec;
returndata.quantifyvec=quantifyvec;

%plot ridges
% clustsrid=ridrefinetab(:,4);
% groups=unique(clustsrid,'stable');
% % groupsele=groups(find(histc(clustsrid,groups)>lengthrid));
% groupsele=groups(cell2mat(arrayfun(@(x)length(find(clustsrid==x)),groups,'UniformOutput',false))>lengthrid);
% % ind=find(ismember(clustsrid,groupsele));
% % ind=find(clustsrid==201)
% ind=1:length(clustsrid);
% mattabhere=ridrefinetab(ind,:);
% clustshere=clustsrid(ind);
% cindallhere=ridrefinetab(ind,2);
% rindallhere=ridrefinetab(ind,3);
% ridvalallhere=ridrefinetab(ind,4);
% fig=figure(), hold on
%     % set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
%     surf(mat,'FaceColor','Interp');
%     ylabel('y')
%     zlabel('z')
%     title(['example'])
%     xlabel('x')
%     scatter3(cindallhere,rindallhere,ridvalallhere,10,clustshere*10,'linewidth',3);
% saveas(fig,strcat(path,'newmeth',num2str(region(1)),'_',num2str(region(2)),'.connected.fig'))
% close(fig);

%plot
%% test figure
% mu=[0 0];
% % Sigma=[.25 .3; .3 1];
% % Sigma=[1 0.9; 0.9 1]
% Sigma=[1 0.9; 0.9 1]
% x1 = -3:.2:3; x2 = -3:.2:3;
% [X1,X2] = meshgrid(x1,x2);
% F = mvnpdf([X1(:) X2(:)],mu,Sigma);
% F = reshape(F,length(x2),length(x1));
% % F=sin(pi*X1)+cos(0.5*X2);
% % mat= -F;%%some examples
% sizef=size(F);
% % F=repmat(0,sizef(1),sizef(2));
% F=F+rand(length(F))/5;
% mat=F;
% %%%%%%%%%%%%%%%%

%% test stuff
% len=0;
% names=fieldnames(segments);
% vec=[];
% for namei = 1:length(names)
%   name=names{namei};
%   len=len+length(segments.(name));
%   vec=[vec segments.(name)];
% end

%plot segments
% groups=unique(clusts,'stable');
% groupsele=groups(find(histc(clusts,groups)>5));
% groupsele=groups(cell2mat(arrayfun(@(x)length(find(clusts==x)),groups,'UniformOutput',false))>5);
% ind=find(ismember(clusts,groupsele));
% % ind=1:length(clusts);
% mattabhere=mattab(ind,:);
% clustshere=clusts(ind);
% cindallhere=cindall(ind);
% rindallhere=rindall(ind);
% ridvalallhere=ridvalall(ind);
% fig=figure(), hold on
%     % set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
%     surf(mat,'FaceColor','Interp');
%     ylabel('y')
%     zlabel('z')
%     title(['example'])
%     xlabel('x')
%     scatter3(cindallhere,rindallhere,ridvalallhere,10,clustshere*11,'linewidth',3);
% saveas(fig,strcat(path,'newmeth',num2str(region(1)),'_',num2str(region(2)),'.segment.fig'))
% close(fig);

%plot connections
% groups=unique(clustsrid,'stable');
% groupsele=groups(find(histc(clustsrid,groups)>5));
% groupsele=groups(cell2mat(arrayfun(@(x)length(find(clustsrid==x)),groups,'UniformOutput',false))>5);
% ind=find(ismember(clustsrid,groupsele));
% % ind=1:length(clusts);
% mattabhere=mattab(ind,:);
% clustshere=clustsrid(ind);
% cindallhere=cindall(ind);
% rindallhere=rindall(ind);
% ridvalallhere=ridvalall(ind);
% fig=figure(), hold on
%     % set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
%     surf(mat,'FaceColor','Interp');
%     ylabel('y')
%     zlabel('z')
%     title(['example'])
%     xlabel('x')
%     scatter3(cindallhere,rindallhere,ridvalallhere,10,clustshere*10,'linewidth',3);
% saveas(fig,strcat(path,'newmeth',num2str(region(1)),'_',num2str(region(2)),'.segment.fig'))
% close(fig);
end
