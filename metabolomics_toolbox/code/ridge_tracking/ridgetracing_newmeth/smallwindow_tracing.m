function [resstr]=smallwindow_tracing(mat,time,ppm,region,ridrefinetab,lentrain,thredseg,vislen)
% This is program to do ridge tracing for small window
% For the small window, ridge tracing are more strigent,
% ridges are assumed with no turn and can overlap
% the bourndary point will be connnected
% start training set for line direction will be lentrain
% argument
%% mat: the data matrix
%% time: time vector
%% ppm: ppm vector
%% ridrefinetab: the ridge table
%% lentrain: the length for training line
%  thredseg: the maximum distance to connect a point default 2 can be changed for ridge with different curvy
%%% vislen: the minimum length for interactive ridge picking
% return: resstr
%%% composed of parameter list (the removed ridges and the region for refine) and
%%% a refined ridrefinetab
%%
% yue wu 11/28/2018

%%%%%test
% newRidges=[-1,-6,-17,-16,-5];
% thredseg=2
% lentrain=4
%%%%%%%%%

%%default parameters
if ~exist('lentrain', 'var')
  lentrain=5;
end
if ~exist('thredseg', 'var')
  thredseg=2;
end
if ~exist('vislen', 'var')
  vislen=5;
end
%% derived parameters
indrag=matchPPMs(region,ppm);
inds=indrag(1):indrag(2);
% mat=mat(:,inds);
sizemat=size(mat);
para=struct();
% remove bad ridges
disp('start solving specific reigon');
disp('1: select bad ridges to delete');
[strres]=interridpickinnernew(ridrefinetab,mat,time,ppm,'select bad ridges to delete',vislen,[],'sw');
if strcmp(strres.clusterreturn,'C')
  refinereturndata='C';
elseif strcmp(strres.clusterreturn,'S')
  refinereturndata='S';
else
  newRidges=strres.clusterreturn;
  para.ridge_remove=ridrefinetab(ismember(ridrefinetab(:,4),-newRidges),:);
  ridrefinetab=ridrefinetab(~ismember(ridrefinetab(:,4),-newRidges),:);
  %%select region with ridge from input matrix plot
  clustsrid=ridrefinetab(:,4);
  groups=unique(clustsrid,'stable');
  % groupsele=groups(find(histc(clustsrid,groups)>vislen));
  groupsele=groups(cell2mat(arrayfun(@(x)length(find(clustsrid==x)),groups,'UniformOutput',false))>vislen);
  ind=find(ismember(clustsrid,groupsele));
  mattabhere=ridrefinetab(ind,:);
  clustshere=clustsrid(ind);
  cindallhere=ridrefinetab(ind,1);
  rindallhere=ridrefinetab(ind,2);
  ridvalallhere=ridrefinetab(ind,3);
  clusters=unique(clustshere,'stable')';
  % nclusters=length(clusters);
  disp('2: select the reigon to refine(make ridge straight)');
  fig=plotRidgesherenew(mat,ppm,time,clustshere,cindallhere,rindallhere,ridvalallhere,clusters,'select the reigon to refine(make ridge straight)');
  %%convert the index
  refinereturndata=ridrefinetab;
  para.refined_region=[];
  while 1
    fig=gcf;
    rect=getrect(fig);%xmin ymin width height; x->ppm y->time
    rectrag=[max(ppm(1),rect(1)) max(time(1),rect(2)); min(ppm(end),rect(1)+rect(3)) min(time(end),rect(2)+rect(4))];%xmin ymin xmax ymax
    para.refined_region=[para.refined_region rectrag];
    rectragind=[matchPPMs(rectrag(:,1)',ppm);matchPPMs(rectrag(:,2)',time)];
    matselec=mat(rectragind(2):rectragind(4),rectragind(1):rectragind(3));
    %%find out the ridges that cross boundary
    clust_ref=[];
    subindex=find((refinereturndata(:,1)>=rectragind(1,1)&refinereturndata(:,1)<=rectragind(1,2))&(refinereturndata(:,2)>=rectragind(2,1)&refinereturndata(:,2)<=rectragind(2,2)));
    subtab=refinereturndata(subindex,:);
    clustersnew=unique(intersect(subtab(:,4),clustshere),'stable');
    if length(clustersnew)==0
      break;
    end
    % nclustersnew=length(clustersnew);
    %% pick out bourndaries
    for clustnewsinglei=1:length(clustersnew)
      clustnewsingle=clustersnew(clustnewsinglei);
      temptab=subtab(find(subtab(:,4)==clustnewsingle),:);
      temptabfull=refinereturndata(find(refinereturndata(:,4)==clustnewsingle),:);
      temptimefull=temptabfull(:,2);
      temptime=temptab(:,2);
      boundtime=[];
      direction=[];
      if (max(temptimefull)>max(temptime) && min(temptimefull)<min(temptime))||(max(temptimefull)==max(temptime) && min(temptimefull)==min(temptime))
        warning('the box should be at the end of the ridge Not middle or cover the whole ridge');
        continue;
      end
      if max(temptimefull)>max(temptime)
        boundtime=max(temptime);
        direction= -1; % the tracing part decrease time
      else
        boundtime=min(temptime);
        direction=1;% the tracing part increase time
      end
      tempres=[temptab(find(temptab(:,2)==boundtime),:) direction];   %cind(ppm) rind(time) intensity ridge_group direction
      clust_ref=[clust_ref; tempres];
    end
    %%re-trace the ridges one by one
    for clusti=1:size(clust_ref,1)
      startpoint=clust_ref(clusti,:);
      temptab=refinereturndata(find(refinereturndata(:,4)==startpoint(1,4)),:);
      if startpoint(1,5)==1
        seqtimeind=(startpoint(:,2)+1):sizemat(1);
        startx=(startpoint(:,2)-lentrain+1):startpoint(:,2);
      else
        seqtimeind=(startpoint(:,2)-1):-1:1;
        startx=(startpoint(:,2)+lentrain-1):-1:startpoint(:,2);
      end
      starty=[];
      for startxelei=1:length(startx)
        startxele=startx(startxelei);
        startyeleind=find(temptab(:,2)==startxele);
        if length(startyeleind)>1
          disp('error in ind matching!');
        end
        starty=[starty temptab(startyeleind,1)];
        % starty
      end
      if length(startx)~=length(starty)
        warning('Some ridges cannot be extended as the prior ridges is too short');
        continue;
      end
      prex=startx;
      prey=starty;
      temptabnew=[];
      for timeindelei=1:length(seqtimeind)
        %%linear extrapolation
        timeindele=seqtimeind(timeindelei);
        midppmind=round(interp1(prex,prey,timeindele,'linear','extrap'));
        ppmrag=[midppmind-thredseg midppmind+thredseg];
        ppmind=ppmrag(1):ppmrag(2);
        indrange=[1,size(mat,2)];
        if any(ppmind<indrange(1)|ppmind>indrange(2))
          warning('the extension run out of the current window');
          break;
        end
        intenvec=mat(timeindele,ppmind);
        peakind=islocalmax(intenvec);
        % [tempmax, peakind]=max(intenvec);
        peakppmind=ppmind(peakind);
        % find(peakind)

        %%% a next try for local maximum in a small reigon
        % if length(peakppmind)==0
        %   thredseg=thredseg*2
        %   ppmrag=[midppmind-thredseg midppmind+thredseg];
        %   ppmind=ppmrag(1):ppmrag(2);
        %   intenvec=mat(timeindele,ppmind);
        %   peakind=islocalmax(intenvec);
        %   peakppmind=ppmind(peakind);
        % end
        if length(peakppmind)==0
          % newppm=midppmind;%alternative method1
          % newppm=prey(end);%alternative method2
          [tempmax, peakind]=max(intenvec);
          newppm=ppmind(peakind);
        else
          [temp,newppmind]=min(abs(peakppmind-midppmind));
          newppm=peakppmind(newppmind);
        end
        temptabnew=[temptabnew; newppm timeindele mat(timeindele,newppm) startpoint(1,4)];
        prex=[prex(2:lentrain) timeindele];
        prey=[prey(2:lentrain) newppm];
      end
      refinereturndata(find(refinereturndata(:,4)==startpoint(1,4)),:)=[];
      addontab=temptab(~(ismember(temptab(:,2),seqtimeind)),:);
      addontab=[addontab; temptabnew];
      addontab=sortrows(addontab,2);
      refinereturndata=[refinereturndata; addontab];
    end
    %% replot the ridges
    close all;
    clustsrid=refinereturndata(:,4);
    groups=unique(clustsrid,'stable');
    % groupsele=groups(find(histc(clustsrid,groups)>vislen));
    groupsele=groups(cell2mat(arrayfun(@(x)length(find(clustsrid==x)),groups,'UniformOutput',false))>vislen);
    ind=find(ismember(clustsrid,groupsele));
    mattabhere=refinereturndata(ind,:);
    clustshere=clustsrid(ind);
    cindallhere=refinereturndata(ind,1);
    rindallhere=refinereturndata(ind,2);
    ridvalallhere=refinereturndata(ind,3);
    clusters=unique(clustshere,'stable')';
    disp('plot refined ridges');
    fig=plotRidgesherenew(mat,ppm,time,clustshere,cindallhere,rindallhere,ridvalallhere,clusters,'refined ridges');
  end
end
resstr.refinereturndata=refinereturndata;
resstr.para=para;
resstr.choice=strres.choice;
% function [fig]=plotRidgesherenew(mat,ppm,time,clustshere,cindallhere,rindallhere,ridvalallhere,clusters)
% % this is derived as the last function
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
%   for clus=clusters
%       indsub=find(clustshere==clus);
%       rindallsub=rindallhere(indsub);
%       [rindallsub indsort]=sort(rindallsub);
%       indsub=indsub(indsort);
%       cindallsub=cindallhere(indsub);
%       rindallsub=rindallhere(indsub);
%       ridvalallsub=ridvalallhere(indsub);
%       fig(clus)=plot3(ppm(cindallsub),time(rindallsub)'',ridvalallsub,'linewidth',5);
%     end
% end
