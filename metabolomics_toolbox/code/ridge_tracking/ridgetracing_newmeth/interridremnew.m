function [resstr]=interridremnew(input,mat,time,ppm,titlehere)
%% this function plot the ridge and use box to remove bad points
%% argument:
%%% input: the input table of ridge information
%%% mat: the intensity matrix of spectral
%%% time: time vector
%%% ppm: ppm vector
%%% titlehere: the title of the figure shown
%% return: resstr
%%% composed of the ridge table after removing selected region(tabtochange) and
%%% the parameters (the points removed)
%% by YW 11/28/2018

%%test
% input=ridrefinetab;

%% intializaiton
para=[];
tabtochange=input;
clustsrid=input(:,4);
groups=unique(clustsrid,'stable');
ind=find(ismember(clustsrid,groups));
% ind=1:length(clustsrid);
mattabhere=input(ind,:);
clustshere=clustsrid(ind);
cindallhere=input(ind,1);
rindallhere=input(ind,2);
ridvalallhere=input(ind,3);
sizes=size(mat);
rowconverter=sizes(1)/(time(end)-time(1));
colconverter=sizes(2)/(ppm(end)-ppm(1));
%% Do interactive stuff
answer=0;
clusters=unique(clustshere,'stable')';
while ~or(answer==1,answer==2)
  close all;
  fig=plotRidgesherenew(mat,ppm,time,clustshere,cindallhere,rindallhere,ridvalallhere,clusters,titlehere);
  answer=menu('Interactive Ridge remove','yes','no');
  switch answer
    case 1
      width=1;
      while 1
        fighere=gcf;
        %%get the line that are needed to deal with
        title('select lines you want deal with (end by ENTER)');
        global clickedRidges;
        clickedRidges=[];
        global lineNumber;
        lineNumber=1;
        selectLine(fighere);
        pause();
        inds=clickedRidges;
        inds=inds(inds~=1)-1; % shift (surface plot = 1)
        oddRidges=inds(find(mod(sum(inds==inds'),2)));%only ones with odd click time will be added
        clusterreturn=clusters(oddRidges);
        seleccluserindex=find(ismember(tabtochange(:,4),clusterreturn));
        %%get the box draw
        while 1
          title('use box to remove bad region (end by clicking)');
          vec=getrect(fighere);
          rectrag=[max(ppm(1),vec(1)) max(time(1),vec(2)); min(ppm(end),vec(1)+vec(3)) min(time(end),vec(2)+vec(4))];%xmin ymin xmax ymax
          rectragind=[matchPPMs(rectrag(:,1)',ppm);matchPPMs(rectrag(:,2)',time)];
          % xmin=(vec(1)-ppm(1))*colconverter;
          % ymin=(vec(2)-time(1))*rowconverter;
          width=vec(3)*colconverter;
          % height=vec(4)*rowconverter;
          if width==0
            break;
          end
          rowind=rectragind(2,1):rectragind(2,2);
          colind=rectragind(1,1):rectragind(1,2);
          rowtabind=find(ismember(tabtochange(:,1),colind)&ismember(tabtochange(:,2),rowind));
          indfinal=intersect(seleccluserindex,rowtabind);
          para=[para; tabtochange(indfinal,:)];
          tabtochange(indfinal,:)=[];
        end
        title('you can end by typing space or go on by typing other');
        waitforbuttonpress;
        chinput=fighere.CurrentCharacter;
        if chinput==' '
          break;
        end
        close all;
        fig=plotRidgesherenew(mat,ppm,time,clustshere,cindallhere,rindallhere,ridvalallhere,clusters,titlehere);
      end
    case 2
      ;
    end
    close(gcf);
end
close all;
if answer == 1
    % Plot and save the selected ridges:
    clustsrid2=tabtochange(:,4);
    groups2=unique(clustsrid2,'stable');
    ind2=find(ismember(clustsrid2,groups2));
    mattabhere2=tabtochange(ind2,:);
    clustshere2=clustsrid2(ind2);
    cindallhere2=tabtochange(ind2,1);
    rindallhere2=tabtochange(ind2,2);
    ridvalallhere2=tabtochange(ind2,3);
    clusters2=unique(clustshere2,'stable')';
    plotRidgesherenew(mat,ppm,time,clustshere2,cindallhere2,rindallhere2,ridvalallhere2,clusters2,'refined ridges');
end
resstr.tabtochange=tabtochange;
resstr.para=para;
end
