function [strres]=interridpickinnernew(input,mathere,timehere,ppmhere,titlehere,vislen,defaultinput,flagstage)
%% this function plot the ridge and give interactive picking
%% click the ridge once(odd times), you will get the ridge recorded. If you click it even times, it will not be selected
%% arguement:
%%% input: the input table of ridge information
%%% mat: the intensity matrix of spectral
%%% time: time vector
%%% ppm: ppm vector
%%% titlehere: the title of the figure shown
%%% vislen: the minimum length for interactive ridge picking
%%% defaultinput: the default input for compound name and quantificability
%%% flagstage: the flag to distinguish small window retracking step ('sw'), final tracking step ('ft'), skipped from small window retracking ('ssw'). default 'sw'.
%% return: strres
%%% contains the cluster number of the selected ridges (clusterreturn), names of spectral (namevec), and whether it is good to quantify (quantifyvec), the menu choice.
%% the ridges can be removed by delte ridges, in which condition, the return value will be negative cluster number
%% this is derived from ridgeTracing_clusterPeaks_interactive_2
%% MJ&YW 11/28/2018

%%%%test%%
% input=ridrefinetab;
% mathere=mat;
% timehere=time;
% ppmhere=ppm(inds);
%%%%%%%%%%

if ~exist('titlehere', 'var')
  titlehere='figure';
end
if ~exist('defaultinput', 'var') || length(defaultinput)==0
  defaultinput=struct();
  defaultinput.compd='unknown';
  defaultinput.quan='N';
end
if ~exist('flagstage', 'var')
  flagstage='ft';
end
%% intializaiton
clustsrid=input(:,4);
groups=unique(clustsrid,'stable');
% groupsele=groups(find(histc(clustsrid,groups)>vislen));
groupsele=groups(cell2mat(arrayfun(@(x)length(find(clustsrid==x)),groups,'UniformOutput',false))>vislen);
ind=find(ismember(clustsrid,groupsele));
% ind=1:length(clustsrid);
mattabhere=input(ind,:);
clustshere=clustsrid(ind);
cindallhere=input(ind,1);
rindallhere=input(ind,2);
ridvalallhere=input(ind,3);
%% Do interactive stuff
answer=0;
clusters=unique(clustshere,'stable')';
clusterreturn=[];
namevec={};
quantifyvec=[];
close all;
choices={'Pick Final Clusters','Delete Ridge(s)','Cancel',''};
choice='start';
while ~((strcmp(choice,'Pick Final Clusters'))||(strcmp(choice,'Delete Ridge(s)'))||(strcmp(choice,'Cancel')))%
  fig=plotRidgesherenew(mathere,ppmhere,timehere,clustshere,cindallhere,rindallhere,ridvalallhere,clusters,titlehere);
  % answer = menu('Interactive Ridge Picking Menu','Pick Clusters to Join','Pick Final Clusters','Delete Ridge(s)','Cancel','box_remove');
  if ~strcmp(flagstage,'ssw')
    answer=menu('Interactive Ridge Picking Menu','Pick Final Clusters','Delete Ridge(s)','Cancel');
  else
    answer=1;
  end
  if strcmp(flagstage,'ft') && answer==2
    warning('Delete Ridge(s) only occur once');
    close(gcf);
    continue;
  end
  % Globals MUST be cleaned up after each iteration and
  % initialized ONCE in the code.
  global clickedRidges;
  clickedRidges=[];
  global lineNumber;
  lineNumber=1;
  choice=choices{answer};
  if strcmp(choice,'Pick Final Clusters') && strcmp(flagstage,'sw')
    % Quit without plotting
    close(gcf);
    clusterreturn='S';
    break;
  elseif strcmp(choice,'Pick Clusters to Join')
    warning('sorry we cannot do join ridge now');
  elseif strcmp(choice,'Pick Final Clusters')
    % Pick the final ridges
    while 1
      title('Select the ridges by clicking on them, then hit Return')
      clickedRidges=[];
      lineNumber=1;
      selectLine(gcf);
      pause();
      % Save them
      inds=clickedRidges;
      inds=inds(inds~=1)-1; % shift (surface plot = 1)
      % Count the number of odd ridges
      oddRidges=inds(find(mod(sum(inds==inds'),2)));%only ones with odd click time will be added
      clusterreturn=[clusterreturn clusters(oddRidges)];
      prompt={'ridge/compound names','quantifiable?(Y/N)'};
      titlelocal='Input';
      dims=[1 35];
      definput={defaultinput.compd,defaultinput.quan};%default result unknow compound and not quantifiable
      answerdialog=inputdlg(prompt,titlelocal,dims,definput);
      namevec=[namevec,repmat(answerdialog(1),[1,length(oddRidges)])];
      quantifyvec=[quantifyvec,repmat(answerdialog(2),[1,length(oddRidges)])];
      title('you can end by typing space or go on by typing other');
      fighere=gcf;
      waitforbuttonpress;
      chinput=fighere.CurrentCharacter;
      if chinput==' '
        break;
      end
      close(gcf);
      fig=plotRidgesherenew(mathere,ppmhere,timehere,clustshere,cindallhere,rindallhere,ridvalallhere,clusters,titlehere);
    end
  elseif strcmp(choice,'Delete Ridge(s)')
    % Delete ridge(s):
    title('Select the ridges by clicking on them, then hit Return')
    clickedRidges=[];
    lineNumber=1;
    selectLine(gcf);
    pause();
    inds=clickedRidges;
    inds=inds(inds~=1)-1; % shift (surface plot = 1)
    % Count the number of odd ridges
    oddRidges=inds(find(mod(sum(inds==inds'),2)));
    clusterreturn=clusters(oddRidges);
    clusterreturn= -clusterreturn;
    % clusters(oddRidges)=[];
  elseif strcmp(choice,'Cancel')
    % Quit without plotting
    close(gcf);
    clusterreturn='C';
    break;
  end
  clear('clickedRidges','lineNumber');
  close(gcf);
end
if strcmp(choice,'Pick Final Clusters') && ~strcmp(flagstage,'sw')
    % Plot and save the selected ridges:
    plotRidgesherenew(mathere,ppmhere,timehere,clustshere,cindallhere,rindallhere,ridvalallhere,clusters,'selected results');
end
strres.clusterreturn=clusterreturn;
strres.namevec=namevec;
strres.quantifyvec=quantifyvec;
strres.choice=choice;
