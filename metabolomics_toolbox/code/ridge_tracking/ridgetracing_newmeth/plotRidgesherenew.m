function [fig]=plotRidgesherenew(mat,ppm,time,clustshere,cindallhere,rindallhere,ridvalallhere,clusters,titlehere)
%% this function plot ridges in surface plot for visualization
%% it is derived from plotRidgeshere
%% argument:
%%% mat: the nmr spectra matrix
%%% ppm: the ppm vector
%%% time: the time vector
%%% clustshere: the cluster vector for plotting. from the ridrefinetab table
%%% cindallhere: the colume index vector for plotting. from the ridrefinetab table
%%% rindallhere: the row index vector for plotting. from the ridrefinetab table
%%% ridvalallhere: the intensity vector for plotting. from the ridrefinetab table
%%% clusters: the unique(clustshere,'stable')
%%% titlehere: the title of the figure
%%return: the figure handler
%% MJ&YW 11/28/2018

%%%%%%test%%%%%
% ppm=ppmhere;
% mat=mathere;
% time=timehere;
%%%%%%%%%%%%%%%

if ~exist('titlehere', 'var')
  titlehere='example';
end
fig=figure('units','normalized','outerposition',[0 0 1 1]); hold on; %% full screen figure
  surf(ppm,time',mat,'FaceColor','Interp');
  shading interp;
  set(gca,'xdir','reverse');
  ylabel('time(h)');
  zlabel('intensity');
  title(titlehere);
  xlabel('ppm');
  set(gcf,'InvertHardCopy','off');
  set(gca,'fontsize',20);
  set(gca,'box','off');
  for clus=clusters
      indsub=find(clustshere==clus);
      rindallsub=rindallhere(indsub);
      [rindallsub indsort]=sort(rindallsub);
      indsub=indsub(indsort);
      cindallsub=cindallhere(indsub);
      rindallsub=rindallhere(indsub);
      ridvalallsub=ridvalallhere(indsub);
      fig(clus)=plot3(ppm(cindallsub),time(rindallsub)'',ridvalallsub,'linewidth',5);
  end
