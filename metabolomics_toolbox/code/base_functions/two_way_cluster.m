function [sample_order,variable_order]=two_way_cluster(X,algorithm,cluster_metric,Y)

% [sample_order,variable_order]=two_way_cluster(X,algorithm,cluster_metric,Y)
% 
% Two-way hierarchical clustering of data matrix 'X' using clustering
% algorithm 'algorithm' and distance metric 'cluster_metric'.  This
% function outputs order vectors for observations and variables along with
% an image of the clustered matrix.  See linkage.m for more information
% about 'algorithm' and 'cluster_metric' parameters.
%
% Arguments:
% 
% X                Data matrix of spectra 
% algorithm        Either 'average', 'centroid', 'complete', 'median', 
%                  'single', 'ward', or 'weighted'.
% cluster_metric   See pdist.m for more options.  Use 'euclidean',
%                  'correlation', or 'spearman' as defaults. 
% Y                Optional: Vector of classes

if isempty(ver('stats'))==1
    error('This function requires the Matlab Statistics Toolbox')
end

m1=linkage(X,algorithm,cluster_metric);
z=figure; h=dendrogram(m1,0,'ORIENTATION','left');
sample_order=str2num(get(gca,'YTickLabel'));
close(z)

m2=linkage(X',algorithm,cluster_metric);
z=figure; k=dendrogram(m2,0);
variable_order=str2num(get(gca,'XTickLabel'));
close(z)

%figure
ax(1) = subplot(3,3,2:3);
dendrogram(m2,0);
ax(2) = subplot(3,3,[4 7]);
dendrogram(m1,0,'ORIENTATION','left');
if exist('Y')==1
    set(ax(2),'YTickLabel',Y(sample_order))
end
ax(3) = subplot(3,3,[5 6 8 9]);
imagesc(X(sample_order,variable_order))
caxis([-3 3])
set(ax(2),'YDir','reverse')
set(ax(3),'XTick',1:length(variable_order));set(ax(3),'XTickLabel',variable_order)
set(ax(3),'YTick',1:length(sample_order));set(ax(3),'YTickLabel',sample_order)

hlink1 = linkprop([ax(3) ax(1)],'xlim');
hlink2 = linkprop([ax(3) ax(2)],'ylim');

setappdata(ax(3),'graphics_linkaxes',hlink1);setappdata(ax(1),'graphics_linkaxes',hlink1);
setappdata(ax(3),'graphics_linkaxes',hlink2);setappdata(ax(2),'graphics_linkaxes',hlink2);

% linkaxes([ax(3) ax(2) ax(1)],'xy');
% set(ax(2),'XlimMode','auto');
% set(ax(1),'YlimMode','[0 1]');


% linkaxes([ax(3) ax(1)],'x');
% linkaxes([ax(3) ax(2)],'y');