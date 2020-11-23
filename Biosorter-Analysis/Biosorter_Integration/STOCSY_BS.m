function [corr,covar]=STOCSY_BS(dp,X,ppm,binsize,bslocation,UL,varargin)

% Author: Sicong Zhang, modified from function STOCSY.m
% Date: 10/07/2018

% STOCSY_BS(5.19,Jointdata,Jointppm,1,15601,1600,'thresh',0.9,'figdisply',...
% 'all','plotType','plot');
%
% Description:
%       Plots correlation projection of NMR spectrum/MS pseudo-peaks and/or
%       worm distribution to target driver peak or response vector
% 
% Input:
%       dp            3 options:
%                        1) Ppm of the driver peak 
%                        2) Response vector with length equal to size(X,1)
%                        3) an empty matrix ([])if you want to choose driver
%                        peak from the spectrum
%       X             Data matrix of spectra
%       ppm           Chemical shift vector corresponding to X
%       binsize       Bin size used for counting
%       bslocation    Location of the first biosorter index
%       UL            Upper limit of TOF and Ext value for counting
%       Optional:
%           thresh    Threshold of correlation to display, default is 0
%           figdisply 'all' (default) Show both 1D and 2D STOCSY respond
%                     figures
%                     '1d' Show only the 1D STOCSY respond figure
%                     '2d' Show only the 2D STOCSY respond figure 
% Output:
%       1st figure: Plots correlation/covariance projection of MS peaks/NMR 
%                   spectrum to target chemical shift or response vector
%       2nd figure: Plots correlation of worm distribution to driver peak 
%                   or response vector on 2D biosorter map
%       corr: Correlation matrix to target driver peak or response vector
%       covar: Covariance matrix to target driver peak or response vector

if size(dp,1)==1
    x=dp;
else
figure,
plot(ppm,X,'k'), hold,
set(gca,'xdir','rev')
zoom on;
pause();
zoom off;
[x,~]=ginput(1);
% It will open a figure with the ZOOM tool enabled:
% Choose the peak, click space and point out the exact peak to be the driver
end
p = inputParser;
p.CaseSensitive=1;
addRequired(p,'x');
addRequired(p,'X',@(x) validateattributes(x,{'numeric'}, {'2d'}));
addRequired(p,'ppm',@(x) validateattributes(x,{'numeric'}, {'row'}));
addRequired(p,'binsize',@(x) validateattributes(x,{'numeric'}, {'scalar'}));
addRequired(p,'bslocation',@(x) validateattributes(x,{'numeric'}, {'scalar'}));
addRequired(p,'UL',@(x) validateattributes(x,{'numeric'}, {'scalar'}));

validfigtypes = {'1d','2d','all'};
addParameter(p, 'figdisply', 'all', @(x) any(validatestring(x,validfigtypes)));
validPlotTypes = {'plot','stem'};
addParamValue(p, 'plotType', 'plot', @(x) any(validatestring(x,validPlotTypes)));
addParameter(p, 'thresh',0,@(x) (x>=0) && (x<=1) && isnumeric(x) && isscalar(x));

parse(p,x,X,ppm,binsize,bslocation,UL,varargin{:})

cellfun(@(f) evalin('caller',[f ' = p.Results.' f ';']), fieldnames(p.Results))

if length(x)==1     % When you input a driver peak
    [~,k]=min(abs(ppm-x));
    target_vect=X(:,k);
else                % When you input a response vector
    target_vect=x;
end
    
if size(target_vect,1)~=size(X,1) && size(target_vect,2)==size(X,1)
    target_vect=target_vect';
end

corr=(zscore(target_vect')*zscore(X))./(size(X,1)-1);   
covar=(target_vect-mean(target_vect))'*(X-repmat(mean(X),size(X,1),1))./(size(X,1)-1);

corr(abs(corr)<thresh) = 0;     % Thresholding
cmap=jet(100);

lines=NaN(size(cmap,1),size(corr,2));
ind=1;
for k=-1:2/size(cmap,1):.99
lines(ind,find(corr>k))=covar(find(corr>k));
ind=ind+1;
end

%% 1D plot
if any(strcmp(figdisply,{'1d','all'}))
    if strcmp(plotType,'stem')
        figure, stem(ppm(1:bslocation-1),lines(1,1:bslocation-1),'Color',cmap(1,:),'marker','none')
        hold on
        for k=2:size(cmap,1)
            stem(ppm(1:bslocation-1),lines(k,1:bslocation-1),'Color',cmap(k,:),'marker','none');
        end
        
    elseif strcmp(plotType,'plot')
        figure, plot(ppm(1:bslocation-1),lines(1,1:bslocation-1),'Color',cmap(1,:),'LineWidth',1.3)
        hold on
        for k=2:size(cmap,1)
            plot(ppm(1:bslocation-1),lines(k,1:bslocation-1),'Color',cmap(k,:),'LineWidth',1.3);
        end
    end
    
    set(gca,'XDir','rev')
    xlabel('Chemical Shift (ppm)')
    if length(x)==1
        ylabel(['Covariance with target'])
    else
        ylabel('Covariance with Y vector')
    end
    
    
    colormap(jet);
    t=colorbar;
    if length(x)==1
        set(get(t,'ylabel'),'String', ['Correlation with target']);
    else
        set(get(t,'ylabel'),'String', ['Correlation with Y vector']);
    end
    caxis([-1 1])
end
%% 2D plot
if any(strcmp(figdisply,{'2d','all'}))
    lines_1=fliplr(lines(:,bslocation:end));
    lines_2=reshape(lines_1,size(cmap,1),UL/binsize,UL/binsize);
    lines_3=NaN(100,UL/binsize,UL/binsize);
    [X,Y]=meshgrid(binsize:binsize:UL,binsize:binsize:UL);
    lines_3(find(lines_2))=lines_2(find(lines_2));
    
    figure
    for k=1:size(cmap,1)
        [~,ind_findx,ind_findy]=ind2sub(size(lines_3(k,:,:)),find(~isnan(lines_3(k,:,:))));
        plot(X(1,ind_findx),log10(Y(ind_findy,1)),'.','Color',cmap(k,:),'MarkerSize',6)
        hold on
    end
    hold off
    
    xlabel('TOF');
    ylabel('log10 EXT');
    
    colormap(jet);
    t=colorbar;
    if length(x)==1
        set(get(t,'ylabel'),'String', ['Correlation with target']);
    else
        set(get(t,'ylabel'),'String', ['Correlation with Y vector']);
    end
    caxis([-1 1])
end
