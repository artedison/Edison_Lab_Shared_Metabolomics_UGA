function [target_vect]=STOCSY_BSnodecal(hhh,data,gsize,bslocation,UL,varargin)

% Author: Sicong Zhang, modified from function STOCSY.m
% Date: 05/07/2019

% STOCSY_BSnodecal(hhh,Jointdata,0,1,15601,1600);
%
% Description:
%       Calculate total worm mass in region defined by input 'hhh'. Also,
%       plots correlation projection of worm distribution to this total
%       worm mass
% 
% Input:
%       hhh           Positions of vertices of the region
%       data          Data matrix of jointed dataset
%       gsize         Bin size used for counting
%       bslocation    Location of the first biosorter index
%       UL            Upper limit of TOF and Ext value for counting
%       Optional:
%           thresh    Threshold (absolute value) of correlation to display, default is 0
%               
% Output:
%       target_vect   Total worm mass for designated region

p = inputParser;
p.CaseSensitive=1;
addRequired(p,'hhh',@(x) validateattributes(x,{'numeric'}, {'2d'}));
addRequired(p,'data',@(x) validateattributes(x,{'numeric'}, {'2d'}));
addRequired(p,'gsize',@(x) validateattributes(x,{'numeric'}, {'scalar'}));
addRequired(p,'bslocation',@(x) validateattributes(x,{'numeric'}, {'scalar'}));
addRequired(p,'UL',@(x) validateattributes(x,{'numeric'}, {'scalar'}));

addParameter(p, 'thresh',0,@(x) (x>=-1) && (x<=1) && isnumeric(x) && isscalar(x));

parse(p,hhh,data,gsize,bslocation,UL,varargin{:})

cellfun(@(f) evalin('caller',[f ' = p.Results.' f ';']), fieldnames(p.Results))
%%
sample=size(data,1); 
[X,Y]=meshgrid(gsize:gsize:UL,gsize:gsize:UL);
data1=fliplr(data(:,bslocation:end)); % biosorter part
data_m1=reshape(data1,sample,UL./gsize,UL./gsize);
% { 
% Option for drawing region in figure that plot worm reads from all samples
% figure 
% for ii=1:sample
%    %plot3(X,log(Y),squeeze(data_m2(ii,:,:)),'.','Color','black');
%    [~,ind_findx2,ind_findy2]=ind2sub(size(data_m1(ii,:,:)),find(data_m1(ii,:,:)));
% plot(X(1,ind_findx2),log10(Y(ind_findy2,1)),'.','Color','black')
%  hold on
% end
% hold off 
% grid on
% hh=impoly;
% hhh=getPosition(hh);
% }
%% Total worm mass calculation

target_vect=zeros(1,sample);
for ii=1:sample
    data_m1_s=squeeze(data_m1(ii,:,:))';
    B=(inpolygon(X,log10(Y),hhh(:,1),hhh(:,2)) & logical(data_m1_s));
    target_vect(ii)=sum(sum(data_m1_s(B~=0)));
end
%% Plot correlation projection (Optional, block out if not wanted)

if size(target_vect,1)~=size(data,1) && size(target_vect,2)==size(data,1)
    target_vect=target_vect';
end

corr=(zscore(target_vect')*zscore(data))./(size(data,1)-1);   
covar=(target_vect-mean(target_vect))'*(data-repmat(mean(data),size(data,1),1))./(size(data,1)-1);

% Thresholding 
corr(abs(corr)<thresh) = 0;

cmap=jet(100);

lines=NaN(size(cmap,1),size(corr,2));
ind=1;
for k=-1:2/size(cmap,1):.99
lines(ind,find(corr>k))=covar(find(corr>k));
ind=ind+1;
end

cmap=jet(100);
lines_1=fliplr(lines(:,bslocation:end));
lines_2=reshape(lines_1,size(cmap,1),UL./gsize,UL./gsize);
lines_3=NaN(100,UL./gsize,UL./gsize);
[X,Y]=meshgrid(gsize:gsize:UL,gsize:gsize:UL);
lines_3(find(lines_2))=lines_2(find(lines_2));

% Plotting
figure 
for k=1:size(cmap,1)
   [~,ind_findx,ind_findy]=ind2sub(size(lines_3(k,:,:)),find(~isnan(lines_3(k,:,:))));
plot(X(1,ind_findx),log10(Y(ind_findy,1)),'.','Color',cmap(k,:))
 hold on
end
hold off

xlabel('TOF');
ylabel('log10 EXT');
colormap(jet);
t=colorbar;
caxis([-1 1])