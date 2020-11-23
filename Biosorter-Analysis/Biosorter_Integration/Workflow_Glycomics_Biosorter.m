%% Parameters setting
bsize=1; % Bin size for worm counting
UL=1600; % The upper limit of TOF and Ext value for worm counting
load('workspace_gly_bs_begin.mat') % Load matrix to start
%% Hierarchical clustering analysis of glycomics data

% HCA on time point-grouped samples
MS=[MSOGALL_DATA;MSNGALL_DATA];
MSG(:,1)=mean(MS(:,1:7),2);
MSG(:,2)=mean(MS(:,8:14),2);
MSG(:,3)=mean(MS(:,15:21),2);
MSG(:,4)=mean(MS(:,22:28),2);
MSG(:,5)=mean(MS(:,29:end),2);
clustergram(MSG,'ColumnLabels',MSG_label,'Colormap',redbluecmap,'Standardize','Row');

% HCA on all samples
clustergram(MS,'ColumnLabels',MS_sample,'Colormap',redbluecmap,'Standardize','Row');
%% Creat pseudo-peak of MS data
% Concatenate O-glycans and N-glycans matrix
MS_data=[MSOGALL_DATA;MSNGALL_DATA]';
% Sort MS data by Biosorter data order
MS_data=[MSOrder MS_data];
MS_data=sortrows(MS_data);
MS_data(:,1)=[];

% Creat pseudo-peak of MS data, this only affect peak visualizing
MS_g=[];

for j = 1:size(MS_data,2) 
  for ii=1:size(MS_data,1) 
        Add_peak_g(ii,:,j)= create_pseudo_peak(0,2,MS_data(ii,j)); 
  end
    MS_g = [MS_g Add_peak_g(:,1:end-1,j)];
end
%% Biosorter data processing -- Bin worm reads

sample=length(Ext3); % Sample number
%% Plot worm distribution for each sample
figure,
for ind=1:7
    subplot(5,7,ind)
    plot(TOF3{ind},log10(Ext3{ind}),'r.','MarkerSize',3)
    set(gca,'Ylim',[0,log10(1600)],'Xlim',[0,1600],'Xtick',[],'Ytick',[])
end
for ind=8:14
    subplot(5,7,ind)
    plot(TOF3{ind},log10(Ext3{ind}),'b.','MarkerSize',3)
    set(gca,'Ylim',[0,log10(1600)],'Xlim',[0,1600],'Xtick',[],'Ytick',[])
end
for ind=15:21
    subplot(5,7,ind)
    plot(TOF3{ind},log10(Ext3{ind}),'g.','MarkerSize',3)
    set(gca,'Ylim',[0,log10(1600)],'Xlim',[0,1600],'Xtick',[],'Ytick',[])
end
for ind=22:28
    subplot(5,7,ind)
    plot(TOF3{ind},log10(Ext3{ind}),'m.','MarkerSize',3)
    set(gca,'Ylim',[0,log10(1600)],'Xlim',[0,1600],'Xtick',[],'Ytick',[])
end
for ind=29:sample
    subplot(5,7,ind)
    plot(TOF3{ind},log10(Ext3{ind}),'k.','MarkerSize',3)
    set(gca,'Ylim',[0,log10(1600)],'Xlim',[0,1600],'Xtick',[],'Ytick',[])
end
for ind=29
    subplot(5,7,ind)
    plot(TOF3{ind},log10(Ext3{ind}),'k.','MarkerSize',3)
    set(gca,'Ylim',[0,log10(1600)],'Xlim',[0,1600],'FontSize',20)
end
%%
for j=1:sample
    pn(j)=length(Ext3{j}); % Number of events in each sample
end
pn_max=max(pn); % Maxium number of events across all samples

% Variable prealocate
nn=zeros(sample,UL./bsize,UL./bsize); % Number of counts in each bin each sample
binx=NaN(sample,pn_max)'; % Binning index of each events, in TOF axis
biny=NaN(sample,pn_max)'; % Binning index of each events, in EXT axis

for j=1:sample
[nn(j,:,:),gsizex,gsizey,binx(1:pn(j),j),biny(1:pn(j),j)]=histcounts2(TOF3{1,j},Ext3{1,j},0:bsize:UL,0:bsize:UL);
end
%% Biosorter data processing -- Worm counts normalization
tic

% Variable prealocate
gxmax=(UL/bsize+1);gymax=(UL/bsize+1);
dafg=NaN(sample,gxmax,gymax);
dafg2=NaN(sample,gxmax,gymax);

% Calculate mean TOF and mean EXT in each bin
for j=1:sample
    for gx=1:gxmax
        for gy=1:gymax
    dafg(j,gx,gy)=mean(Ext3{1,j}(binx(1:pn(j),j)==(gx-1) & biny(1:pn(j),j)==(gy-1)));
    dafg2(j,gx,gy)=mean(TOF3{1,j}(binx(1:pn(j),j)==(gx-1) & biny(1:pn(j),j)==(gy-1)));
        end
    end
end
toc

z=dafg.*dafg2;
nn=nn.*z(:,2:(UL/bsize+1),2:(UL/bsize+1));
nn(isnan(nn))=0;
pn_matrix=repmat(sum(sum(nn,2),3),1,UL/bsize,UL/bsize); % Total worm mass for each sample
nn=nn./pn_matrix;

%% Biosorter data processing -- Matrix linearization
nn=reshape(nn,size(nn,1),size(nn,2)*size(nn,3));
AddpeakBS=nn;
AddpeakBS(35,:)=mean(nn(29:end,:)); 
%% Combine glycomics and biosorter data
scf1=200; % Vertical scaling factor for Biosorter data, only affect visualization
scf2=100000; % Horizontal scaling factor for Biosorter data, only affect visualization

% Concantenate data matrixes
Jointdata=[MS_g fliplr(AddpeakBS.*scf1)]; 

% Concantenate X coordinates of data matrixes
ppmMS_g=(1:length(MS_g))./1000;
extrappmBS=fliplr((1:length(AddpeakBS))./scf2+ceil(size(MS_g,2)/1000));
Jointppm=[ppmMS_g extrappmBS];

figure, plot(Jointppm,Jointdata, 'LineWidth',1.1)
set(gca,'xdir','rev')
title('Glycomics data + BioSorter data')
%% STOCSY on integrated dataset

ind=27; % Row number of interested glycan
[~,~]=STOCSY_BS(ppmMS_g(200*ind-99),Jointdata,Jointppm,bsize,(size(MS_g,2)+1),UL,'figdisply','2d');
title(['STOCSY, target: ',glycans{ind},' (',num2str(ind),')'])
%% Integrated STOCSY for each glycan as a target

for ind=1:length(glycans)
[~,~]=STOCSY_BS(ppmMS_g(200*ind-99),Jointdata,Jointppm,bsize,(size(MS_g,2)+1),UL,'figdisply','2d');
title(['STOCSY, target: ',glycans{ind},' (',num2str(ind),')'])
saveas(gcf,[num2str(ind),'.fig']); % Save figure in .fig file
close
end