%% Parameters setting
bsize=1; % Bin size for worm counting
UL=1600; % The upper limit of TOF and Ext value for worm counting
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
%% For network analysis
%% Retrive binned NMR features of interest  
featuresR.intensities=featuresR.intensities(run_order,:); %Reorder NMR features
featuresNew.bounds=featuresR.bounds(:,idx); % idx: index of interesting binned NMR features, manually picked
featuresNew.intensities=featuresR.intensities(:,idx);
featuresNew.regionsIndicesList=featuresR.regionsIndicesList(:,idx);
featuresNew.peaks.peakmatrix=featuresR.peaks.peakmatrix(:,idx);
featuresNew.peaks.shifts=featuresR.peaks.shifts(:,idx);
featuresNew.peaks.sameNames=featuresR.peaks.sameNames;
%% Draw Biosorter regions and calculate total worm mass in the regions
% Plot worm reads from different TP (in different colors) in one figure
figure,
hold on
for ind=29:sample
    plot(TOF3{ind},log10(Ext3{ind}),'k.','MarkerSize',3)
  set(gca,'Ylim',[0,log10(1600)],'Xlim',[0,1600],'FontSize',20)
end
for ind=22:28
    plot(TOF3{ind},log10(Ext3{ind}),'m.','MarkerSize',3)
end
for ind=15:21
    plot(TOF3{ind},log10(Ext3{ind}),'g.','MarkerSize',3)
end
for ind=8:14
    plot(TOF3{ind},log10(Ext3{ind}),'b.','MarkerSize',3)
end
for ind=1:7
    plot(TOF3{ind},log10(Ext3{ind}),'r.','MarkerSize',3)
end

% Manually pick regions
hh=impoly;

% Get position of each region (right-click on polygon in the figure and
% select 'copy position', then paste here
pos_ws1{1}=[283.500717360115 1.70796647814964;206.982305117169 1.1586887668372;89.143950263032 1.05457730521589;87.6135820181731 1.34896143807615;122.812051649928 1.56795451252098];
pos_ws1{2}=[125.872788139646 1.64693562133715;174.844571975131 1.93772970379668;254.423720707795 2.09569192142902;626.303204208512 2.2895546430687;519.177427068388 1.78335753656508;411.800120409392 1.53221136461886;244.190246839253 1.22903860188147;337.063605930177 1.77976748616434];
pos_ws1{3}=[277.379244380679 2.14236257663857;396.747967479675 2.43674670949883;578.861788617886 2.56957857432603;944.619799139168 2.60188902793264;858.919177427068 2.18903323184812;640.076518412243 1.91259935099153;556.291390728477 1.80566130748003;646.83925346177 2.31094924537568];
pos_ws1{4}=[588.079470198675 2.61784333703591;996.78334910123 2.8669572230869;1426.67928098392 3.00439936711504;1561.95063214931 2.98862953855337;1546.53822998194 2.73895785159317;1409.75316074654 2.51900898450918;1180.49367850692 2.37633944910335;872.245635159542 2.17422427394509;973.696795791487 2.65214973354293];

% Calculate total worm mass in each region
for ind=1:4
BSNode_ws1(:,ind)=STOCSY_BSnodecal(pos_ws1{ind},Jointdata,0,1,(size(MS_g,2)+1),UL);
end

% Show regions in the figure
patch('XData',pos_ws1{1}(:,1),'YData',pos_ws1{1}(:,2),'EdgeColor','w','FaceColor','none','LineWidth',3);
patch('XData',pos_ws1{2}(:,1),'YData',pos_ws1{2}(:,2),'EdgeColor','r','FaceColor','none','LineWidth',3);
patch('XData',pos_ws1{3}(:,1),'YData',pos_ws1{3}(:,2),'EdgeColor','k','FaceColor','none','LineWidth',3);
patch('XData',pos_ws1{4}(:,1),'YData',pos_ws1{4}(:,2),'EdgeColor','k','FaceColor','none','LineWidth',3);
%% Creat correlation table for Cytoscape
clear matrix labels
n2s_NMR=@(x) [num2str(x,'%2.2f')];
n2s_Glycan1=@(x) ['OG_',num2str(x,'%d')];
n2s_Glycan2=@(x) ['NG_',num2str(x,'%d')];
 names_Glycan1=arrayfun(n2s_Glycan1,[1:29],'UniformOutput',0);
 names_Glycan2=arrayfun(n2s_Glycan2,[30:78],'UniformOutput',0);
 names_Glycan=[names_Glycan1 names_Glycan2];
 names_BS=[{'WS1'} {'WS2'} {'WS3'} {'WS4'} ];
 names_NMR=arrayfun(n2s_NMR,featuresNew.peaks.shifts,'UniformOutput',0);
 
matrix(1,:)=[mean(MS_data(1:7,:)) mean(BSNode_ws1(1:7,:)) featuresNew.intensities(1,:)];
matrix(2:29,:) = [MS_data(8:end,:) BSNode_ws1(8:end,:) featuresNew.intensities(2:end,:)];
labels = [names_Glycan, names_BS, names_NMR];


            filename = 'Network_GlycansNMRBS_0.5_ws1.csv';
            Rthresh = 0.5; 
            Pthresh = 1;
            RmatIn = [];
            PmatIn = [];
            % Refer to function correlationNetwork for details about input
            % parameters
            
            [Rmat,Pmat,RmatT,PmatT] = correlationNetwork(matrix,Rthresh,Pthresh,filename,labels,RmatIn,PmatIn,'removeSelfEdges','both',[]);