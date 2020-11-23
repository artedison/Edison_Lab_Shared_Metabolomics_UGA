%% For network analysis
load('workspace_gly_bs_binNMR.mat')
%% Retrive binned NMR features of interest  

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
BSNode_ws1(:,ind)=STOCSY_BSnodecal(pos_ws1{ind},Jointdata,1,(size(MS_g,2)+1),UL);
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
%% stack-plot of STOCSY responds to glycan
Clear
load('workspace_stackstocsy.mat')
% Set the right and left boundries of NMR spectra to be presented
rppm=4.146; 
lppm=4.25;

rb=matchPPMs(rppm,JointPpm);
lb=matchPPMs(lppm,JointPpm);
[~,~,lines_v3]=StackSTOCSYall_v3(JointPpm(size(XALNFinal,2)+20*(order_cluster1-1)+10),JointX,JointPpm,rb,lb);
lines_stack_v3=lines_v3+permute(repmat(0.4*[0:1:size(lines_v3,3)-1],100,1,length(rb:lb)),[1,3,2]); % Change 0.4 for another number if needed

figure, plot(JointPpm(rb:lb),squeeze(lines_stack_v3(1,:,:))','Color',cmap(1,:), 'LineWidth',1.1)
hold on
for k=2:size(cmap,1)
    plot(JointPpm(rb:lb),squeeze(lines_stack_v3(k,:,:))','Color',cmap(k,:), 'LineWidth',1.1);
end
set(gca,'XDir','rev')
xlabel('Chemical Shift (ppm)')
set(gca,'YTickLabel',[]);
t=colorbar('colormap',jet(100));
caxis([-1 1])