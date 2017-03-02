%% Binning workflow
% Francesca V. Ponce
% July 2016
% PROJECT ID HERE
%% Be sure that workspace variables from wormstorter_workflow are loaded before running this script

%%
% ext = cell array raw extinction values from data files with gaps from control
% particles
% Ext = mat variable with empty spaces erased
% nExt = normalized extinction values
%% 
%create new matrix of all events without control particles
numericparameters_nocp = numericparameters;
for l = 1:length(numericparameters_nocp);
    for k = 1:length(numericparameters_nocp{l,1});
    numericparameters_nocp{l,1}{k,1}([intersectt{l,1}],:) = []; 
    end
end
%%
% taking extinction values from  numericparameters_nocp and putting them in
% new matrix (ext)
for l = 1:length(numericparameters_nocp); 
   ext{l,1} = numericparameters_nocp{l,1}{7,1}; 
end
%%
%taking TOF values from numericparameters_nocp into new matrix (tof)
for l = 1:length(numericparameters_nocp); 
   tof{l,1} = numericparameters_nocp{l,1}{6,1}; 
end
%%
% take PHYellow values from numericparameters_nocp into new matrix (phy)

for l = 1:length(numericparameters_nocp); 
   phy{l,1} = numericparameters_nocp{l,1}{15,1}; 
end

%% manipulating value matrices for processing w prexisting code blocks
ext = ext';
for j = 1:length(ext);
Ext{1,j} = cell2mat(ext{1,j});
end
%%
tof = tof';
for j = 1:length(tof);
Tof{1,j} = cell2mat(tof{1,j});
end
%%
phy = phy';
for j = 1:length(phy);
Phy{1,j} = cell2mat(phy{1,j});
end

%% take mean CP extinction values from mean_controlparticles into new double
%with correct dimensions
mean_controlparticles_ext = mean_controlparticles(:,2);
mean_controlparticles_ext = mean_controlparticles_ext';

%%
mean_controlparticles_phy = mean_controlparticles(:,11);
mean_controlparticles_phy = mean_controlparticles_phy';  

%% Normalization of extinction values by dividing all extinctions values in each file by the mean control particle extinction from that same file
for h = 1:length(mean_controlparticles_ext);
    Ext{1,h} = num2cell(Ext{1,h});
    nExt{1,h} = cellfun(@(v) v./mean_controlparticles_ext{1,h}, Ext{1,h},'UniformOutput', false);
end
%%
for h = 1:length(mean_controlparticles_phy);
    Phy{1,h} = num2cell(Phy{1,h});
    nPhy{1,h} = cellfun(@(v) v./mean_controlparticles_phy{1,h}, Phy{1,h},'UniformOutput', false);
end

    %% take log of normalized extinction values
for l = 1:length(nExt);
   nExt{1,l} = cell2mat(nExt{1,l});
   lnExt{1,l} = log(nExt{1,l});    
end
%%
for l = 1:length(nPhy);
   nPhy{1,l} = cell2mat(nPhy{1,l});
   lnPhy{1,l} = log(nPhy{1,l});    
end
%% invoke worm binning script to generate table with number of events in each defined bin and plot histograms

% make_worm_bins

%% plotting histograms
% [rows, columns] = size(m);
% for k = 1:rows;
%     subplot(3,4,k)
%     bar(m(k,:));
% end

%% plot TOF vs ext or phy
% 
% x = Tof{1,7:8};
% y1 = nExt{1,7:8};
% y2 = nPhy{1,7:8};



