%% selects the folder the .fcs files are saved in
%# build a list of file names with absolute path
fPath = uigetdir('.', 'Select directory');
if fPath==0, error('no folder selected'), 
end
fNames = dir( fullfile(fPath,'*.fcs') );
filenames = strcat(fPath, filesep, {fNames.name});
shortfilenames = {fNames.name};
%% import data from fcs files
for z= 1:length(filenames);
    DATA{1,z} = readfcs_v2(filenames{1,z});
end
%% extract the cluster row and ext and tof
[lengthdata, width] = size(DATA{1,1});
for z= 1:length(filenames);
    clusters{1,z} = DATA{1,z}(15,:);
    ext{1,z} = DATA{1,z}(1,:);
    tof{1,z} = DATA{1,z}(2,:);
    green{1,z} = DATA{1,z}(3,:);
    red{1,z} = DATA{1,z}(4,:);
    yellow{1,z} = DATA{1,z}(5,:);
end
%%
for z= 1:length(filenames);
    which_clusters{1,z} = unique(clusters{1,z});
end
%% get number of events in each node per file
number_of_clusters = 7;
for z= 1:length(filenames);
clusters1{1,z} = find_envents_per_clusters(clusters{1,z},number_of_clusters);
end
%%
for z= 1:length(filenames);
    clusters1{1,z} = clusters1{1,z}';
end
clusterss = cell2mat(clusters1);
%%
%make matrix to find idx of each cluster
for z = 1:length(clusters);
maxLength(1,z) = max(length(clusters{1,z}));
end
maxlenght = max(maxLength);
for z = 1:length(filenames);
clusters{1,z}(length(clusters{1,z})+1:maxlenght) = 0;
ext{1,z}(length(ext{1,z})+1:maxlenght) = 0;
green{1,z}(length(green{1,z})+1:maxlenght) = 0;
red{1,z}(length(red{1,z})+1:maxlenght) = 0;
yellow{1,z}(length(yellow{1,z})+1:maxlenght) = 0;
end

clusters = clusters';
ext = ext';
green = green';
red = red';
yellow = yellow';

clusters_matrix = cell2mat(clusters);
ext_matrix = cell2mat(ext);
green_matrix = cell2mat(green);
red_matrix = cell2mat(red);
yellow_matrix = cell2mat(yellow);
%%
clusters_idx = find_idx_clusters(clusters_matrix,number_of_clusters);

%% find parameter values of each cluster
for i = 1:number_of_clusters;
ext_perbin{1,i} = ext_matrix(clusters_idx{1,i});
green_perbin{1,i} = green_matrix(clusters_idx{1,i});
red_perbin{1,i} = red_matrix(clusters_idx{1,i});
yellow_perbin{1,i} = yellow_matrix(clusters_idx{1,i});
end
%%
for i = 1:number_of_clusters;
    mean_ext_perbin(1,i) = mean(ext_perbin{1,i}); 
    mean_green_perbin(1,i) = mean(green_perbin{1,i});
    mean_red_perbin(1,i) = mean(red_perbin{1,i});
    mean_yellow_perbin(1,i) = mean(yellow_perbin{1,i});
end
%%
% num_bins = (1:number_of_clusters);
% mean_ext_perbin2 = [mean_ext_perbin;num_bins]';
% mean_ext_perbin2 = sortrows(mean_ext_perbin2);
% order_bins = mean_ext_perbin2(:,2); 
%%
load('order_bins.mat');
for z = 1:length(filenames);
    clusters2{1,z} = [order_bins clusters1{1,z}];
    clusters3{1,z} = sortrows(clusters2{1,z});
    clusters3{1,z}(:,1) = [];
end
%%  
clusters1 = clusters3;      
plot_w_09_spadebins
%%
clusters1 = clusters3;  
cluster2_matrix = cell2mat(clusters1);
numworms = cluster2_matrix';
%% normalize by mass
% change the mass normalization here
l1 = 100;% 15 l1's make the same mass as one adult...
l3 = 46;
l4 = 12;
A = 1;


numworms(7,:) = mean(numworms(1:7,:));
numworms(1:6,:) =[];

numworms_last = mean(numworms(23:28,:));
numworms = [numworms;numworms_last];


% numworms(:,4) = numworms(:,1)/l1;
% numworms(:,5) = numworms(:,2)/l3;
% numworms(:,6) = numworms(:,3)/l4;
% numworms(:,7) = numworms(:,4)/A;

load('run_order.mat')
numworms = [run_order numworms];
numworms_spade = sortrows(numworms);
numworms_spade(:,1) = [];