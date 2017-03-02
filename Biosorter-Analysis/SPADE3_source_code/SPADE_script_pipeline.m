% default result file names, do not change
parameter_filename = 'SPADE_parameters.mat';
pooled_downsampled_filename = 'SPADE_pooled_downsampled_data.mat';
cluster_mst_upsample_filename='SPADE_cluster_mst_upsample_result.mat';

% set durrent directory
directoryname = cd;

% read all_fcs_filenames all_markers all_overlapping_markers
all_fcs_filenames = getfilenames(directoryname);
all_fcs_filenames = all_fcs_filenames(isInListEnd(all_fcs_filenames,'.fcs'));
all_markers = cell(0);
all_overlapping_markers = cell(0);
for i=1:length(all_fcs_filenames)
    filename = fullfile(directoryname,all_fcs_filenames{i});
    [fcshdr] = fca_readfcs_only_header(filename);
    display(['Get marker names from fcs file:'])
    display(['   ',all_fcs_filenames{i}]);
    marker_names = cell(0);
    for i=1:length(fcshdr.par), 
        marker_names{i,1} = fcshdr.par(i).name2;  
        if isequal(unique(fcshdr.par(i).name2),' ')
            marker_names{i,1} = fcshdr.par(i).name;  
        end
    end
    [C,I] = setdiff(marker_names, all_markers);
    if ~isempty(I)
        all_markers = [all_markers;marker_names(sort(I))];
    end
    if length(all_overlapping_markers)==0
        all_overlapping_markers = marker_names;
    else
        [C,IA,IB] = intersect(all_overlapping_markers, marker_names);
        all_overlapping_markers = all_overlapping_markers(sort(IA));
    end
end



% % setup file annotations 
% % NOTE: NEED TO EDIT
% file_annot = all_fcs_filenames; % make sure that all_fcs_filenames and file_annot are of the same dimension
% for i=1:length(all_fcs_filenames)
%     file_annot{i} = ['File',num2str(i,'%02d')]; % this line assumes that the total number of fcs files <100, which is reasonable 
% end
file_annot = {'unstim'; 'LPS'; 'TPO'; 'SurfFix'};
% algorithm parameters
% % NOTE: NEED TO EDIT
used_markers = {'Cell Length'; '191-DNA'; '115-CD45'; '139-CD45RA'; '142-CD19'; '144-CD11b'; ...
                '145-CD4'; '146-CD8'; '148-CD34'; '147-CD20'; '158-CD33'; 'CD123'; 'CD38'; 'CD90'; '110_114-CD3'};
file_used_to_build_SPADE_tree = {'unstim'; 'SurfFix'};
transformation_option = 1; % 0 means no transformation, 1 means arcsinh, 2 means arcsinh followed by 0-mean 1-var
arcsinh_cofactor = 5;
kernel_width_factor = 5;
density_estimation_optimization_factor = 1.5;
outlier_density = 1;
target_density_mode = 2; % 1 means using target density percentile, 2 means choose a TD such that a fixed number of cells survive downsampling
target_density = 3;
target_cell_number = 20000;
max_allowable_events = 50000;
number_of_desired_clusters = 200;
save(parameter_filename,'all_fcs_filenames','file_annot','all_markers','all_overlapping_markers', ...
    'used_markers', 'transformation_option', 'arcsinh_cofactor', 'kernel_width_factor', ...
    'density_estimation_optimization_factor', 'outlier_density', 'target_density_mode', 'target_density',...
    'target_cell_number', 'max_allowable_events', 'number_of_desired_clusters', 'file_used_to_build_SPADE_tree');



% button 1: Compute local density for each file, if a .mat file already exist, we will skip this fcs file, and use the existing one 
local_density_available = 0;
for i=1:length(all_fcs_filenames)
    fcs_filename = fullfile(directoryname,all_fcs_filenames{i});
    mat_filename = [fcs_filename(1:end-3),'mat'];
    if exist(mat_filename)==2 % if the .mat file that stores the downsampling info already exist
        local_density_available = local_density_available + 1;
    end
end
fprintf('There are %d fcs files in this directory, %d with local density computed already.\nStart to work on the remaining %d files ...\n', length(all_fcs_filenames), local_density_available, length(all_fcs_filenames)-local_density_available);

for i=1:length(all_fcs_filenames)
    fcs_filename = fullfile(directoryname,all_fcs_filenames{i});
    mat_filename = [fcs_filename(1:end-3),'mat'];
    if exist(mat_filename)==2 % if the .mat file that stores the downsampling info already exist
        continue;
    end
    fprintf('Read fcs file ... ');
    fprintf('%s\n',all_fcs_filenames{i});
    [data, marker_names] = readfcs(fcs_filename);
    
    fprintf('Data transformation options in SPADE parameters ... ');
    switch transformation_option
        case 0, 
            fprintf('No transofrmation performed\n');
        case 1, 
            fprintf(['arcsinh transformation with cofactor ',num2str(arcsinh_cofactor),'\n']);
            data = flow_arcsinh(data,arcsinh_cofactor); 
        case 2, display('2');
            fprintf(['arcsinh transformation with cofactor ',num2str(arcsinh_cofactor),', followed by 0-mean-1-var normalization \n']);
            data = SPADE_per_gene_normalization(flow_arcsinh(data,arcsinh_cofactor)); 
        otherwise, 1;
    end
    data = data(:,1:min(end,500000)); %%NOTE: we don't want one single file to be super super huge, a single file normally does not get this big
        
    [C,IA,IB] = intersect(marker_names, used_markers);
    
    fprintf('Compute local density for each cell in this file\n')
    new_data = data(:,1:min(size(data,2),2000));  %% NOTE: this should be 2000, need to alter back later
    fprintf('  calculate median min dist ...')
    tic; [min_dist,NN_ind] = compute_min_dist_downsample(new_data,data);toc
    median_min_dist = median(min_dist);
    kernel_width = median_min_dist*kernel_width_factor;
    optimizaiton_para = median_min_dist*density_estimation_optimization_factor;
    fprintf('  calculate local densities ...')
    tic; [local_density] = compute_local_density(data, kernel_width, optimizaiton_para); toc

    save(mat_filename, 'data', 'marker_names', 'used_markers', 'local_density', 'kernel_width');
    clear('data', 'local_density');
end
fprintf('Done computing local density.\n\n');



% button 2: Pool selected files
% gather the list of files used in  building SPADE tree
[c,used_file_ind,IB] = intersect(file_annot,file_used_to_build_SPADE_tree);
used_file_ind = sort(used_file_ind); 
% define variable for all pooled data
all_data=[];
tube_channel=[];
all_local_density=[];
for i=1:length(used_file_ind)
    % get file name
    fcs_filename = fullfile(directoryname,all_fcs_filenames{used_file_ind(i)});
    mat_filename = [fcs_filename(1:end-3),'mat'];
    % load file
    display(['downsampling and pooling fcs file: ',num2str(i),'/',num2str(length(used_file_ind))]);
    display(all_fcs_filenames{used_file_ind(i)});
    load(mat_filename);
    % used to normalize the local densities for the other files 
    if i==1
        RefDataSize = size(data,2); % used to normalize the local densities for the other files 
        all_marker_names = marker_names;
        used_marker_names = used_markers;
    end
    % remove outliers
    outlier_density_value = prctile(local_density,outlier_density);
    data(:,local_density<=outlier_density_value)=[];
    local_density(local_density<=outlier_density_value)=[];
    % compute target density
    switch target_density_mode
        case 1
            target_density_value = prctile(local_density,target_density);
        case 2
            num_desired_cells = target_cell_number;
            target_density_value = downsample_to_certain_num_cells(data, local_density, num_desired_cells);
        otherwise
            1;
    end
    % downsample
    keep_prob = min(1,(target_density_value./local_density));
    is_keep = rand(1,length(local_density))<keep_prob;  
    is_keep(find(sum(isnan(data))~=0))=0;
    display([num2str(sum(is_keep)),' cells keeped in this fcs file'])
    display(' ');
    data = data(:,is_keep);
    local_density = local_density(is_keep)/length(is_keep)*RefDataSize;
    % pool data
    if isequal(marker_names,all_marker_names)
        all_data = [all_data,data];
    else
        new_marker_names = setdiff(marker_names,all_marker_names);
        all_marker_names = [all_marker_names;new_marker_names];
        all_data = [all_data;repmat(NaN,length(new_marker_names),size(all_data,2))];
        data_tmp = zeros(size(all_data,1),size(data,2))+NaN;
        [C,IA,IB] = intersect(marker_names,all_marker_names);
        data_tmp(IB,:) = data(IA,:);
        all_data = [all_data, data_tmp];
    end
    all_local_density = [all_local_density,local_density];
    tube_channel = [tube_channel,repmat(used_file_ind(i),1,size(data,2))];
end
all_data = [all_data;tube_channel];
all_marker_names{end+1} = 'FileInd';

data = all_data; clear('all_data');
marker_names = all_marker_names;
local_density = all_local_density;
display(['PooledDownsampledData has ', num2str(size(data,2)), ' cells from ', num2str(length(used_file_ind)), ' files']);
if size(data,2)>max_allowable_events
    display(['Since the number of cells exceeds the max number of allowable events ', num2str(max_allowable_events),', uniform downsampling is performed']);
    keep_ind = sort(randsample(1:size(data,2),max_allowable_events));
    data = data(:,keep_ind);
    local_density = local_density(keep_ind);
    fprintf('Now, PooledDownsampledData has %d cells from %d files. \n\n', size(data,2),length(used_file_ind));
end
save(fullfile(directoryname, pooled_downsampled_filename), 'data', 'local_density', 'marker_names', 'used_markers');



% button 3 clustering mst upsample
load(fullfile(directoryname, pooled_downsampled_filename));
[C,IA,IB] = intersect(marker_names,used_markers); IA = sort(IA);
% clustering and mst
[mst_tree, idx] = SPADE_cluster_cells(data(IA,:), number_of_desired_clusters);
data(:,idx==0)=[];
local_density(idx==0)=[];
idx(idx==0)=[];
% layout
disp('Working on visualization of the tree structure ... ');
disp('This may take a few minutes if the number of nodes is large ');
node_positions = arch_layout(mst_tree);
% determine initial node_size
node_size = zeros(1,size(node_positions,2));
for i=1:length(node_size), node_size(i) = sum(local_density(idx==i)); end
node_size = flow_arcsinh(node_size, median(node_size)/2);
node_size = ceil(node_size/max(node_size)*10);
node_size(node_size<5)=5;
node_size = node_size * 1.2;
% initialize annotations
tree_annotations = [];
tree_bubble_contour = [];
save(fullfile(directoryname, cluster_mst_upsample_filename), 'data', 'local_density', 'marker_names', 'used_markers','idx','mst_tree','node_positions','node_size','tree_annotations','tree_bubble_contour');
% upsample
clustered_data = data;
clustered_data_idx = idx;
[C,IA1,IB] = intersect(marker_names,used_markers);
for i=1:length(all_fcs_filenames)
    % get file name
    fprintf('upsampling %d of %d files\n',i, length(all_fcs_filenames));
    fcs_filename = fullfile(directoryname,all_fcs_filenames{i});
    mat_filename = [fcs_filename(1:end-3),'mat'];
    load(mat_filename,'data','marker_names');
    [C,IA2,IB] = intersect(marker_names,used_markers);
    tic;[min_dist,NN_index] = compute_min_dist_upsample(data(IA2,:),clustered_data(IA1,:));toc
    all_assign{i} = clustered_data_idx(NN_index);
end

% build the cell arrary that stores the average protein expression per node per marker per file
fprintf('\nComputing average protein expression per cluster per marker per file for %d files ... %3d',length(all_fcs_filenames)+1,1);
marker_node_average=cell(0); counter = 1;
% get the average from the pooled data
load(fullfile(directoryname, cluster_mst_upsample_filename),'data','marker_names');
for j=1:length(marker_names)
    marker_node_average{counter,1} = 'POOLED';
    marker_node_average{counter,2} = marker_names{j};
    [group_avg, counts, group_idx_values] = SPADE_compute_one_marker_group_mean(data(j,:), idx);
    group_avg(group_idx_values==0)=[];
    counts(group_idx_values==0)=[];
    group_idx_values(group_idx_values==0)=[];
    tmp = zeros(1,max(idx))+NaN;
    tmp(group_idx_values) = group_avg;
    marker_node_average{counter,3} = tmp;
    counter = counter + 1;
end
marker_node_average{counter,1} = 'POOLED';
marker_node_average{counter,2} = 'CellFreq';
[dummy, tmp] = SPADE_compute_one_marker_group_mean(ones(1,length(idx)),idx);    
marker_node_average{counter,3} = tmp(:)';
counter = counter + 1;
% get the average from individual files
for i=1:length(all_fcs_filenames)
    fprintf('\b\b\b%3d',i+1);
    load(fullfile(directoryname, [all_fcs_filenames{i}(1:end-3),'mat']),'data','marker_names');
    for j=1:length(marker_names)
        marker_node_average{counter,1} = file_annot{i};
        marker_node_average{counter,2} = marker_names{j};
        [group_avg, counts, group_idx_values] = SPADE_compute_one_marker_group_mean(data(j,:), all_assign{i}); % the following few lines are for the purpose that: some file may not have any cell belong to one particular node, and therefore, the "group_avg" does not have information for every node
        group_avg(group_idx_values==0)=[];
        counts(group_idx_values==0)=[];
        group_idx_values(group_idx_values==0)=[];
        tmp = zeros(1,max(idx))+NaN;
        tmp(group_idx_values) = group_avg;
        marker_node_average{counter,3} = tmp;
        counter = counter + 1;
    end    
    marker_node_average{counter,1} = file_annot{i};
    marker_node_average{counter,2} = 'CellFreq';
    [group_avg, counts, group_idx_values] = SPADE_compute_one_marker_group_mean(ones(1,length(all_assign{i})),all_assign{i}); 
    marker_node_average{counter,3} = zeros(1,max(idx));
    marker_node_average{counter,3}(group_idx_values)=counts;
    counter = counter + 1;
end
save(fullfile(directoryname, cluster_mst_upsample_filename), 'all_assign','all_fcs_filenames','file_annot', 'marker_node_average', '-append');
fprintf('\b\b\bDone\n\n')



