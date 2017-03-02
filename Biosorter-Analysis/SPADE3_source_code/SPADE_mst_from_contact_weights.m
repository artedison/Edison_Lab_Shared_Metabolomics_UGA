function  [mst_tree,adj2] = SPADE_mst_from_contact_weights(all_idx,all_data,all_local_density, file_ind)

% although this function here computes contact weight file by file, the way
% that this function is called from the main SPADE function says that we
% look at all the downsampled files together. This makes more sense,
% especially when we have a lot of files, and the number of cells from each
% file is small, and we don't have NN edges that connect cells from the
% same file

N = length(all_idx);
if size(all_data,2)==N  % make sure that each row of data is one point/cell
    all_data = all_data';  
end

all_cluster_contact_strength = zeros(max(all_idx),max(all_idx));

for i=1:max(file_ind)
    idx = all_idx(file_ind==i);
    data = all_data(file_ind==i,:);
    local_density = all_local_density(file_ind==i);
    
    N = length(idx);
    
    fprintf('Working on cluster similarity for the %d''th input file ... \n', i);
    % build KNN graph
    fprintf('    constructing nearest neighbor graph for %d downsampled cells ... ',N);
    tic
    KNN_parameter = 5;
    ns = createns(data,'nsmethod','kdtree','Distance','euclidean');
    [idx_tmp, dist_tmp] = knnsearch(ns,data,'k',KNN_parameter+1);
    KNN_edges = unique(sort([reshape(repmat(idx_tmp(:,1),1,size(idx_tmp,2)),prod(size(idx_tmp)),1),idx_tmp(:)],2),'rows');
    [~, ~, components_size,components] = extract_connected_component_from_edge_list(KNN_edges);
    fprintf('%8d ',length(components_size));

    while length(components_size)>500 && KNN_parameter<10 % in cases where KNN=5 results in too many components, increase k till number of components become small or k too large
        KNN_parameter = KNN_parameter + 1;
        [idx_tmp, dist_tmp] = knnsearch(ns,data,'k',KNN_parameter+1);
        KNN_edges = unique(sort([reshape(repmat(idx_tmp(:,1),1,size(idx_tmp,2)),prod(size(idx_tmp)),1),idx_tmp(:)],2),'rows');
        [~, ~, components_size,components] = extract_connected_component_from_edge_list(KNN_edges);
        fprintf('\b\b\b\b\b\b\b\b\b%8d ',length(components_size));
    end
    
    while length(components_size)~=1
        [~,ind] = min(components_size);
        indicator_inside_smallest_component = find(components(:,ind));
        indicator_outside_smallest_component = setdiff((1:N)',indicator_inside_smallest_component);
        ns = createns(data(indicator_outside_smallest_component,:),'nsmethod','kdtree','Distance','cityblock');
        [idx_tmp,dist_tmp] = knnsearch(ns,data(indicator_inside_smallest_component,:),'k',1);
        [min_dist_tmp,ind] = min(dist_tmp);
        KNN_edges = [KNN_edges; sort([indicator_inside_smallest_component(ind), indicator_outside_smallest_component(idx_tmp(ind))])];
        components_to_merge = find(sum(components(KNN_edges(end,:),:),1)==1);
        components(:,components_to_merge(1)) = components(:,components_to_merge(1)) + components(:,components_to_merge(2));
        components_size(components_to_merge(1)) = components_size(components_to_merge(1)) + components_size(components_to_merge(2));
        components(:,components_to_merge(2))=[];
        components_size(:,components_to_merge(2))=[];
        
        %[~, ~, components_size2,components2] = extract_connected_component_from_edge_list(KNN_edges);
        %[isequal(components,components2), isequal(components_size,components_size2)]
        fprintf('\b\b\b\b\b\b\b\b\b%8d ',length(components_size));
    end
    toc


    fprintf('    compute pairwise cluster contact strength ... ');
    tic
    density_weight = local_density(KNN_edges(:,1))'+local_density(KNN_edges(:,2))';
    cluster_connection = sort([idx(KNN_edges(:,1))',idx(KNN_edges(:,2))'],2); 
    [B,~,J] = unique(cluster_connection,'rows');
    cluster_contact_strength = zeros(max(idx),max(idx));
    for i=1:size(B,1)
        if B(i,1)~=B(i,2)
            cluster_contact_strength(B(i,1),B(i,2)) = sum(density_weight(J==i));
            cluster_contact_strength(B(i,2),B(i,1)) = cluster_contact_strength(B(i,1),B(i,2));
        end
    end
    toc
    
    all_cluster_contact_strength = all_cluster_contact_strength + cluster_contact_strength;
end


fprintf('Build spanning tree to connect the clusters ... \n')
tic
[adj,adj2, cost_value] = max_spanning_tree_from_weight_matrix(all_cluster_contact_strength);
mst_tree = adj;
toc





% tic
% KNN_parameter = 5;
% ns = createns(data,'nsmethod','kdtree','Distance','euclidean');
% [idx_tmp, dist_tmp] = knnsearch(ns,data,'k',KNN_parameter+1);
% toc
% tic
% A = sparse(N,N);
% for i=1:N
%     A(sort(idx_tmp(i,:)),i)=1;
% end
% A = (A+A'>=1);
% toc
% [~, ~, components_size,components] = extract_connected_component(A);
% while length(components_size)~=1
%     [~,ind] = min(components_size);
%     indicator_inside_smallest_component = find(components(:,ind));
%     indicator_outside_smallest_component = setdiff((1:N)',indicator_inside_smallest_component);
%     ns = createns(data(indicator_outside_smallest_component,:),'nsmethod','kdtree','Distance','cityblock');
%     [idx_tmp,dist_tmp] = knnsearch(ns,data(indicator_inside_smallest_component,:),'k',1);
%     [min_dist_tmp,ind] = min(dist_tmp);
%     A(indicator_inside_smallest_component(ind) + (indicator_outside_smallest_component(idx_tmp(ind))-1)*N)=1;
%     A((indicator_inside_smallest_component(ind)-1)*N + indicator_outside_smallest_component(idx_tmp(ind)))=1;
%     [~, ~, components_size,components] = extract_connected_component(A);
% end
% toc
% 
% fprintf('Compute pairwise contact strength ...');
% tic
% density_weighted_A = bsxfun(@times, A, local_density(:)) + bsxfun(@times, A, local_density(:)');
% idx_matrix = sparse(N,max(idx));
% idx_matrix(sort((idx-1)*N+(1:N)))=1;
% cluster_contact_strength = (idx_matrix')* density_weighted_A * idx_matrix;
% cluster_contact_strength = cluster_contact_strength - diag(diag(cluster_contact_strength));
% toc
% 
% 
% fprintf('Build spanning tree to connect the clusters ... \n')
% tic
% [adj,adj2, cost_value] = max_spanning_tree_from_weight_matrix(cluster_contact_strength);
% mst_tree = adj;
% toc





function [adj,adj2, cost_value] = max_spanning_tree_from_weight_matrix(weight_matrix)
% [adj,adj2, cost_value] = max_spanning_tree_from_weight_matrix(weight_matrix)
% Max Spanning Tree based on weight_matrix

cost_value = 0;
[nX] = size(weight_matrix,1);
weight_matrix = weight_matrix - diag(diag(weight_matrix)); % make sure the diag element are 0, so that we don't have self-loops

components = []; active_components =[];
adj = sparse(nX,nX); 
adj2 = sparse(nX,nX); 

fprintf('creating a total of %d edges ... %6d',nX-1,0); count = 0;
Xmst = [];
for i=1:nX
    weights = weight_matrix(i,:);
    [Wmax,Wind] = max(weights);
    if Wmax==0 
        continue; 
    end
    Xmst = [Xmst; [i Wind]];
    if adj(i,Wind)==0 && adj(Wind,i)==0
        adj(i,Wind)=1;     adj(Wind,i)=1;
        adj2(i,Wind)=Wmax; adj2(Wind,i)=Wmax;
        cost_value = cost_value + Wmax;
        count = count +1; fprintf('\b\b\b\b\b\b%6d',count);
    end
    if isempty(components)
        components = sparse(zeros(nX,1)); 
        components([i,Wind],1) = 1; 
        active_components=1;
    else
        [existing_comp1] = find(components(i,:)==1 & active_components==1);
        [existing_comp2] = find(components(Wind,:)==1 & active_components==1);
        if isempty(existing_comp1) && isempty(existing_comp2)
            components = [components,zeros(nX,1)]; components([i,Wind],end) = 1; active_components = [active_components,1];
        elseif ~isempty(existing_comp1) && isempty(existing_comp2)
            components([i,Wind],existing_comp1)=1;
        elseif isempty(existing_comp1) && ~isempty(existing_comp2)
            components([i,Wind],existing_comp2)=1;
        elseif ~isempty(existing_comp1) && ~isempty(existing_comp2) && existing_comp1~=existing_comp2
            components = [components, components(:,existing_comp1)+components(:,existing_comp2)];
            active_components = [active_components,1];
            active_components([existing_comp1, existing_comp2])=0;
        end
    end
end
 
while sum(active_components)>1
%     sum(active_components)
    components_sizes = sum(components);   components_sizes(active_components==0) = max(components_sizes+1);  % make sure that inactive components are bigger than all active components
    [dummy, existing_comp1] = min(components_sizes);                                                         % find out the smallest compoent, which has to be active because of the line above
    ind1 = find(components(:,existing_comp1)==1); ind1 = ind1(:)';                                           % elements in the smallest component
    ind2 = setdiff(1:size(components,1),ind1); ind2 = ind2(:)';                                              % elements outside the smallest component
    weights = weight_matrix(ind1,ind2);
    [Wmax,Wind] = max(reshape(weights,length(ind1)*length(ind2),1));
    j = ceil(Wind/length(ind1));
    i = Wind - (j-1)*length(ind1);
    Xmst = [Xmst; [ind1(i),ind2(j)]];
    adj(ind1(i),ind2(j))=1; adj(ind2(j),ind1(i))=1; 
    adj2(ind1(i),ind2(j))=Wmax; adj2(ind2(j),ind1(i))=Wmax; 
    cost_value = cost_value + Wmax;
    [existing_comp2] = find(components(ind2(j),:)==1 & active_components==1);
    components(:,existing_comp1) = components(:,existing_comp1) + components(:,existing_comp2);
    active_components(existing_comp2)=0;
    count = count +1; fprintf('\b\b\b\b\b\b%6d',count);
end
fprintf('\n');

return




