function clusters_idx = find_idx_clusters(clusters_matrix,number_of_clusters);

% for z = 1:length(clusters);
% maxLength(1,z) = max(length(clusters{1,z}));
% end
% maxlenght = max(maxLength);
% for z = 1:length(filenames);
% clusters{1,z}(length(clusters{1,z})+1:maxlenght) = 0;
% ext{1,z}(length(ext{1,z})+1:maxlenght) = 0;
% end
% 
% clusters_matrix = cell2mat(clusters);
% ext_matrix = cell2mat(ext);

for i = 1:number_of_clusters;
    clusters_idx{1,i} = find(clusters_matrix==i);
end

