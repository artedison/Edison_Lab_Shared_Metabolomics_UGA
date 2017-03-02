function [positions] = TPE_2D_contour_only(DIST, real_positions)
% function [positions] = TPE_2D(DIST)

% figure(10);


N = size(DIST,1); % total number of points to be visualized

[adj,adj2, cost_value] = mst_from_dist_matrix(DIST); % build an mst according to the dist matrix


adj = triu(adj,1); % figure out the order of the edges during the mst construction
edges = SPADE_find_matrix_big_element(adj,1);
weights = diag(DIST(edges(:,1),edges(:,2)));  % for i=1:size(edges,1),  weights(i,1) = dist(edges(i,1),edges(i,2)); end
[~,I] = sort(weights,'ascend');
edges = edges(I,:);
weights = weights(I);

% initialize the visulation 
positions = zeros(2,N);
idx = 1:N;  % cluster assignments

% grow the visualizatoin according to the mst construction process
fprintf('merging points %d times for visualization ... %5d', N, 0);
for i=1:size(edges,1)
    fprintf('\b\b\b\b\b%5d', i);
    cluster1_id = idx(edges(i,1)); cluster1_members = find(idx==cluster1_id); cluster1_size = length(cluster1_members); 
    cluster2_id = idx(edges(i,2)); cluster2_members = find(idx==cluster2_id); cluster2_size = length(cluster2_members); 
    if cluster1_size<cluster2_size  % make sure that cluster 1 is bigger
        cluster1_id = idx(edges(i,2)); cluster1_members = find(idx==cluster1_id); cluster1_size = length(cluster1_members); 
        cluster2_id = idx(edges(i,1)); cluster2_members = find(idx==cluster2_id); cluster2_size = length(cluster2_members); 
    end    
%     subplot(1,1,1)
%     subplot(3,3,2); plot(real_positions(1,:),real_positions(2,:),'g.',real_positions(1,cluster1_members),real_positions(2,cluster1_members),'ob',real_positions(1,cluster2_members),real_positions(2,cluster2_members),'^r');  axis2 = axis;
%     for kk=1:size(edges,1)
%         line(real_positions(1,edges(kk,:)),real_positions(2,edges(kk,:)),'color',[0.8,0.8,0.8]);
%     end
%     for kk=[cluster1_members,cluster2_members]
%         subplot(3,3,5); text(real_positions(1,kk),real_positions(2,kk),num2str(kk)); axis(axis2);
%     end

    
    if cluster1_size==1 && cluster2_size==1
        new_cluster_id = max(idx)+1;
        positions(:,cluster2_members) = [weights(i);0]; % [DIST(cluster1_members,cluster2_members),0];
        positions(:,[cluster1_members,cluster2_members]) = positions(:,[cluster1_members,cluster2_members]) - repmat(mean(positions(:,[cluster1_members,cluster2_members]),2),1,length(cluster1_members)+length(cluster2_members));
        idx(cluster1_members)=new_cluster_id;
        idx(cluster2_members)=new_cluster_id;
    end
    if cluster1_size~=1 && cluster2_size==1
        [current_guess,bubble_contour,score] = set_plus_one_contour_search(positions(:,cluster1_members), DIST(cluster1_members,cluster2_members), weights(i));
        new_cluster_id = max(idx)+1;
        positions(:,cluster2_members) = current_guess;
        positions(:,[cluster1_members,cluster2_members]) = positions(:,[cluster1_members,cluster2_members]) - repmat(mean(positions(:,[cluster1_members,cluster2_members]),2),1,length(cluster1_members)+length(cluster2_members));
        idx(cluster1_members)=new_cluster_id;
        idx(cluster2_members)=new_cluster_id;
    end
    
    if cluster1_size~=1 && cluster2_size~=1  % clusters 1 and 2 both have multiple members
        
        merged_positions = set_plus_set_contour_based_search(positions(:,cluster1_members), positions(:,cluster2_members), DIST(cluster1_members,cluster2_members), weights(i),0);
        
        new_cluster_id = max(idx)+1;
        positions(:,[cluster1_members,cluster2_members]) = merged_positions;
        positions(:,[cluster1_members,cluster2_members]) = positions(:,[cluster1_members,cluster2_members]) - repmat(mean(positions(:,[cluster1_members,cluster2_members]),2),1,length(cluster1_members)+length(cluster2_members));
        idx(cluster1_members)=new_cluster_id;
        idx(cluster2_members)=new_cluster_id;        
    end
    
%     subplot(3,3,1); plot(positions(1,:),positions(2,:),'g.',positions(1,cluster1_members),positions(2,cluster1_members),'ob',positions(1,cluster2_members),positions(2,cluster2_members),'^r');  axis equal; axis1 = axis;
%     subplot(3,1,2);
%     for kk=[cluster1_members,cluster2_members]
%         subplot(3,3,4); text(positions(1,kk),positions(2,kk),num2str(kk));           axis(axis1);
%         subplot(3,3,5); text(real_positions(1,kk),real_positions(2,kk),num2str(kk)); axis(axis2);
%     end
%     edges_to_draw = edges(intersect(find(ismember(edges(:,1),[cluster1_members,cluster2_members])), find(ismember(edges(:,2),[cluster1_members,cluster2_members]))),:); 
%     for kk=1:size(edges_to_draw,1)
%         subplot(3,3,4); line(positions(1,edges_to_draw(kk,:)),positions(2,edges_to_draw(kk,:)),'color',[0.8,0.8,0.8]);
%         subplot(3,3,5); line(real_positions(1,edges_to_draw(kk,:)),real_positions(2,edges_to_draw(kk,:)),'color',[0.8,0.8,0.8]);
%     end
%     drawnow
    
    F = (positions(1,cluster1_members)'*ones(1,length(cluster2_members)) - ones(length(cluster1_members),1)*positions(1,cluster2_members)).^2 + (positions(2,cluster1_members)'*ones(1,length(cluster2_members)) - ones(length(cluster1_members),1)*positions(2,cluster2_members)).^2;
    % [sqrt(F),DIST(cluster1_members,cluster2_members)]
    % [min(sqrt(F(:))), weights(i), min(sqrt(F(:)))-weights(i)]
    
%     if weights(i)*0.99>min(sqrt(F(:)))
%         1
%     end
    
    if weights(i)<min(sqrt(F(:)))
        weights(i:end) = weights(i:end)/weights(i)*min(sqrt(F(:)));
    end
    
%     subplot(3,3, 3)
%     plot(pdist(positions(:,[cluster1_members,cluster2_members])'),squareform(DIST([cluster1_members,cluster2_members],[cluster1_members,cluster2_members])),'.')
    
end
fprintf('\n');





function merged_positions = set_plus_set_contour_based_search(positions_a, positions_b, Dab, merge_dist, isdisplay)

Na = size(positions_a,2);
Nb = size(positions_b,2);
% get the contour of set a, and the point in a that is closest to each
% point on the contour
[bubble_contour_a] = get_annotation_contour(positions_a, merge_dist);
dist_contour_to_a = sqrt((repmat(bubble_contour_a(1,:)',1,size(positions_a,2)) - repmat(positions_a(1,:),size(bubble_contour_a,2),1)).^2 + (repmat(bubble_contour_a(2,:)',1,size(positions_a,2)) - repmat(positions_a(2,:),size(bubble_contour_a,2),1)).^2);
[~,ind] = min(dist_contour_to_a,[],2);
points_in_a = positions_a(:,ind);
unique_point_ind = unique(ind);
good_ones_in_unique_point_ind = unique_point_ind(mean(Dab(unique_point_ind,:),2)<=min(mean(Dab(unique_point_ind,:),2)));
bubble_contour_a = bubble_contour_a(:,ismember(ind,good_ones_in_unique_point_ind));
points_in_a = points_in_a(:,ismember(ind,good_ones_in_unique_point_ind));
% do the same for set b
[bubble_contour_b] = get_annotation_contour(positions_b, merge_dist);
dist_contour_to_b = sqrt((repmat(bubble_contour_b(1,:)',1,size(positions_b,2)) - repmat(positions_b(1,:),size(bubble_contour_b,2),1)).^2 + (repmat(bubble_contour_b(2,:)',1,size(positions_b,2)) - repmat(positions_b(2,:),size(bubble_contour_b,2),1)).^2);
[~,ind] = min(dist_contour_to_b,[],2);
points_in_b = positions_b(:,ind);
unique_point_ind = unique(ind);
good_ones_in_unique_point_ind = unique_point_ind(mean(Dab(:,unique_point_ind),1)<=prctile(mean(Dab(:,unique_point_ind),1),50));
bubble_contour_b = bubble_contour_b(:,ismember(ind,good_ones_in_unique_point_ind));
points_in_b = points_in_b(:,ismember(ind,good_ones_in_unique_point_ind));
% get the mirror of b
bubble_contour_b_mirror = bubble_contour_b; bubble_contour_b_mirror(1,:) = -bubble_contour_b_mirror(1,:);
points_in_b_mirror = points_in_b;  points_in_b_mirror(1,:) = -points_in_b_mirror(1,:);
positions_b_mirror = positions_b;  positions_b_mirror(1,:) = -positions_b_mirror(1,:);
 
Y = Dab.^2;
score = zeros(size(bubble_contour_a,2), size(bubble_contour_b,2),2) + Inf;
best_score = Inf;
for i=1:size(bubble_contour_a,2)
    sin_theta_a = (points_in_a(2,i) - bubble_contour_a(2,i))/norm(points_in_a(:,i) - bubble_contour_a(:,i));
    cos_theta_a = (points_in_a(1,i) - bubble_contour_a(1,i))/norm(points_in_a(:,i) - bubble_contour_a(:,i));
    for j=1:size(bubble_contour_b,2)
        sin_theta_b = (bubble_contour_b(2,j) - points_in_b(2,j))/norm(bubble_contour_b(:,j) - points_in_b(:,j));
        cos_theta_b = (bubble_contour_b(1,j) - points_in_b(1,j))/norm(bubble_contour_b(:,j) - points_in_b(:,j));
        sin_a_minus_b = sin_theta_a*cos_theta_b - cos_theta_a*sin_theta_b;
        cos_a_minus_b = cos_theta_a*cos_theta_b + sin_theta_a*sin_theta_b;
        rotation_matrix = [cos_a_minus_b, -sin_a_minus_b; sin_a_minus_b, cos_a_minus_b];
        positions_b2      = rotation_matrix*(positions_b + repmat(- points_in_b(:,j),1,size(positions_b,2))) + repmat(bubble_contour_a(:,i),1,size(positions_b,2));
        bubble_contour_b2 = rotation_matrix*(bubble_contour_b + repmat(- points_in_b(:,j),1,size(bubble_contour_b,2))) + repmat(bubble_contour_a(:,i),1,size(bubble_contour_b,2));
        
        F = (positions_a(1,:)'*ones(1,Nb) - ones(Na,1)*positions_b2(1,:)).^2 + (positions_a(2,:)'*ones(1,Nb) - ones(Na,1)*positions_b2(2,:)).^2;
        score(i,j,1) = sum(sum(abs(F-Y)));
        if min(min(F))>=min(sum((bubble_contour_a - points_in_a).^2)) && sum(sum(abs(F-Y)))<best_score
            best_score = sum(sum(abs(F-Y)));
            best_positions_b = positions_b2;
        end
        
        sin_theta_b = (bubble_contour_b_mirror(2,j) - points_in_b_mirror(2,j))/norm(bubble_contour_b_mirror(:,j) - points_in_b_mirror(:,j));
        cos_theta_b = (bubble_contour_b_mirror(1,j) - points_in_b_mirror(1,j))/norm(bubble_contour_b_mirror(:,j) - points_in_b_mirror(:,j));
        sin_a_minus_b = sin_theta_a*cos_theta_b - cos_theta_a*sin_theta_b;
        cos_a_minus_b = cos_theta_a*cos_theta_b + sin_theta_a*sin_theta_b;
        rotation_matrix = [cos_a_minus_b, -sin_a_minus_b; sin_a_minus_b, cos_a_minus_b];
        positions_b2_mirror = rotation_matrix*(positions_b_mirror + repmat(- points_in_b_mirror(:,j),1,size(positions_b_mirror,2))) + repmat(bubble_contour_a(:,i),1,size(positions_b_mirror,2));
        bubble_contour_b2_mirror = rotation_matrix*(bubble_contour_b_mirror + repmat(- points_in_b_mirror(:,j),1,size(bubble_contour_b_mirror,2))) + repmat(bubble_contour_a(:,i),1,size(bubble_contour_b_mirror,2));
        
        F = (positions_a(1,:)'*ones(1,Nb) - ones(Na,1)*positions_b2_mirror(1,:)).^2 + (positions_a(2,:)'*ones(1,Nb) - ones(Na,1)*positions_b2_mirror(2,:)).^2;
        score(i,j,2) = sum(sum(abs(F-Y)));
        if min(min(F))>=min(sum((bubble_contour_a - points_in_a).^2)) && sum(sum(abs(F-Y)))<best_score
            best_score = sum(sum(abs(F-Y)));
            best_positions_b = positions_b2_mirror;
        end
        
%         if isdisplay==1
% 
%             subplot(3,3,9);
%             plot(points_in_b(1,j),points_in_b(2,j),'*',bubble_contour_b(1,j),bubble_contour_b(2,j),'*');
%             line(bubble_contour_b(1,:),bubble_contour_b(2,:));
%             for k=1:size(positions_b,2)
%                 text(positions_b(1,k),positions_b(2,k),num2str(k));
%             end
% 
%             subplot(3,3,7);
%             plot(positions_a(1,:),positions_a(2,:),'o',positions_b2(1,:),positions_b2(2,:),'.');
%             for k=1:size(positions_b2,2)
%                 text(positions_b2(1,k),positions_b2(2,k),num2str(k));
%             end
%             line(bubble_contour_a(1,:),bubble_contour_a(2,:));
%             line(bubble_contour_b2(1,:),bubble_contour_b2(2,:));
%             drawnow
% 
%             subplot(3,3,8);
%             plot(positions_a(1,:),positions_a(2,:),'o',positions_b2_mirror(1,:),positions_b2_mirror(2,:),'.');
%             for k=1:size(positions_b2,2)
%                 text(positions_b2_mirror(1,k),positions_b2_mirror(2,k),num2str(k));
%             end
%             line(bubble_contour_a(1,:),bubble_contour_a(2,:));
%             line(bubble_contour_b2_mirror(1,:),bubble_contour_b2_mirror(2,:));
%             drawnow
%             isdisplay=0;
%         end
        
    end
end
merged_positions = [positions_a, best_positions_b];
return







function [current_guess,bubble_contour,score]= set_plus_one_contour_search(positions_a, Dab, merge_dist)

bubble_contour = get_annotation_contour(positions_a, merge_dist);
dist_contour_to_a = sqrt((repmat(bubble_contour(1,:)',1,size(positions_a,2)) - repmat(positions_a(1,:),size(bubble_contour,2),1)).^2 + (repmat(bubble_contour(2,:)',1,size(positions_a,2)) - repmat(positions_a(2,:),size(bubble_contour,2),1)).^2);

% score = sum(abs(dist_contour_to_a - repmat(Dab(:)',size(dist_contour_to_a,1),1)),2); 
% [~,ind] = min(score);  % find the best point on contour
% initial_guess = bubble_contour(:,ind);
% [~,ind] = min(sqrt(sum((repmat(initial_guess,1,size(positions_a,2))-positions_a).^2,1)));
% initial_guess = (initial_guess-positions_a(:,ind))/norm(initial_guess-positions_a(:,ind))*merge_dist*1.001 + positions_a(:,ind);  % scale it away a bit if it is too close
% current_guess = initial_guess; 
% 
% subplot(3,3,7); plot(positions_a(1,:),positions_a(2,:),'o',bubble_contour(1,:),bubble_contour(2,:),initial_guess(1,:),initial_guess(2,:),'*', current_guess(1,:), current_guess(2,:),'^'); axis equal;
% subplot(3,3,8); plot(score);

score = sum(abs(dist_contour_to_a - repmat(Dab(:)',size(dist_contour_to_a,1),1)).*repmat(exp(-(Dab(:)/min(Dab(:))).^2)',size(dist_contour_to_a,1),1),2); 
[~,ind] = min(score);  % find the best point on contour
initial_guess = bubble_contour(:,ind);
[~,ind] = min(sqrt(sum((repmat(initial_guess,1,size(positions_a,2))-positions_a).^2,1)));
initial_guess = (initial_guess-positions_a(:,ind))/norm(initial_guess-positions_a(:,ind))*merge_dist*1.001 + positions_a(:,ind);  % scale it away a bit if it is too close
current_guess = initial_guess; 

% subplot(3,3,9); plot(positions_a(1,:),positions_a(2,:),'o',bubble_contour(1,:),bubble_contour(2,:),initial_guess(1,:),initial_guess(2,:),'*', current_guess(1,:), current_guess(2,:),'^'); axis equal;
% subplot(3,3,8); hold on; plot(score,'r'); hold off;
% 
% 1;




function bubble_contour = get_annotation_contour(node_positions, merge_dist)

canvas_axis = [min(node_positions(1,:))-merge_dist*1.05, max(node_positions(1,:))+merge_dist*1.05, min(node_positions(2,:))-merge_dist*1.05, max(node_positions(2,:))+merge_dist*1.05];
resolution=100;
delta_x = (canvas_axis(2)-canvas_axis(1))/resolution;
delta_y = (canvas_axis(4)-canvas_axis(3))/resolution;
x=canvas_axis(1):delta_x:canvas_axis(2);
y=canvas_axis(3):delta_y:canvas_axis(4);

grid_position_x = repmat(x(:),1,length(y));
grid_position_y = repmat(y,length(x),1);
grid_position = [grid_position_x(:),grid_position_y(:)]';
dist_grid_to_nodes = sqrt((repmat(grid_position(1,:)',1,size(node_positions,2)) - repmat(node_positions(1,:),size(grid_position,2),1)).^2 + (repmat(grid_position(2,:)',1,size(node_positions,2)) - repmat(node_positions(2,:),size(grid_position,2),1)).^2);

indicators = (sum(dist_grid_to_nodes<=merge_dist,2)~=0);  % grid points that are close to at least one node by merge_dist
outer_points = get_outer_points(grid_position', indicators, delta_x, delta_y);
bubble_contour = outer_points';







function outer_points = get_outer_points(grid_position, indicators, delta_x, delta_y)

patch_squareform = reshape(indicators,sqrt(length(indicators)), sqrt(length(indicators)));
contour_squareform = zeros(size(patch_squareform));

for i=1:size(patch_squareform,1)
    ind_tmp = (patch_squareform(i,:)~=[-1,patch_squareform(i,1:end-1)] | patch_squareform(i,:)~=[patch_squareform(i,2:end),-1]) & patch_squareform(i,:)==1;
    contour_squareform(i,ind_tmp)=1;
end
for j=1:size(patch_squareform,2)
    ind_tmp = (patch_squareform(:,j)~=[-1;patch_squareform(1:end-1,j)] | patch_squareform(:,j)~=[patch_squareform(2:end,j);-1]) & patch_squareform(:,j)==1;
    contour_squareform(ind_tmp,j)=1;
end

contour_squareform = contour_squareform(:);
outer_points = grid_position(contour_squareform==1,:);

% order them
ordered_points = outer_points(1,:); 
flag_outer_points = zeros(size(outer_points,1),1);flag_outer_points(1)=1;
while sum(flag_outer_points)~=length(flag_outer_points)
    dist = sum((repmat(ordered_points(end,:),size(outer_points,1),1) - outer_points).^2,2);
    dist(flag_outer_points==1) = max(dist) + 1;
    [dummy,I] = min(dist);
    if dummy > sqrt(delta_x^2+delta_y^2)*2;
        flag_outer_points(I)=1; continue;
    end
    ordered_points = [ordered_points; outer_points(I,:)]; flag_outer_points(I)=1;
end
outer_points = ordered_points;



