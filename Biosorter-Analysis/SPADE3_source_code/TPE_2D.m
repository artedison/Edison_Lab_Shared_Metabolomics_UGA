function [positions] = TPE_2D(DIST, real_positions)
% function [positions] = TPE_2D(DIST)

figure(10);

% total number of points to be visualized
N = size(DIST,1);

% build an mst according to the dist matrix
[adj,adj2, cost_value] = mst_from_dist_matrix(DIST);

% figure out the order of the edges during the mst construction
adj = triu(adj,1);
edges = SPADE_find_matrix_big_element(adj,1);
weights = diag(DIST(edges(:,1),edges(:,2)));  % for i=1:size(edges,1),  weights(i,1) = dist(edges(i,1),edges(i,2)); end
[~,I] = sort(weights,'ascend');
edges = edges(I,:);
weights = weights(I);

% initialize the visulation 
positions = zeros(2,N);
idx = 1:N;  % cluster assignments

% grow the visualizatoin according to the mst construction process
for i=1:size(edges,1)
    cluster1_id = idx(edges(i,1)); cluster1_members = find(idx==cluster1_id); cluster1_size = length(cluster1_members); 
    cluster2_id = idx(edges(i,2)); cluster2_members = find(idx==cluster2_id); cluster2_size = length(cluster2_members); 
    if cluster1_size<cluster2_size  % make sure that cluster 1 is bigger
        cluster1_id = idx(edges(i,2)); cluster1_members = find(idx==cluster1_id); cluster1_size = length(cluster1_members); 
        cluster2_id = idx(edges(i,1)); cluster2_members = find(idx==cluster2_id); cluster2_size = length(cluster2_members); 
    end    
    subplot(1,1,1)
    subplot(3,3,2); plot(real_positions(1,:),real_positions(2,:),'g.',real_positions(1,cluster1_members),real_positions(2,cluster1_members),'ob',real_positions(1,cluster2_members),real_positions(2,cluster2_members),'^r');  axis2 = axis;
    for kk=[cluster1_members,cluster2_members]
        subplot(3,3,5); text(real_positions(1,kk),real_positions(2,kk),num2str(kk)); axis(axis2);
    end

    
    if cluster1_size==1 && cluster2_size==1
        new_cluster_id = max(idx)+1;
        positions(:,cluster2_members) = [weights(i);0]; % [DIST(cluster1_members,cluster2_members),0];
        positions(:,[cluster1_members,cluster2_members]) = positions(:,[cluster1_members,cluster2_members]) - repmat(mean(positions(:,[cluster1_members,cluster2_members]),2),1,length(cluster1_members)+length(cluster2_members));
        idx(cluster1_members)=new_cluster_id;
        idx(cluster2_members)=new_cluster_id;
    end
    if cluster1_size~=1 && cluster2_size==1
        [current_guess,bubble_contour,score] = set_plus_one_contour_initial_gradient_search(positions(:,cluster1_members), DIST(cluster1_members,cluster2_members), weights(i));
        new_cluster_id = max(idx)+1;
        positions(:,cluster2_members) = current_guess;
        positions(:,[cluster1_members,cluster2_members]) = positions(:,[cluster1_members,cluster2_members]) - repmat(mean(positions(:,[cluster1_members,cluster2_members]),2),1,length(cluster1_members)+length(cluster2_members));
        idx(cluster1_members)=new_cluster_id;
        idx(cluster2_members)=new_cluster_id;
    end
    
    if cluster1_size~=1 && cluster2_size~=1  % clusters 1 and 2 both have multiple members
        current_guess = set_plus_one_contour_initial_gradient_search(positions(:,cluster1_members), mean(DIST(cluster1_members,cluster2_members),2), weights(i));         % gradient_based_search
        cosine = current_guess(1)/norm(current_guess);
        sine = current_guess(2)/norm(current_guess);
        rotation_matrix = [cosine, sine; ...
                           -sine,  cosine];
        positions(:,cluster1_members) = rotation_matrix*positions(:,cluster1_members);
        current_guess1 = rotation_matrix*current_guess;
        
        current_guess = set_plus_one_contour_initial_gradient_search(positions(:,cluster2_members), mean(DIST(cluster2_members,cluster1_members),2), weights(i));         % gradient_based_search
        cosine = current_guess(1)/norm(current_guess);
        sine = current_guess(2)/norm(current_guess);
        rotation_pi = [cos(pi), sin(pi); ...
                       -sin(pi),cos(pi)];
        rotation_matrix = [cosine, sine; ...
                           -sine,  cosine];
        positions(:,cluster2_members) = rotation_pi*rotation_matrix*positions(:,cluster2_members);
        current_guess2 = rotation_pi*rotation_matrix*current_guess;

        F1 = (positions(1,cluster1_members)'*ones(1,length(cluster2_members)) - ones(length(cluster1_members),1)*positions(1,cluster2_members)).^2 + (positions(2,cluster1_members)'*ones(1,length(cluster2_members)) - ones(length(cluster1_members),1)*positions(2,cluster2_members)).^2;
        F2 = (positions(1,cluster1_members)'*ones(1,length(cluster2_members)) - ones(length(cluster1_members),1)*positions(1,cluster2_members)).^2 + (positions(2,cluster1_members)'*ones(1,length(cluster2_members)) + ones(length(cluster1_members),1)*positions(2,cluster2_members)).^2;
        if sum(sum(abs(sqrt(F1)-DIST(cluster1_members,cluster2_members))))>sum(sum(abs(sqrt(F2)-DIST(cluster1_members,cluster2_members))))
            positions(2,cluster2_members)  = -positions(2,cluster2_members);
        end
        
        x = max(abs(current_guess1(1)),abs(current_guess2(1)))*10;
        y = 0; 
        theta = 0;
        
        initial_guess = [x;y;theta];
        current_guess = set_plus_set_gradient_based_search(positions(:,cluster1_members), positions(:,cluster2_members), DIST(cluster1_members,cluster2_members), initial_guess, weights(i));
        
        new_cluster_id = max(idx)+1;
        rotation_matrix = [cos(current_guess(3)), -sin(current_guess(3)); ...
                           sin(current_guess(3)),  cos(current_guess(3))];
        positions(:,cluster2_members) = rotation_matrix*positions(:,cluster2_members) + repmat(current_guess(1:2),1,length(cluster2_members));
        positions(:,[cluster1_members,cluster2_members]) = positions(:,[cluster1_members,cluster2_members]) - repmat(mean(positions(:,[cluster1_members,cluster2_members]),2),1,length(cluster1_members)+length(cluster2_members));
        idx(cluster1_members)=new_cluster_id;
        idx(cluster2_members)=new_cluster_id;        
    end
    
    subplot(3,3,1); plot(positions(1,:),positions(2,:),'g.',positions(1,cluster1_members),positions(2,cluster1_members),'ob',positions(1,cluster2_members),positions(2,cluster2_members),'^r');  axis equal; axis1 = axis;
    subplot(3,1,2);
    for kk=[cluster1_members,cluster2_members]
        subplot(3,3,4); text(positions(1,kk),positions(2,kk),num2str(kk));           axis(axis1);
        subplot(3,3,5); text(real_positions(1,kk),real_positions(2,kk),num2str(kk)); axis(axis2);
    end
    drawnow
    
    F = (positions(1,cluster1_members)'*ones(1,length(cluster2_members)) - ones(length(cluster1_members),1)*positions(1,cluster2_members)).^2 + (positions(2,cluster1_members)'*ones(1,length(cluster2_members)) - ones(length(cluster1_members),1)*positions(2,cluster2_members)).^2;
    % [sqrt(F),DIST(cluster1_members,cluster2_members)]
    % [min(sqrt(F(:))), weights(i), min(sqrt(F(:)))-weights(i)]
    
%     if weights(i)*0.99>min(sqrt(F(:)))
%         1
%     end
    
    if weights(i)<min(sqrt(F(:)))
        weights(i:end) = weights(i:end)/weights(i)*min(sqrt(F(:)));
    end
    
    subplot(3,3,9)
    plot(pdist(positions(:,[cluster1_members,cluster2_members])'),squareform(DIST([cluster1_members,cluster2_members],[cluster1_members,cluster2_members])),'.')
    
end



function current_guess = set_plus_set_gradient_based_search(positions_a, positions_b, Dab, initial_guess, merge_dist)

Na = size(positions_a,2);
Nb = size(positions_b,2);
current_guess = initial_guess;
current_positions_b = [cos(current_guess(3)), -sin(current_guess(3)); sin(current_guess(3)),  cos(current_guess(3))] * positions_b + repmat(current_guess(1:2),1,Nb);

subplot(3,3,7); plot(positions_a(1,:),positions_a(2,:),'ob', current_positions_b(1,:),current_positions_b(2,:),'.'); 
for kk=1:size(current_positions_b,2), text(current_positions_b(1,kk),current_positions_b(2,kk),num2str(kk)); end
drawnow

lambda = 0.01;
v = 1.5;
while 1
    Y = Dab.^2;
    Y = Y(:);
    current_positions_b = [cos(current_guess(3)), -sin(current_guess(3)); sin(current_guess(3)),  cos(current_guess(3))] * positions_b + repmat(current_guess(1:2),1,Nb);
    F = (positions_a(1,:)'*ones(1,Nb) - ones(Na,1)*current_positions_b(1,:)).^2 + (positions_a(2,:)'*ones(1,Nb) - ones(Na,1)*current_positions_b(2,:)).^2;
    F = F(:);
    Jx = 2 * (current_guess(1) + ones(Na,1)*positions_b(1,:)*cos(current_guess(3)) + ones(Na,1)*positions_b(2,:)*sin(current_guess(3)) - positions_a(1,:)'*ones(1,Nb));
    Jy = 2 * (current_guess(2) - ones(Na,1)*positions_b(1,:)*sin(current_guess(3)) - ones(Na,1)*positions_b(2,:)*cos(current_guess(3)) - positions_a(2,:)'*ones(1,Nb));
    Jtheta = 2 * (current_guess(1) + ones(Na,1)*positions_b(1,:)*cos(current_guess(3)) + ones(Na,1)*positions_b(2,:)*sin(current_guess(3)) - positions_a(1,:)'*ones(1,Nb)) .* (-ones(Na,1)*positions_b(1,:)*sin(current_guess(3)) + ones(Na,1)*positions_b(2,:)*cos(current_guess(3))) + ...
             2 * (current_guess(2) - ones(Na,1)*positions_b(1,:)*sin(current_guess(3)) + ones(Na,1)*positions_b(2,:)*cos(current_guess(3)) - positions_a(2,:)'*ones(1,Nb)) .* (-ones(Na,1)*positions_b(1,:)*cos(current_guess(3)) - ones(Na,1)*positions_b(2,:)*sin(current_guess(3)));
    J = [Jx(:),Jy(:),Jtheta(:)];     
    current_obj = sum((Y-F).^2);
    % decide lambda
    delta1 = pinv(J'*J+lambda*(J'*J))*(J'*(Y-F));
    current_positions_b = [cos(current_guess(3)+delta1(3)), -sin(current_guess(3)+delta1(3)); sin(current_guess(3)+delta1(3)),  cos(current_guess(3)+delta1(3))] * positions_b + repmat(current_guess(1:2)+delta1(1:2),1,Nb);
    F = (positions_a(1,:)'*ones(1,Nb) - ones(Na,1)*current_positions_b(1,:)).^2 + (positions_a(2,:)'*ones(1,Nb) - ones(Na,1)*current_positions_b(2,:)).^2;
    F = F(:);
    obj1 = sum((Y-F).^2);
    delta2 = pinv(J'*J+(lambda/v)*(J'*J))*(J'*(Y-F));
    current_positions_b = [cos(current_guess(3)+delta2(3)), -sin(current_guess(3)+delta2(3)); sin(current_guess(3)+delta2(3)),  cos(current_guess(3)+delta2(3))] * positions_b + repmat(current_guess(1:2)+delta2(1:2),1,Nb);
    F = (positions_a(1,:)'*ones(1,Nb) - ones(Na,1)*current_positions_b(1,:)).^2 + (positions_a(2,:)'*ones(1,Nb) - ones(Na,1)*current_positions_b(2,:)).^2;
    F = F(:);
    obj2 = sum((Y-F).^2);
    while obj1>current_obj && obj2>current_obj
        lambda = lambda*v;
        delta1 = pinv(J'*J+lambda*(J'*J))*(J'*(Y-F));
        current_positions_b = [cos(current_guess(3)+delta1(3)), -sin(current_guess(3)+delta1(3)); sin(current_guess(3)+delta1(3)),  cos(current_guess(3)+delta1(3))] * positions_b + repmat(current_guess(1:2)+delta1(1:2),1,Nb);
        F = (positions_a(1,:)'*ones(1,Nb) - ones(Na,1)*current_positions_b(1,:)).^2 + (positions_a(2,:)'*ones(1,Nb) - ones(Na,1)*current_positions_b(2,:)).^2;
        F = F(:);
        obj1 = sum((Y-F).^2);
        delta2 = pinv(J'*J+(lambda/v)*(J'*J))*(J'*(Y-F));
        current_positions_b = [cos(current_guess(3)+delta2(3)), -sin(current_guess(3)+delta2(3)); sin(current_guess(3)+delta2(3)),  cos(current_guess(3)+delta2(3))] * positions_b + repmat(current_guess(1:2)+delta2(1:2),1,Nb);
        F = (positions_a(1,:)'*ones(1,Nb) - ones(Na,1)*current_positions_b(1,:)).^2 + (positions_a(2,:)'*ones(1,Nb) - ones(Na,1)*current_positions_b(2,:)).^2;
        F = F(:);
        obj2 = sum((Y-F).^2);
    end
    if obj2<current_obj
        lambda = lambda/v;
    end
    delta = pinv(J'*J+lambda*(J'*J))*(J'*(Y-F));
    current_positions_b = [cos(current_guess(3)+delta(3)), -sin(current_guess(3)+delta(3)); sin(current_guess(3)+delta(3)),  cos(current_guess(3)+delta(3))] * positions_b + repmat(current_guess(1:2)++delta(1:2),1,Nb);
    F = (positions_a(1,:)'*ones(1,Nb) - ones(Na,1)*current_positions_b(1,:)).^2 + (positions_a(2,:)'*ones(1,Nb) - ones(Na,1)*current_positions_b(2,:)).^2;
    F = F(:);
    if sum(F<merge_dist^2)~=0
        break;
    end
%     while sum(F<merge_dist^2)~=0 && norm(delta)>=merge_dist*(1e-10)
%         delta = delta/2;
%         current_positions_b = [cos(current_guess(3)+delta(3)), -sin(current_guess(3)+delta(3)); sin(current_guess(3)+delta(3)),  cos(current_guess(3)+delta(3))] * positions_b + repmat(current_guess(1:2)++delta(1:2),1,Nb);
%         F = (positions_a(1,:)'*ones(1,Nb) - ones(Na,1)*current_positions_b(1,:)).^2 + (positions_a(2,:)'*ones(1,Nb) - ones(Na,1)*current_positions_b(2,:)).^2;
%         F = F(:);
%     end
    current_guess = current_guess + delta;
    current_positions_b = [cos(current_guess(3)), -sin(current_guess(3)); sin(current_guess(3)),  cos(current_guess(3))] * positions_b + repmat(current_guess(1:2),1,Nb);
    sum((Y-F).^2)
    if norm(delta)<merge_dist*(1e-10)
        break;
    end
    subplot(3,3,7); plot(positions_a(1,:),positions_a(2,:),'ob',current_positions_b(1,:),current_positions_b(2,:),'.r'); 
    for kk=1:size(current_positions_b,2), text(current_positions_b(1,kk),current_positions_b(2,kk),num2str(kk)); end
    drawnow
end






function [current_guess,bubble_contour,score]= set_plus_one_contour_initial_gradient_search(positions_a, Dab, merge_dist)

bubble_contour = get_annotation_contour(positions_a, merge_dist);
dist_contour_to_a = sqrt((repmat(bubble_contour(1,:)',1,size(positions_a,2)) - repmat(positions_a(1,:),size(bubble_contour,2),1)).^2 + (repmat(bubble_contour(2,:)',1,size(positions_a,2)) - repmat(positions_a(2,:),size(bubble_contour,2),1)).^2);
score = sum(abs(dist_contour_to_a - repmat(Dab(:)',size(dist_contour_to_a,1),1)),2); 
[~,ind] = min(score);
initial_guess = bubble_contour(:,ind);
[~,ind] = min(sqrt(sum((repmat(initial_guess,1,size(positions_a,2))-positions_a).^2,1)));
initial_guess = (initial_guess-positions_a(:,ind))/norm(initial_guess-positions_a(:,ind))*merge_dist*1.001 + positions_a(:,ind);
current_guess = set_plus_one_gradient_based_search(positions_a, Dab, initial_guess, merge_dist);

subplot(3,3,7); plot(positions_a(1,:),positions_a(2,:),'o',bubble_contour(1,:),bubble_contour(2,:),initial_guess(1,:),initial_guess(2,:),'*', current_guess(1,:), current_guess(2,:),'^'); axis equal;
subplot(3,3,8); plot(score);







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



function current_guess = set_plus_one_gradient_based_search(positions_a, Dab, initial_guess, merge_dist)

Na = size(positions_a,2);
current_guess = initial_guess;
% plot(positions_a(1,:),positions_a(2,:),'ob',current_guess(1),current_guess(2),'^r'); 
% drawnow
lambda = 0.01;
v = 1.5;
while 1
    Y = Dab.^2;
    F = sum((positions_a - repmat(current_guess,1,Na)).^2,1)';
    J = 2 * (repmat(current_guess,1,Na)-positions_a)';
    current_obj = sum((Y-F).^2);
    % decide lambda
    delta1 = pinv(J'*J+lambda*(J'*J))*(J'*(Y-F));
    obj1 = sum((Y - sum((positions_a - repmat(current_guess+delta1,1,Na)).^2,1)').^2);
    delta2 = pinv(J'*J+(lambda/v)*(J'*J))*(J'*(Y-F));
    obj2 = sum((Y - sum((positions_a - repmat(current_guess+delta2,1,Na)).^2,1)').^2);
    while obj1>current_obj && obj2>current_obj
        lambda = lambda*v;
        delta1 = pinv(J'*J+lambda*(J'*J))*(J'*(Y-F));
        obj1 = sum((Y - sum((positions_a - repmat(current_guess+delta1,1,Na)).^2,1)').^2);
        delta2 = pinv(J'*J+(lambda/v)*(J'*J))*(J'*(Y-F));
        obj2 = sum((Y - sum((positions_a - repmat(current_guess+delta2,1,Na)).^2,1)').^2);
    end
    if obj2<current_obj
        lambda = lambda/v;
    end
    delta = pinv(J'*J+lambda*(J'*J))*(J'*(Y-F));
    if sum(sum((positions_a - repmat(current_guess+delta,1,Na)).^2,1)'<merge_dist^2)~=0
        break;
    end
    current_guess = current_guess + delta;
%     plot(positions_a(1,:),positions_a(2,:),'ob',current_guess(1),current_guess(2),'^r'); 
%     drawnow
    if norm(delta)<merge_dist*(1e-10)
        break;
    end
end

        



