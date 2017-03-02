function dist = weighted_tree_distance(weighted_adj_matrix)
% dist = graph_distance(weighted_adj_matrix) 
% adj_matrix can be weighted

adj_matrix = (abs(weighted_adj_matrix) + abs(weighted_adj_matrix') + eye(size(weighted_adj_matrix))) >0; 
weights = (abs(weighted_adj_matrix) + abs(weighted_adj_matrix'))/2;   
weights = weights - diag(diag(weights));

nNodes=size(adj_matrix,1);
shortest_hop  = zeros(nNodes);
shortest_dist = zeros(nNodes);
e = zeros(nNodes,1);
e(1)=1;

while sum(e==0)~=0  % e vector serves as a flag, see whether all nodes are "included/esamined"
    ind = find(sum(adj_matrix(e==1,:),1)'~=0 & e==0); % ind of all the nodes that are connected to one of the elements in e
    ind_new = ind(1); % take one of them
    ind_exist = find(adj_matrix(:,ind(1))==1 & e==1); % the one node in e that the new node connects to
    shortest_hop(ind_new,e==1) = shortest_hop(ind_exist,e==1)+1;
    shortest_hop(e==1,ind_new) = shortest_hop(e==1,ind_exist)+1;
    shortest_dist(ind_new,e==1) = shortest_dist(ind_exist,e==1)+weights(ind_exist,ind_new);
    shortest_dist(e==1,ind_new) = shortest_dist(e==1,ind_exist)+weights(ind_new,ind_exist);
    e(ind_new)=1;
end

shortest_hop = shortest_hop + eye(size(shortest_hop));
dist = shortest_dist;