function shortest_hop = tree_shortest_hop(adj_matrix)
% shortest_hop = tree_shortest_hop(adj_matrix)

adj_matrix = (abs(adj_matrix) + abs(adj_matrix') + eye(size(adj_matrix))) >0;
nNodes=size(adj_matrix,1);
shortest_hop = zeros(nNodes);
e = zeros(nNodes,1);
e(1)=1;

while sum(e==0)~=0  % e vector serves as a flag, see whether all nodes are "included/esamined"
    ind = find(sum(adj_matrix(e==1,:),1)'~=0 & e==0); % ind of all the nodes that are connected to one of the elements in e
    ind_new = ind(1); % take one of them
    ind_exist = find(adj_matrix(:,ind(1))==1 & e==1); % the one node in e that the new node connects to
    shortest_hop(ind_new,e==1) = shortest_hop(ind_exist,e==1)+1;
    shortest_hop(e==1,ind_new) = shortest_hop(e==1,ind_exist)+1;
    e(ind_new)=1;
end

shortest_hop = shortest_hop + eye(size(shortest_hop));