function clusterss = find_envents_per_clusters(clusters,number_of_clusters)

for i = 1:number_of_clusters;
    clusterss(1,i) = length(find(clusters==i));
end

