function pops = get_cluster_populations(idx)
n_clust = numel(unique(idx));
pops = zeros(n_clust, 1);
for ii = 1 : n_clust
    pops(ii) = sum( (idx == ii) );
end
end



