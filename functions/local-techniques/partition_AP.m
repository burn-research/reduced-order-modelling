% This subroutine performs the partition the data into clusters
    
function [nz_X_k, nz_idx_clust, k] = partition_AP(X, idx, k)

[rows, columns] = size(X);

idx_clust = cell(k, 1);
n_points = zeros(k, 1); 

for j = 1 : k
    idx_clust{j} = find(idx == j);
    n_points(j) = size(idx_clust{j}, 1);
    if (n_points(j) < columns)
        fprintf('\nNo points in the cluster n. %d, cluster removed \n', j);
    end
end
nz_idx = find(n_points > columns);  
k_new = size(nz_idx, 1);
k = k_new;
nz_X_k = cell(k, 1);
nz_idx_clust = cell(k, 1);
for j = 1 : k
    nz_X_k{j} = zeros(n_points(j), columns);
    nz_idx_clust{j} = idx_clust{nz_idx(j)};
    nz_X_k{j} = X(nz_idx_clust{j}, :);
end
end