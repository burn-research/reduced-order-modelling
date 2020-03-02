function [clusters, clusters_idx, k_new] = get_partition(X, idx)
% This function partitions the dataset to clusters given by `idx` vector. It returns
% a cell array of data set divided into clusters and a cell array of indices of original observations
% diveded into clusters. If less than `n_vars` observations are assigned to a particular cluster,
% this cluster will be removed.
%
% Input:
% ------------
% - X
%     the raw data set.
%
% - idx
%     a vector specifying division to clusters.
%
% Output:
% ------------
% - clusters
%     a cell array containing the original data set X divided into clusters.
%     Each cell is linked to a particular cluster and contains all original data set observations in that cluster.
%
% - clusters_idx
%     a cell array containing the indices of the original data set X divided into clusters.
%     Each cell is linked to a particular cluster and contains indices of observations in that cluster.
%
% - k_new
%     updated number of clusters, different from the original k only if clusters have been removed.

%% get_partition()
[n_obs, n_vars] = size(X);

k = numel(unique(idx));

idx_clust = cell(k, 1);
n_points = zeros(k, 1);
clusters = cell(k, 1);
clusters_idx = cell(k, 1);

for j = 1:1:k

    idx_clust{j} = find(idx==j);
    n_points(j) = size(idx_clust{j}, 1);

    if (n_points(j) < n_vars)
        fprintf('\nToo few points in cluster k.%d, cluster will be removed.\n', j);
    end

end

nz_idx = find(n_points >= n_vars);
k_new = size(nz_idx, 1);

for j = 1:1:k_new

    % Assign observations to clusters:
    clusters{j} = zeros(n_points(j), n_vars);
    clusters_idx{j} = idx_clust{nz_idx(j)};
    clusters{j} = X(clusters_idx{j},:);

end
