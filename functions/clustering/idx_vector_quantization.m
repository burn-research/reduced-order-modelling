function [idx] = idx_vector_quantization(X, idx, k, is_remove, M, varargin)
% This function partitions the data into `k` clusters according to the Vector Quantization (VQ) algorithm.
%
% Input:
% ------------
% - X
%       the raw data set.
%
% - k
%       number of bins (clusters) to partition the data set.
%
% - idx_0
%       is the user-provided division to clusters that will be used as the initial centroids estimation.
%
% Output:
% ------------
% - idx
%       a vector specifying division to clusters.

%% Input
[rows, columns] = size(X);
if ~exist('is_remove', 'var') || isempty(is_remove)
    is_remove = false;
end
if ~exist('M', 'var') || isempty(M)
    M = columns;
end
%% Main
idx_clust = cell(k, 1);
n_points = zeros(k, 1);
% Evaluate clusters' populations
for j = 1 : k
    idx_clust{j} = find(idx == j); % Rows for cluster j
    n_points(j) = size(idx_clust{j}, 1);
    if is_remove
        % Remove clusters that do not have enough population
        if (n_points(j) < M)
            fprintf('\nNo points in the cluster n. %d, cluster removed \n', j);
        end
    end
end
% Find clusters with too little a population (nx_idx)
if ~is_remove
    nz_idx = find(n_points > -1); % Every cluster has a population > -1: no clusters will be removed
else
    nz_idx = find(n_points > M); % Only this clusters will be kept
end
k_new = size(nz_idx, 1); % New number of clusters
k = k_new;
nz_X_k = cell(k, 1);
nz_idx_clust = cell(k, 1);
% Re-assing data to clusters: the array "nz_idx" contains the indeces of
% the clusters to keep. "idx_clust" contains the row indeces corresponding
% to each cluster.
for j = 1 : k
    nz_X_k{j} = zeros(n_points(j), columns);
    nz_idx_clust{j} = idx_clust{nz_idx(j)};
    nz_X_k{j} = X(nz_idx_clust{j}, :);
end
% Optional output
idx = degrade_clusters(nz_idx);
end
