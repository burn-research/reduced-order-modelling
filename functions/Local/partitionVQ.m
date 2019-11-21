function [nz_X_k, nz_idx_clust, k, idx] = partitionVQ(X, idx, k, is_remove, M, varargin)
%% Description
% nz_X_k: cell-array of clusters
% nz_idx_clust: cell-array of row indeces of the clusters
% k: integer, number of clusters
%

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
idx = removeHoles(nz_idx, true);
end

function [nz_X_k, nz_idx_clust, k, varargout] = partitionVQ_old(X, idx, varargin)
[n_samples, n_variables] = size(X); % n_samples; n_variables
threshold = n_variables;
% Remove cluster if n_samples < threshold
isremove = true; % Default
if ~isempty(varargin)
    isremove = varargin{1}; % User-provided input
end
if ~isa(isremove, 'logical')
    isremove = false; % User was trying to provide an input: interpret as false
end
k = numel(unique(idx)); % Number of clusters = number of unique integers in IDX
idx_clust = cell(k, 1); % Initialize
n_points = zeros(k, 1); % Initialize
for j = 1 : k
    idx_clust{j} = find(idx == j); % Rows of X belonging to cluster j
    n_points(j) = size(idx_clust{j}, 1); % Clusters' population
    if isremove
        % Remove clusters that do not have enough population
        if (n_points(j) < threshold)
            fprintf('\nNo points in the cluster n. %d, cluster removed \n', j);
        end
    end
end
if ~isremove
    nz_idx = find(n_points > -1); % Every cluster has a population > -1: no clusters will be removed
else
    nz_idx = find(n_points > threshold); % Only this clusters will be kept
end
k = size(nz_idx, 1); % New number of clusters
nz_X_k = cell(k, 1); % Initialize
nz_idx_clust = cell(k, 1); % Initialize 
for j = 1 : k
    nz_X_k{j} = zeros(n_points(j), n_variables); % Cluster j
    nz_idx_clust{j} = idx_clust{nz_idx(j)}; % New IDX and IDX_CLUST: in case the number of clusters changed 
    nz_X_k{j} = X(nz_idx_clust{j}, :); % Cluster j
end
% Optional output
if nargout > 3
    idx = removeHoles(nz_idx, true);
    varargout{1} = idx;
end
end

function [nz_X_k, nz_idx_clust, k] = partitionVQ_original(X, idx, k, varargin)
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
nz_idx = find(n_points >= columns);  
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

 
       