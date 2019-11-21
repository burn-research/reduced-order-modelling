function [clusters, centroids] = get_clusters(X, idx, is_rows)
%% Description
% INPUT
% X: data-matrix (observations x variables).
% idx: vector of cluster assignments.
% [optional] is_rows: TRUE if idx is meant for the rows of X, FALSE
%                     otherwise.
% OUTPUT
% clusters: cell-array of clusters of X

%% Inputs
if nargin < 3
    is_rows = true;
end
%% Main
n_clust = numel(unique(idx)); % Number of clusters
clusters = cell(n_clust, 1); % Initialize cell-array
if is_rows
    % Clusters are along the rows of the data-matrix
    for ii = 1 : n_clust
        clusters{ii} = X(idx == ii, :);
    end
    centroids = get_centroids(X, idx);
else
    % Clusters are along the columns of the data-matrix
    for ii = 1 : n_clust
        clusters{ii} = X(:, idx == ii);
    end
end
end
