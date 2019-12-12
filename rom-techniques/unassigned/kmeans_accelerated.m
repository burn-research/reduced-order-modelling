function [idx, clusters, clusters_idx] = kmeans_accelerated(Y, ktot, n_clust_in, varargin)
%% Description
%
%

%% Input
% Total number of clusters to find
if nargin < 2
    ktot = 2; 
end
% Number of clusters of one run
if nargin < 3 || isempty(n_clust_in)
    n_clust_in = 3;
end

%% Main
y_dim = size(Y, 2);
% Initialization
n_clust = 1;
clusters = {Y};
idx = ones(y_dim, 1);
% Split each cluster until max num of clusters is reached
cond = true;
while(cond)
    % Split each cluster
    idx_c = cell(n_clust, 1);
    parfor ii = 1 : n_clust
        n_obj = size(clusters{ii},2);
        if n_obj > n_clust_in
            [idx_c{ii}] = kmeans(clusters{ii}', n_clust_in, 'MaxIter', 2000);
        else
            idx_c{ii} = ones(n_obj, 1);
        end
    end
    % Re-set global idx
    offset = 0;
    idx_old = idx;
    for ii = 1 : n_clust
        I = (idx_old == ii);
        idx(I) = idx_c{ii} + offset;
        % Update offset
        offset = max(idx);
    end
    idx = removeHoles(idx, true);
    % Update clusters
    n_clust_old = n_clust;
    n_clust = max(idx); % Update number of clusters found
    clusters = cell(n_clust,1);
    clusters_idx = cell(n_clust,1);
    for ii = 1 : n_clust
        clusters{ii} = Y(:, idx == ii);
        clusters_idx{ii} = idx(idx == ii);
    end
    % Convergence checks
    if n_clust <= n_clust_old
        % Check the number of clusters has increased
        cond = false;
    elseif n_clust >= ktot
        % Check goal has been reached
        cond = false;
    end
end
end

