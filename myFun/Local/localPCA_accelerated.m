function [idx, clusters, clusters_idx] = localPCA_accelerated(Y, ktot, n_modes, n_clust_in, idx_0, local_cent_crit, local_scal_crit, varargin)
%% Description
%
%

%% Input
% Total number of clusters to find
if nargin < 2
    ktot = 2; 
end
% Dimension of the cluster
if nargin < 3
    n_modes = 2;
end
% Number of clusters of one run
if nargin < 4 || isempty(n_clust_in)
    n_clust_in = 3;
end
% idx initialization
if nargin < 5
    idx_0 = 0;
end
% Local PCA centering and scaling criteria
if nargin < 6
    local_cent_crit = 1;
end
if nargin < 7
    local_scal_crit = 0;
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
    k = zeros(n_clust,1);
    for ii = 1 : n_clust
        clusters{ii} = clusters{ii}';
        evalc('[idx_c{ii}, c] = localPCA(clusters{ii}, n_modes, n_clust_in, idx_0, false, local_cent_crit, local_scal_crit);');
        clusters{ii} = clusters{ii}';
        k(ii) = length(c);
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

