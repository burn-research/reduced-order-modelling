function lpca_X = recoverLpca(idx, A, Z, q, centroids, local_gamma, ave, gamma)
%% Description

%% Input
% Get info
n_clust = length(A);
n_samples = length(idx);
dim = size(A{1}, 1);
% Handle optional inputs
if ~exist('q','var') || isempty(q) 
    q = size(A{1}, 2);
end
if ~exist('centroids', 'var') || isempty(centroids)
    centroids = zeros(n_clust, dim);
elseif iscell(centroids)
    tmp = zeros(n_clust, dim);
    for ii = 1 : length(centroids)
        a = centroids{ii};
        tmp(ii,:) = a(:)';
    end
    centroids = tmp;
end
if ~exist('local_gamma', 'var') || isempty(local_gamma)
    local_gamma = ones(n_clust, dim);
elseif iscell(local_gamma)
    tmp = zeros(n_clust, dim);
    for ii = 1 : length(local_gamma)
        a = local_gamma{ii};
        tmp(ii,:) = a(:)';
    end
    local_gamma = tmp;
end
if ~exist('ave','var') || isempty(ave) 
    ave = [];
end
if ~exist('gamma','var') || isempty(gamma) 
    gamma = [];
end
% Dummy check
if length(A) ~= length(Z)
    error('No coherence found.');
end
%% Main
% Initialize matrix
lpca_X = zeros(n_samples, dim);
population = zeros(n_clust, 1);
cluster_size = zeros(n_clust, 1);
% Recover data from Local PCA
for ii = 1 : n_clust
    % Check app. order
    if q > size(Z{ii}, 2)
        qi = size(Z{ii}, 2);
    else
        qi = q;
    end
    % Check consistency of sizes
    population(ii) = size(Z{ii}, 1); % Population of i-th cluster
    cluster_size(ii) = sum(idx == ii); % Number of assignments to i-th cluster
    if population(ii) ~= cluster_size(ii)
        error('Cluster %i has %i assignments but a population of %i.',...
            ii, cluster_size(ii), population(ii));
    end
    % Recover data from one cluster
    C = Z{ii}(:,1:qi) * A{ii}(:,1:qi)';
    C = unscale(C, local_gamma(ii,:)); % Unscale
    C = uncenter(C, centroids(ii,:)); % Uncenter
    lpca_X(idx == ii,:) = C;
end
% This data is still centered-scaled. We have recovered the centered-scaled
% data that was passed as input to the clustering/local PCA procedure.
if ~isempty(ave) && ~isempty(gamma)
    lpca_X = unscale(lpca_X, gamma);
    lpca_X = uncenter(lpca_X, ave);
end
end

