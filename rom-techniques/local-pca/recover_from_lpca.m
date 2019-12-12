function [X_app_lpca] = recover_from_lpca(idx, eigenvectors, scores, q, centroids, local_scalings, centerings, scalings)
% This function is used to recover the approximated data from Local Principal Component Analysis.
%
% Input:
% ------------
% - idx
%       a vector specifying division to clusters.
%
% - eigenvectors
%           a cell array of eigenvectors (Principal Components) in each cluster.
%
% - scores
%           a cell array of PC-scores in each cluster. It corresponds to
%           computing the equation scores = X*eigenvectors' in each cluster.
%
% - q
%           the number of Principal Components (PCs) to retain in the approximation.
%
% - centroids
%           a matrix of each cluster's centerings. Each row corresponds to each cluster and each
%           column corresponds to each variable.
%
% - local_scalings
%           a vector of local scalings that were applied to scale the data set.
%
% - centerings
%           a vector of centerings that were applied to center the data set.
%
% - scalings
%           a vector of scalings that were applied to scale the data set.
%
% Output:
% ------------
% - X_app_lpca
%           an approximated data set.

%% Input
% Get info
n_clust = length(eigenvectors);
n_samples = length(idx);
dim = size(eigenvectors{1}, 1);
% Handle optional inputs
if ~exist('q','var') || isempty(q)
    q = size(eigenvectors{1}, 2);
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
if ~exist('local_scalings', 'var') || isempty(local_scalings)
    local_scalings = ones(n_clust, dim);
elseif iscell(local_scalings)
    tmp = zeros(n_clust, dim);
    for ii = 1 : length(local_scalings)
        a = local_scalings{ii};
        tmp(ii,:) = a(:)';
    end
    local_scalings = tmp;
end
if ~exist('centerings','var') || isempty(centerings)
    centerings = [];
end
if ~exist('scalings','var') || isempty(scalings)
    scalings = [];
end
% Dummy check
if length(eigenvectors) ~= length(scores)
    error('No coherence found.');
end
%% Main
% Initialize matrix
X_app_lpca = zeros(n_samples, dim);
population = zeros(n_clust, 1);
cluster_size = zeros(n_clust, 1);
% Recover data from Local PCA
for ii = 1 : n_clust
    % Check app. order
    if q > size(scores{ii}, 2)
        qi = size(scores{ii}, 2);
    else
        qi = q;
    end
    % Check consistency of sizes
    population(ii) = size(scores{ii}, 1); % Population of i-th cluster
    cluster_size(ii) = sum(idx == ii); % Number of assignments to i-th cluster
    if population(ii) ~= cluster_size(ii)
        error('Cluster %i has %i assignments but a population of %i.',...
            ii, cluster_size(ii), population(ii));
    end
    % Recover data from one cluster
    C = scores{ii}(:,1:qi) * eigenvectors{ii}(:,1:qi)';
    C = unscale(C, local_scalings(ii,:)); % Unscale
    C = uncenter(C, centroids(ii,:)); % Uncenter
    X_app_lpca(idx == ii,:) = C;
end
% This data is still centered-scaled. We have recovered the centered-scaled
% data that was passed as input to the clustering/local PCA procedure.
if ~isempty(centerings) && ~isempty(scalings)
    X_app_lpca = unscale(X_app_lpca, scalings);
    X_app_lpca = uncenter(X_app_lpca, centerings);
end
end
