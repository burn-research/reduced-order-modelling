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
%           It is equivalent to specifying the rank of the approximation.
%
% - centroids
%           a matrix of each cluster's centerings. Each row corresponds to each cluster and each
%           column corresponds to each variable.
%
% - local_scalings
%           a vector of local scalings that were applied to scale the data
%           set. Used when observations in each cluster were additionally
%           scaled apart from global data scaling.
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

%% recover_from_lpca()
% Get dimensions:
n_clust = length(eigenvectors);
n_samples = length(idx);
dim = size(eigenvectors{1}, 1);

% Checks:
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
    centerings = 0;
end
if ~exist('scalings','var') || isempty(scalings)
    scalings = 0;
end

% Initialize matrices:
X_app_lpca = zeros(n_samples, dim);
population = zeros(n_clust, 1);
cluster_size = zeros(n_clust, 1);

% Recover data from Local PCA:
for ii = 1:1:n_clust

    % Checks:
    population(ii) = size(scores{ii}, 1);
    cluster_size(ii) = sum(idx == ii);

    if q > size(scores{ii}, 2)
        qi = size(scores{ii}, 2);
    else
        qi = q;
    end

    if population(ii) ~= cluster_size(ii)
        error('Cluster %i has %i assignments but a population of %i.', ii, cluster_size(ii), population(ii));
    end

    % Recover data from a single cluster:
    C = scores{ii}(:,1:qi) * eigenvectors{ii}(:,1:qi)';

    % Unscale and uncenter observations in a single cluster:
    C = unscale(C, local_scalings(ii,:));
    C = uncenter(C, centroids(ii,:));

    % Append to the global data matrix:
    X_app_lpca(idx == ii,:) = C;

end

% Uncenter and unscale the data set globally:
if ~isempty(centerings) && ~isempty(scalings)

    X_app_lpca = unscale(X_app_lpca, scalings);
    X_app_lpca = uncenter(X_app_lpca, centerings);

end

end
