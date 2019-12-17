function [eigenvectors, scores, eigenvalues, centroids, local_scales] = lpca(clusters, q, cent_crit, scal_crit)
% This function performs Principal Component Analysis in each of the local clusters.
%
% Input:
% ------------
% - clusters
%           a cell array with division of data to clusters. This is the output
%           `clusters` from the function `get_clusters()`.
%           It has the following structure:
%
%           {[data samples from cluster k_1]; [from cluster k_2]; ...; [from cluster k_k]}
%
% - q
%           the number of eigenvectors (PCs) to retain in PCA. The returned variables will
%           be trimmed to show only q components. If this parameter is not given,
%           all eigenvectors will be retained.
%
% - cent_crit
%           centering criteria as per function `center()`. Data samples from each cluster
%           will be centered with this centering criteria.
%
% - scal_crit
%           scaling criteria as per function `scale()`. Data samples from each cluster
%           will be scaled with this scaling criteria. If the scaling is
%           not provided, data samples within the clusters will not be
%           scaled.
%
% Output:
% ------------
% - eigenvectors
%           a cell array of eigenvectors (Principal Components) in each cluster.
%
% - scores
%           a cell array of PC-scores in each cluster. It corresponds to
%           computing the equation Z = X*A' in each cluster.
%
% - eigenvalues
%           a cell array of eigenvalues in each cluster.
%
% - centroids
%           a matrix of each cluster's centers. Each row corresponds to each cluster and each
%           column corresponds to each variable.
%
% - local_scales
%           a matrix of each cluster's scalings. Each row corresponds to each cluster and each
%           column corresponds to each variable.

%% lpca()
% Number of clusters:
k = length(clusters);

% Number of variables:
n_vars = size(clusters{1}, 2);

% Checks:
if ~exist('cent_crit', 'var') || isempty(cent_crit)
    cent_crit = 1;
end
if ~exist('scal_crit', 'var') || isempty(scal_crit)
    scal_crit = 0;
end
if ~exist('q', 'var') || isempty(q)
    q = n_vars;
end

% Initialize outputs:
eigenvectors = cell(k, 1);
scores = cell(k, 1);
eigenvalues = cell(k, 1);
centroids = zeros(k, n_vars);
local_scales = zeros(k, n_vars);

for j = 1:1:k
    
    % Center and scale data samples from this cluster:
    [X, centroids(j,:)] = center(clusters{j}, cent_crit);
    [X, local_scales(j,:)] = scale(X, clusters{j}, scal_crit);

    % Perform PCA on data samples from this cluster:
    [eigvec, sc, eigvals] = pca(X, 'Centered', false, 'Algorithm', 'svd');

    % Retain only q first modes:
    eigenvectors{j} = eigvec(:, 1:q);
    scores{j} = sc(:, 1:q);
    eigenvalues{j} = eigvals(1:q);

end

end
