function [eigvec, u_scores, eigenvalues, centroids, scales] = lpca(clusters, cent_crit, scal_crit)
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
% - cent_crit
%           centering criteria as per function `center()`. Data samples from each cluster
%           will be centered with this centering criteria.
%
% - scal_crit
%           scaling criteria as per function `scale()`. Data samples from each cluster
%           will be scaled with this scaling criteria.
%
% Output:
% ------------
% - eigvec
%           a cell array of eigenvectors (Principal Components) in each cluster.
%
% - u_scores
%           a cell array of PC-scores in each cluster. It corresponds to
%           computing the equation Z = X*A' in each cluster.
%
% - eigenvalues
%           a cell array of eigenvalues in each cluster.
%
% - centroids
%           a matrix of centers. Each row corresponds to each cluster and each
%           column corresponds to each variable.
%
% - scales
%           a matrix of scalings. Each row corresponds to each cluster and each
%           column corresponds to each variable.

%% lpca()
% Checks:
if ~exist('cent_crit', 'var') || isempty(cent_crit)
    cent_crit = 1;
end
if ~exist('scal_crit', 'var') || isempty(scal_crit)
    scal_crit = 0;
end

% Number of clusters:
k = length(clusters);

% Number of variables:
n_vars = size(clusters{1}, 2);

% Initialize outputs:
eigvec = cell(k, 1);
u_scores = cell(k, 1);
eigenvalues = cell(k, 1);
centroids = zeros(k, n_vars);
scales = zeros(k, n_vars);

for j = 1:1:k

    % Center and scale data samples from this cluster:
    [X, centroids(j,:)] = center(clusters{j}, cent_crit);
    [X, scales(j,:)] = scale(X, clusters{j}, scal_crit);

    % Perform PCA on data samples from this cluster:
    [eigvec{j}, u_scores{j}, eigenvalues{j}] = pca(X, 'Centered', false, 'Algorithm', 'svd');

end

end
