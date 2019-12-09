function [eigvec, u_scores, eigenvalues, centroids, gamma, eps_rec] = lpca(clusters, cent_crit, scal_crit, is_cpca, idx, cpca_options)
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
% - is_cpca
%           boolean specifying if Local Constrained PCA will be used.
%           If set to true, an idx must be provided in cpca_options.
%
% - idx
%           idx specifying division to clusters for Local Constrained PCA.
%
% - cpca_options
%           options for Local Constrained PCA.
%
% Output:
% ------------
% - eigvec
%           a cell array of eigenvectors (Principal Components) in each cluster.
%
% - n_eig
%           number of eigenvectors
%
% - gamma
%           a matrix of scalings. Each row corresponds to each cluster and each
%           column corresponds to each variable.
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
% - eps_rec
%           reconstruction errors when using Local Constrained PCA.
%

%% lpca()
% Checks:
if ~exist('cent_crit', 'var') || isempty(cent_crit)
    cent_crit = 1;
end
if ~exist('scal_crit', 'var') || isempty(scal_crit)
    scal_crit = 0;
end
if ~exist('is_cpca', 'var') || isempty(is_cpca)
    is_cpca = false;
end
if ~exist('idx', 'var') || isempty(idx)
    idx = [];
end
if nargin < 6
    cpca_options = [];
end
if is_cpca
    cpca_options.idx = idx;
end
if is_cpca && isempty(idx)
    error('For Local Constrained PCA, you must provide idx.');
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
gamma = zeros(k, n_vars);
sq_rec_err = zeros(length(idx), 1);
eps_rec = [];

% Tolerance for division:
a_tol = 1e-16;

for j = 1:1:k

    % Center and scale data samples from this cluster:
    [X, centroids(j,:)] = center(clusters{j}, cent_crit);
    [X, gamma(j,:)] = scale(X, clusters{j}, scal_crit);

    % Perform PCA on data samples from this cluster:
    [eigvec{j}, u_scores{j}, eigenvalues{j}] = pca(X, 'Centered', false, 'Algorithm', 'svd');

    % Reconstruction errors:
    if ~isempty(idx)
        try
            D = spdiags(gamma(j,:) + a_tol, 0, n_vars, n_vars);
        catch
           try
               D = diag(gamma(j,:) + a_tol);
           catch
               D = 1;
           end
        end

        C_mat = repmat(centroids(j,:), size(X, 1), 1);
        rec_err_os = (clusters{j} - C_mat) - (clusters{j} - C_mat) * D^-1 * eigvec{j} * eigvec{j}' * D; % Get error
        sq_rec_err(idx == j, :) = sum(rec_err_os.^2, 2);

        eps_rec = mean(sq_rec_err);

    end
end

end
