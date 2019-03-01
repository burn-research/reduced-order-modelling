function [u_scores, X_rec, r2_rec, sq_rec_err, varargout] = lpca_projection(X, idx, modes, centroids, neta, varargin)
%% Description:
% Projects matrix X on the Local PCA manifolds MODES. It is implied that
% any pre-processing operations on X has happened prior to the use of this
% function.
%
% Input:
% X (n_obs x n_vars): data matrix
% idx (n_obs x 1): vector of cluster assignments
% modes (cell, n_clusters x 1): Local PCs
% centroids (n_clusters x n_vars): cluster centroids
% neta (optional, integer): number of local PCs
%
% Output:
% u_scores (cell, n_clusters x 1): each cell contains the local scores
% X_rec (n_obs x n_vars): X matrix recovered from the Local PCA compression
% 

%% Input
if nargin < 5
    neta = [];
end
% Number of clusters
if iscell(modes)
    n_clust = length(modes);
else
    error('Input MODES must be a cell-array.');
end
% Check on CENTROIDS
if size(centroids, 1) ~= n_clust
    error('You have %i clusters but provided %i centroids.', n_clust, size(centroids, 1));
end
% Check on IDX
if max(idx) > n_clust
    error('Rows of X have been assigned to cluster %i, but you provided only %i clusters.', max(idx), n_clust);
end

%% Main
% Evaluate local scores and recover data matrix
X_rec = zeros(size(X));
u_scores = cell(n_clust,1);
sq_rec_err = zeros(size(X,1), 1);
for ii = 1 : n_clust
    I = (idx == ii); % Get rows of ii-th cluster
    n_rows = sum(I); % Number of elements in cluster ii
    % Subtract centroid ii to rows assigned to cluster ii
    C_mat = repmat(centroids(ii,:), n_rows, 1);
    X0 = X(I,:) - C_mat;
    if isempty(neta)
        temp = modes{ii};
    else
        % Check on NETA
        if neta > size(modes{ii}, 2)
            neta_temp = size(modes{ii},2);
        else
            neta_temp = neta;
        end
        temp = modes{ii}(:,1:neta_temp);
    end
    % Project on PCs and get local scores
    u_scores{ii} = X0 * temp;
    % Recover matrix
    X_rec(I,:) = u_scores{ii} * temp' + C_mat;
    % Get orthogal distances to local manifold
    rec_err_os = X0 - u_scores{ii} * temp';
    sq_rec_err(I) = sqrt(sum(rec_err_os.^2, 2));
end
% Get R2 values
try
    r2_rec = rSquared(X, X_rec);
catch
    r2_rec = [];
end

end


% rec_err_os = (scal_X - C_mat - (scal_X - C_mat) * D^-1 * eigvec{j} * eigvec{j}' * D);
% sq_rec_err(:, j) = sum(rec_err_os.^2, 2);





