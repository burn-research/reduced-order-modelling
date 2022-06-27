function [idx_Y] = idx_vqpca_based_on_centroids(X, idx_X, Y, n_eigs, cent_crit, scal_crit)

% Get the training data dimensions:
[n_obs_X, n_vars_X] = size(X);

% Get the target data dimensions:
[n_obs_Y, n_vars_Y] = size(Y);

% Center and scale the training data:
[scal_X, centers_X] = center(X, cent_crit);
[scal_X, scales_X] = scale(scal_X, X, scal_crit);

% Apply the same centering and scaling to the target data:
[scal_Y, ~, ~] = center(Y, cent_crit, centers_X);
[scal_Y, ~, ~] = scale(scal_Y, Y, scal_crit, scales_X);

% Compute the centroids based on the existing VQPCA clustering solution:
[C] = get_centroids(scal_X, idx_X);

% Compute the training data partitoning based on the existing VQPCA clustering solution:
[nz_X_k, ~, k] = get_partition(scal_X, idx_X);

% Initialization of cell arrays:
eigvec = cell(k, 1);

% Perform PCA in local clusters:
for j = 1:1:k

    [centered_scaled_local_X, ~] = center(nz_X_k{j}, cent_crit);
    [modes] = pca(centered_scaled_local_X, 'Centered', false, 'Algorithm', 'svd');
    eigvec{j} = modes(:,1:n_eigs);
end

% Initialize the reconstruction error matrix:
sq_rec_err = zeros(n_obs_Y, k);

% Reconstruct the target data from the low-dimensional representation of the training data:
for j = 1:1:k

    C_mat = repmat(C(j, :), n_obs_Y, 1);
    rec_err_os = (scal_Y - C_mat - (scal_Y - C_mat) * eigvec{j} * eigvec{j}');
    sq_rec_err(:, j) = sum(rec_err_os.^2, 2);

end

% Assign the observations to clusters based on the minimum reconstruction error:
[~, idx_Y] = min(sq_rec_err, [], 2);

% Degrade clusters if needed:
if max(idx_Y) ~= length(unique(idx_Y))
    idx_Y = degrade_clusters(idx_Y);
end
