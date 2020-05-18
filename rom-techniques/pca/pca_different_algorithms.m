function [eigenvalues, eigenvectors] = pca_different_algorithms(X, cent_crit, scal_crit, pca_algorithm)

[centered_scaled_X, centerings] = center(X, cent_crit);
[centered_scaled_X, scalings] = scale(centered_scaled_X, X, scal_crit);

[n_obs, n_vars] = size(X);

% PCA with eigendecomposition of the covariance matrix:
if strcmp(pca_algorithm, 'cov-eig')
    
    S = 1/(n_obs - 1) .* centered_scaled_X' * centered_scaled_X;
    [eigenvectors, eigenvalues] = eig(S);

end

% PCA with SVD of the covariance matrix:
if strcmp(pca_algorithm, 'cov-svd')
    
    S = 1/(n_obs - 1) .* centered_scaled_X' * centered_scaled_X;
    [U, eigenvalues, eigenvectors] = svd(S);
    
end

% PCA with SVD of the data matrix:
if strcmp(pca_algorithm, 'X-svd')
    
    [U, eigenvalues, eigenvectors] = svd(centered_scaled_X, 'econ');
    
end

% PCA builtin Matlab function that does SVD:
if strcmp(pca_algorithm, 'builtin')
    
    [eigenvectors, ~, eigenvalues] = pca(centered_scaled_X, 'Centered', false);
    
end

end