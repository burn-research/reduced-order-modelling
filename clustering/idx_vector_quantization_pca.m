function [idx] = idx_vector_quantization_pca(X, n_eigs, k, cent_crit, scal_crit, idx_0)
% This function partitions the data into `k` clusters according to
% Vector Quantization Principal Component Analysis (VQPCA) algorithm.
%
% Input:
% ------------
% - X
%     the raw data set. Centering and scaling will be performed within this
%     function.
%
% - n_eigs
%     the number of eigenvector (Principal Components) to retain in the VQPCA algorithm.
%     When a low-rank approximation of the data is computed, it will be approximated
%     using only `n_eigs` Principal Components.
%
% - k
%     the number of clusters to partition the data to.
%
% - cent_crit
%     centering criterion (as per `center()` function).
%
% - scal_crit
%     scaling criterion (as per `scale()` function).
%
% - idx_0
%     a vector specifying division to clusters that will initialize the cluster centroids.
%     If not provided, a uniform initialization of cluster centroids is performed.
%
% Output:
% ------------
% - idx
%       a vector specifying division to clusters.

%% idx_vector_quantization()
% Get data dimensions:
[n_obs, n_vars] = size(X);

% Checks:
if exist('idx_0', 'var')
    if length(unique(idx_0)) ~= k
        error('The number of clusters in a user-provided idx is not equal to k.')
    end

    if length(idx_0) ~= n_obs
        error('The number of elements in a user-provided idx does not match the number of observations in the data set.')
    end
end

% Initialize parameters:
convergence = 0;
iter = 0;
iter_max = 600;
eps_rec = 1.0;
eps_rec_min = 1.0e-02;
a_tol = 1.0e-16;
r_tol = 1.0e-08;
cpu_start_time = cputime;
n_eigs_max = n_eigs;
eigvec = cell(k, 1);
gamma = cell(k, 1);

for j = 1:1:k
    eigvec{j} = eye(n_vars, n_eigs);
    gamma{j} = ones(1, n_vars);
end

% Center and scale the data:
[scal_X, ~] = center(X, cent_crit);
[scal_X, ~] = scale(scal_X, X, scal_crit);

% Initialization of cluster centroids:
if exist('idx_0', 'var')

    % If there is a user provided initial idx_0, find the initial centroids:
    [C] = get_centroids(scal_X, idx_0);

else

    % Initialize centroids automatically in a uniform way:
    C_int = linspace(1, n_obs, k+2);
    C = scal_X(round(C_int(2:k+1)), :);

end

tic;

% VQPCA algorithm:
while ((convergence == 0) && (iter < iter_max))

    fprintf('\nIteration n. %d, convergence %d \n', iter, convergence);

    % Initialize the reconstruction error matrix:
    sq_rec_err = zeros(n_obs, k);

    % Initialize the convergence of the cluster centroids:
    C_convergence = 0;

    % Initialize the convergence of the reconstruction error:
    eps_rec_convergence = 0;

    % Reconstruct the data from the low-dimensional representation, evaluate the mean squared reconstruction error:
    for j = 1:1:k

        D = diag(gamma{j});
        C_mat = repmat(C(j, :), n_obs, 1);

        rec_err_os = (scal_X - C_mat - (scal_X - C_mat) * D^-1 * eigvec{j} * eigvec{j}' * D);
        sq_rec_err(:, j) = sum(rec_err_os.^2, 2);

    end

    % Assign the observations to clusters based on the minimum reconstruction error:
    [rec_err_min, idx] = min(sq_rec_err, [], 2);
    rec_err_min_rel = rec_err_min;

    % Evaluate the global mean reconstruction error:
    eps_rec_new = mean(rec_err_min_rel);

    % Partition the data into clusters:
    [nz_X_k, nz_idx_clust, k] = get_partition(scal_X, idx, k);
    fprintf('\nClusters dimension \n');
    disp(nz_X_k);

    % Evaluate the relative recontruction errors in each cluster:
    rec_err_min_rel_k = cell(k, 1);

    for j = 1:1:k
        rec_err_min_rel_k{j} = rec_err_min_rel(nz_idx_clust{j}, 1);
    end

    % Evaluate the mean reconstruction error in each cluster:
    eps_rec_new_clust = zeros(k, 1);
    size_clust = zeros(k, 1);

    for j = 1:1:k
        eps_rec_new_clust(j) = mean(rec_err_min_rel_k{j});
        size_clust(j) = size(nz_X_k{j}, 1);
    end

    fprintf('\nGlobal mean recontruction error at iteration n. %d equal to %d \n', iter, eps_rec_new);
    fprintf('\nLocal mean recontruction error at iteration n. %d \n', iter);
    for j = 1:1:k
        fprintf('%d \t', eps_rec_new_clust(j));
    end
    fprintf('\n');

    % Find the new cluster centroids:
    C_new = zeros(k, n_vars);

    for j = 1 : k
        C_new(j, :) = mean(nz_X_k{j}, 1);
    end

    eps_rec_var = abs((eps_rec_new - eps_rec) / eps_rec_new);
    fprintf('\nReconstruction error variance equal to %d \n', eps_rec_var);

    if ((eps_rec_var < r_tol) && (eps_rec_new > eps_rec_min) && (n_eigs < n_eigs_max))
        n_eigs = n_eigs + 1;
        fprintf('\n Cluster %d dimension increased to %d \n', j,  n_eigs);
    end

    % Judge the convergence:
    if (eps_rec_var < r_tol)
        eps_rec_convergence = 1;
    end

    if (size(C) == size(C_new))

        C_var = abs((C_new - C) ./ (C_new + a_tol));

        if (C_var(:, :) < r_tol)
            C_convergence = 1;
        end

    end

    % If the convergence of centroids and reconstruction error is reached, the algorithm stops:
    if ((iter > 1) && (C_convergence == 1) && (eps_rec_convergence == 1))
        convergence = 1;
        fprintf('\nConvergence reached in %d iteration \n', iter);
    end

    % Update recontruction error and cluster centroids:
    C = C_new;
    eps_rec = eps_rec_new;

    % Initialization of cell arrays:
    eigvec = cell(k, 1);
    U_scores = cell(k, 1);
    gamma = cell(k, 1);

    % Perform PCA in local clusters:
    for j = 1:1:k
        [centered_scaled_local_X, ~] = center(nz_X_k{j}, cent_crit);
        [centered_scaled_local_X, gamma{j}] = scale(centered_scaled_local_X, nz_X_k{j}, 0); % don't scale in local clusters
        [modes] = pca(centered_scaled_local_X, 'Centered', false, 'Algorithm', 'svd');
        eigvec{j} = modes(:,1:n_eigs);
        u_scores{j} = centered_scaled_local_X * eigvec{j};
    end

    % Increment the iteration counter:
    iter = iter + 1;

end

% Evaluate total CPU time:
overall_cpu_time = cputime - cpu_start_time;

disp(['CPU time used: ', num2str(overall_cpu_time), ' seconds.']);
toc;

if (convergence == 0)
    fprintf('\nConvergence not reached in %d iterations \n', iter);
end

% Degrade clusters if needed:
if max(idx) ~= length(unique(idx))
    idx = degrade_clusters(idx);
end
