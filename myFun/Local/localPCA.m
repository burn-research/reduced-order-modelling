function [idx, nz_X_k, u_scores, eigvec, nz_idx_clust, eps_rec, convergence] = localPCA(scal_X, n_eigs, k, idx_0, iter_max, cent_crit, scal_crit, varargin)
%% Description:
% INPUT
% scal_X: (n_obs x n_vars) - Data to be clustered.
% n_eigs: integer - Number of PCs.
% k: integer - Number of clusters (to look for).
% idx_0: User-supplied centroids initialization
% stop_rule: bool - if TRUE, the function stops early.
% cent_crit: integer - centering criterion inside the clusters.
% scal_crit: integer - scaling criterion inside the clusters.
%
% OUTPUT
% idx: (n_obs x 1), integer - Cluster assignments.
% nz_X_k: cell array - nz_X_k{ii} contains the objects of cluster ii.
% u_scores: cell array - u_scores{ii} contains the PCA scores evaluated in
%           cluster ii.
% eigvec: cell array - eigvec{ii} contains the PCA modes found in cluster
%         ii.
% nz_idx_clust: cell array - nz_idx_clust{ii} contains the indeces of the
%               objects assigned to cluster ii (row indeces).
%

%% Input
% Starting guess for IDX (thus, centroids initialization)
if ~exist('idx_0', 'var') || isempty(idx_0)
    % Empty: kmeans - Scalar: uniform - Array: initial solotuion
    idx_0 = [];
end
% Stop rule
stop_rule = false;
% Centering and scaling info
if ~exist('cent_crit', 'var') || isempty(cent_crit)
    cent_crit = 1; % Inside the cluster
end
if ~exist('scal_crit', 'var') || isempty(scal_crit)
    scal_crit = 0; % Inside the cluster
end
% Maximum number of iterations
if ~exist('iter_max', 'var') || isempty(iter_max) || islogical(iter_max)...
        || ~isscalar(iter_max)
    iter_max = 1500;
end

%% Parameters
% Convergence indicators initialization
convergence = 0;    
iter = 0;                   
eps_rec = 1.0;      
eps_rec_min = 1.0e-02;
a_tol = 1.0e-16;    
r_tol = 1.0e-8;
% Dimensions
[rows, columns] = size(scal_X);
if rows > columns
    isremove = true;
else
    isremove = false;
end
% DUMMY: Check that k < rows
if k >= rows
    warning('We cannot look for %k clusters among %d objects.', k, rows);
    k = max([1, round(rows * .5) - 1]); % E.g.: 100 objects, look for 49 clusters
    fprintf('\nNumber of clusters changed to %d.\n', k);
end

%% 1) INITIALIZATION
% Centroids initialization
if ~isempty(idx_0)
    % User-supplied centroids initialization
    if size(idx_0, 1) == rows 
        % Get clusters
        nz_X_k = cell(k, 1); 
        for j = 1 : k
            nz_X_k{j} = scal_X(idx_0 == j, :);
        end
        % Get cluster centroids
        C = zeros(k, columns);        
        for j = 1 : k
            C(j, :) = mean(nz_X_k{j}, 1);
        end
    else
        % Ignore idx_0 and get the centroids directly (uniform)
        C_int = linspace(1, rows, k+2);
        C = scal_X(round(C_int(2:k+1)), :);
    end
    % Initialization of eigenvectors
    eigvec = cell(k,1); 
    gamma = cell(k,1);
    parfor j = 1 : k
        eigvec{j} = eye(columns, n_eigs);
        gamma{j} = ones(1, columns);
    end
else
    % KMEANS initialization
    fprintf('\nCentroids and cluster assignments initialization. \n');
    try
        [idx_0, C] = kmeans(scal_X, k, 'MaxIter', 1e8); % Use kmeans to initialize 
    catch ME
        [idx_0, C] = kmeans(scal_X, k, 'MaxIter', 1e8); % Use kmeans to initialize 
    end
    k = length(unique(idx_0)); % Be sure k clusters were found 
    % Initialize clusters
    nz_X_k = cell(k,1); 
    for j = 1 : k
        nz_X_k{j} = scal_X(idx_0 == j, :);
    end
    % Initialization of eigenvectors
    [eigvec, n_eig, gamma] = LPCA_here(nz_X_k, n_eigs, cent_crit, scal_crit);
end
% Check if the user wants to stop here
if isempty(stop_rule) || stop_rule
    idx = idx_0;
    uncentered_nz_X_k = {};
    return
end

%% STEPS 2-3-4
% 2) PARTITION
n_eigs_max = min(size(scal_X));
while (convergence == 0 && iter < iter_max) && (k ~= 1)
    C_convergence = 0; % Convergence on the centroids  
    eps_rec_convergence = 0; % Convergence on the recovered variance
    fprintf('\nIteration n. %d, convergence %d \n', iter, convergence);
    % Print time info
    fprintf('%s - %s \n', date, datestr(now, 'HH:MM:SS'));
    % Reconstruction errors
    sq_rec_err = zeros(rows, k);
    parfor j = 1 : k
        D = get_D(gamma{j}, eigvec{j});
        C_mat = repmat(C(j, :), rows, 1);     
        % Squared mean reconstruction error
        rec_err_os = (scal_X - C_mat - (scal_X - C_mat) * D^-1 * eigvec{j} * eigvec{j}' * D);
        sq_rec_err(:, j) = sum(rec_err_os.^2, 2); 
    end
    % Evalutate the recovered minimum error, so assign an object to the
    % cluster whose PCA-manifold works best for it
    [rec_err_min, idx] = min(sq_rec_err, [], 2);
    rec_err_min_rel = (rec_err_min);
    % Evaluate the global mean error
    eps_rec_new = mean(rec_err_min_rel);
    % Partition the data into clusters
    [nz_X_k, nz_idx_clust, k, ~] = partitionVQ(scal_X, idx, k, isremove);
    fprintf('\nThe current number of clusters is %d.\n', length(nz_X_k));
    % Evaluate the relative recontruction errors in each cluster
    rec_err_min_rel_k = cell(k, 1);
    for j = 1 : k
        rec_err_min_rel_k{j} = rec_err_min_rel(nz_idx_clust{j}, 1);
    end
    % Evaluate the mean error in each cluster
    eps_rec_new_clust = zeros(k, 1);
    size_clust = zeros(k, 1);
    for j = 1 : k
        eps_rec_new_clust(j) = mean(rec_err_min_rel_k{j});
        size_clust(j) = size(nz_X_k{j}, 1);
    end
    fprintf('\nGlobal mean recontruction error at iteration n. %d equal to %d', iter, eps_rec_new);
% 3-4) EVALUATE NEW CLUSTERS' CENTROIDS & PERFORM LOCAL PCA
    [eigvec, n_eig, gamma, u_scores, C_new] = LPCA_here(nz_X_k, n_eigs, cent_crit, scal_crit);
    eps_rec_var = abs((eps_rec_new  - eps_rec) / eps_rec_new);
    fprintf('\nReconstruction error variance equal to %d \n', eps_rec_var);
    try
        if ((eps_rec_var < r_tol) && (eps_rec_new > eps_rec_min) && (n_eigs < n_eigs_max)) 
            n_eigs = n_eigs + 1;
            fprintf('\n Clusters dimension increased to %d \n', n_eigs);
        end
    catch
    end
    % Judge convergence: clusters centroid and relative reconstruction
    % error
    if (eps_rec_var < r_tol)
        eps_rec_convergence = 1;
    end
    if (size(C) == size(C_new))
        C_var = abs((C_new - C) ./ (C_new + a_tol));
        if (C_var(:, :) < r_tol)
            C_convergence = 1;
        end
    end
    if ((iter > 1) && (C_convergence == 1) && (eps_rec_convergence == 1))
        convergence = 1;
        fprintf('\nConvergence reached in %d iteartion \n', iter);
    end
    % Update recontruction error and cluster centroids
    C = C_new;
    eps_rec = eps_rec_new;
    iter = iter + 1; 
end

%% OUTPUT
% Return the clusters
if k == 1
    idx = ones(rows,1);
    [nz_X_k, nz_idx_clust] = partitionVQ(scal_X, idx, k);
else
    [nz_X_k, nz_idx_clust] = partitionVQ(scal_X, idx, k, isremove);
    [eigvec, n_eig, gamma, u_scores] = LPCA_here(nz_X_k, n_eigs, cent_crit, scal_crit);
end
% Message about convergence
if (convergence == 0)
    fprintf('\nConvergence not reached in %d iterations \n', iter);
end
fprintf('\n\n');

end

% f LPCA
function [eigvec, n_eig, gamma, u_scores, centroids, latent] = LPCA_here(nz_X_k, n_eigs, cent_crit, scal_crit)
% Centering and scaling criterion
if nargin < 3
    cent_crit = 1;
end
if nargin < 4
    scal_crit = 0;
end
% Number of clusters
k = length(nz_X_k);
% Initialization of cell arrays
ydim = size(nz_X_k{1}, 2);
centroids = zeros(k, ydim);
eigvec = cell(k, 1);
u_scores = cell(k, 1);
n_eig = cell(k, 1);
gamma = cell(k, 1);
latent = cell(k, 1);
% Apply PCA in each cluster
parfor j = 1 : k
    % Center and scale, then do PCA
    [X, centroids(j,:)] = center(nz_X_k{j}, cent_crit);
    [X, gamma{j}] = scale(X, nz_X_k{j}, scal_crit);
    [modes, ~, eigvals] = pca(X, 'Centered', false, 'Algorithm', 'svd'); 
    % Check n_eigs does not exceed the found number of modes
    n_modes = n_eigs;
    if n_eigs > size(modes, 2)
        n_modes = size(modes, 2);
    end
    % Outputs (and gamma)
    n_eig{j} = n_eigs;
    eigvec{j} = modes(:, 1:n_modes);  
    u_scores{j} = X * eigvec{j};
end
end

function D = get_D(gamma, eigvec)
gamma(gamma < eps) = eps;
try
    D = spdiags(gamma, 0, size(eigvec',2), size(eigvec',2));
catch
   try
       D = diag(gamma);
   catch
       D = 1;
   end
end
end



