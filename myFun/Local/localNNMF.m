function [idx, clusters, W, H, T, eps_rec, convergence, C, gamma] = localNNMF(X, nF, k, idx_0, iter_max, cent_crit, scal_crit, alg, options, replicates, varargin)
%% Description:
% INPUT
% X: (n_obs x n_vars) - Data to be clustered.
% nF: integer - Low-rank approximation order.
% k: integer - Number of clusters (to look for).
% idx_0: User-supplied centroids initialization
% iter_max: bool - if TRUE, the function stops early.
% cent_crit: integer - centering criterion inside the clusters.
% scal_crit: integer - scaling criterion inside the clusters.
%
% OUTPUT
% idx: (n_obs x 1), integer - Cluster assignments.
% clusters: cell array - clusters{ii} contains the objects of cluster ii.
% W: cell array - W{ii} contains the NNMF scores evaluated in
%           cluster ii.
% H: cell array - H{ii} contains the non-negative factors found in cluster ii.
% T: 
%
% NOTES
% This function might be useless if size(X, 2) == nF.
%

%% Input
% Starting guess for IDX (thus, centroids initialization)
if ~exist('idx_0', 'var') || isempty(idx_0)
    % Empty: kmeans - Scalar: uniform - Array: initial solotuion
    idx_0 = [];
end
% Centering and scaling info
if ~exist('cent_crit', 'var') || isempty(cent_crit)
    cent_crit = 0; % Inside the cluster (no centering, instead of minimum)
end
if ~exist('scal_crit', 'var') || isempty(scal_crit)
    scal_crit = 0; % Inside the cluster
end
% Maximum number of iterations
if ~exist('iter_max', 'var') || isempty(iter_max) || islogical(iter_max)...
        || ~isscalar(iter_max)
    iter_max = 1500;
end
% Options
if ~exist('options','var') || isempty(options)
    options = statset();
    options.Display = 'off';
    options.MaxIter = 300;
    options.TolFun = 1e-6;
end
% Algorithm
if ~exist('alg','var') || isempty(alg)
    alg = 'mult'; 
end
% Replicates
if ~exist('replicates','var') || isempty(replicates)
    replicates = 5; 
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
[rows, columns] = size(X);
is_remove = (rows > columns);
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
        clusters = cell(k, 1); 
        for j = 1 : k
            clusters{j} = X(idx_0 == j, :);
        end
        % Get cluster centroids
        C = zeros(k, columns);        
        for j = 1 : k
            C(j, :) = mean(clusters{j}, 1);
        end
    else
        % Ignore idx_0 and get the centroids directly
        C_int = linspace(1, rows, k+2);
        C = X(round(C_int(2:k+1)), :);
    end
    % Initialization of eigenvectors
    H = cell(k,1);
    T = cell(k,1);
    gamma = ones(k, columns);
    parfor j = 1 : k
        H{j} = eye(columns, nF);
        T{j} = eye(columns, columns);
    end
else
    % KMEANS initialization
    fprintf('\nCentroids and cluster assignments initialization. \n');
    try
        [idx_0, C] = kmeans(X, k, 'MaxIter', 1e8); % Use kmeans to initialize 
    catch ME
        [idx_0, C] = kmeans(X, k, 'MaxIter', 1e8); % Use kmeans to initialize 
    end
    k = length(unique(idx_0)); % Be sure k clusters were found 
    % Initialize clusters
    clusters = cell(k,1); 
    for j = 1 : k
        clusters{j} = X(idx_0 == j, :);
    end
    % Initialization of eigenvectors
    [H, nFs, gamma, W, T] = NNMF_here(clusters, nF, cent_crit, scal_crit, alg, options, replicates);
end

%% STEPS 2-3-4
% 2) PARTITION
nIC_max = min(size(X));
while (convergence == 0 && iter < iter_max) && (k ~= 1)
    C_convergence = 0; % Convergence on the centroids  
    eps_rec_convergence = 0; % Convergence on the recovered variance
    fprintf('\nIteration n. %d, convergence %d \n', iter, convergence);
    % Print time info
    fprintf('%s - %s \n', date, datestr(now, 'HH:MM:SS'));
    % Reconstruction errors
    sq_rec_err = zeros(rows, k);
    parfor j = 1 : k
        % Scaling and centering
        try
            D = spdiags(gamma(j,:) + a_tol, 0, columns, columns);
        catch
           try
               D = diag(gamma(j,:) + a_tol);
           catch
               D = 1;
           end
        end
        C_mat = repmat(C(j, :), rows, 1);     
        % Squared mean reconstruction error
        rec_err_os = (X - C_mat - (X - C_mat) * D^-1 * H{j} * H{j}' * D);
        sq_rec_err(:, j) = sum(rec_err_os.^2, 2);
    end
    % Evalutate the recovered minimum error, so assign an object to the
    % cluster whose PCA-manifold works best for it
    [rec_err_min, idx] = min(sq_rec_err, [], 2);
    rec_err_min_rel = (rec_err_min);
    % Evaluate the global mean error
    eps_rec_new = mean(rec_err_min_rel);
    % Partition the data into clusters
    [clusters, nz_idx_clust, k, ~] = partitionVQ(X, idx, k, is_remove);
    fprintf('The current number of clusters is %d.\n', length(clusters));
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
        size_clust(j) = size(clusters{j}, 1);
    end
    fprintf('Global mean recontruction error at iteration n. %d equal to %d \n', iter, eps_rec_new);
% 3-4) EVALUATE NEW CLUSTERS' CENTROIDS & PERFORM LOCAL ICA
    [H, nFs, gamma, W, T, C_new] = NNMF_here(clusters, nF, cent_crit, scal_crit, alg, options, replicates);
    eps_rec_var = abs((eps_rec_new  - eps_rec) / eps_rec_new);
    fprintf('Reconstruction error variance equal to %d \n', eps_rec_var);
    % Increase cluster dimension
%     try
%         if ((eps_rec_var < r_tol) && (eps_rec_new > eps_rec_min) && (nIC < nIC_max)) 
%             nIC = nIC + 1;
%             fprintf('\n Clusters dimension increased to %d \n', nIC);
%         end
%     catch
%     end
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
    % Iteration counter
    iter = iter + 1; 
end

%% OUTPUT
% Return the clusters
if k == 1
    idx = ones(rows,1);
    [clusters, nz_idx_clust] = partitionVQ(X, idx, k);
else
    [clusters, nz_idx_clust] = partitionVQ(X, idx, k, is_remove);
    [H, nFs, gamma, W, T] = NNMF_here(clusters, nF, cent_crit, scal_crit, alg, options, replicates);
end
% Message about convergence
if (convergence == 0)
    fprintf('Convergence not reached in %d iterations \n', iter);
end
fprintf('\n\n');

end

function [H, nFs, gamma, W, T, centroids] = NNMF_here(clusters, nF, cent_crit, scal_crit, alg, options, replicates)
%% Description

%% Parameters
% Centering and scaling criterion
if ~exist('cent_crit', 'var') || isempty(cent_crit)
    cent_crit = 3;
end
if ~exist('scal_crit', 'var') || isempty(scal_crit)
    scal_crit = 0;
end
% Options
if ~exist('options','var') || isempty(options)
    options = statset();
    options.Display = 'off';
    options.MaxIter = 400;
    options.TolFun = 1e-6;
end
% Algorithm
if ~exist('alg','var') || isempty(alg)
    alg = 'mult'; 
end
% Replicates
if ~exist('replicates','var') || isempty(replicates)
    replicates = 5; 
end
%% Main
% Number of clusters
k = length(clusters);
% Initialization of cell arrays
ydim = size(clusters{1}, 2);
centroids = zeros(k, ydim);
W = cell(k, 1);
H = cell(k, 1);
T = cell(k, 1);
nFs = zeros(k, 1);
gamma = ones(k, ydim);
% Apply NNMF in each cluster
parfor j = 1 : k
    % Center and scale, then do NNMF
    [X, centroids(j,:)] = center(clusters{j}, cent_crit);
    [X, gamma(j,:)] = scale(X, clusters{j}, scal_crit);
    if nF <= min(size(X))
        nFs(j) = nF;
    else
        nFs(j) = min(size(X));
    end
    fprintf('Size(cluster{%i}) = %i x %i; nF = %i\n', j, size(X,1), size(X,2), nFs(j));
    [scores, modes, T{j}] = nnmf(X, nFs(j), 'algorithm', alg, 'options', options, 'replicates', replicates); 
    W{j} = scores;
    H{j} = modes';
end
end

