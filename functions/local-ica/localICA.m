function [idx, clusters, S, W, T, eps_rec, convergence, C, gamma] = localICA(X, nIC, k, idx_0, iter_max, cent_crit, scal_crit, varargin)
%% Description:
% INPUT
% X: (n_obs x n_vars) - Data to be clustered.
% nIC: integer - Number of PCs.
% k: integer - Number of clusters (to look for).
% idx_0: User-supplied centroids initialization
% iter_max: bool - if TRUE, the function stops early.
% cent_crit: integer - centering criterion inside the clusters.
% scal_crit: integer - scaling criterion inside the clusters.
%
% OUTPUT
% idx: (n_obs x 1), integer - Cluster assignments.
% clusters: cell array - clusters{ii} contains the objects of cluster ii.
% S: cell array - S{ii} contains the ICA scores evaluated in
%           cluster ii.
% W: cell array - W{ii} contains the ICs found in cluster ii.
% T: 
%
% NOTES
% This function might be useless if size(X, 2) == nIC.
% I think this can work only when the parameters cent_crit = 1.
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
[rows, columns] = size(X);
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
    W = cell(k,1);
    T = cell(k,1);
    gamma = ones(k, columns);
    parfor j = 1 : k
        W{j} = eye(columns, nIC);
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
    [W, nICs, gamma, S, T] = ICA_here(clusters, nIC, cent_crit, scal_crit);
end
% Check if the user wants to stop here
if isempty(stop_rule) || stop_rule
    idx = idx_0;
    uncentered_clusters = {};
    return
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
        % Squared mean reconstruction error [Z = (S * W') / T;] [S = Z * T * W;]
        X0 = center(X, [], C(j,:)); % Subtract centroid
        X_ica = scale(X0, X, [], gamma(j,:)); % Scale
        S_tmp = X_ica * T{j} * W{j}; % Get ICA scores
        X_ica = (S_tmp * W{j}') / T{j}; % Reconstruct
        X_ica = unscale(X_ica, gamma(j,:)); % Unscale
        rec_err_os = X0 - X_ica; % Get error
        sq_rec_err(:, j) = sum(rec_err_os.^2, 2); 
    end
    % Evalutate the recovered minimum error, so assign an object to the
    % cluster whose PCA-manifold works best for it
    [rec_err_min, idx] = min(sq_rec_err, [], 2);
    rec_err_min_rel = (rec_err_min);
    % Evaluate the global mean error
    eps_rec_new = mean(rec_err_min_rel);
    % Partition the data into clusters
    [clusters, nz_idx_clust, k, ~] = partitionVQ(X, idx, k, isremove);
    fprintf('\nThe current number of clusters is %d.\n', length(clusters));
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
    fprintf('\nGlobal mean recontruction error at iteration n. %d equal to %d', iter, eps_rec_new);
% 3-4) EVALUATE NEW CLUSTERS' CENTROIDS & PERFORM LOCAL ICA
    [W, nICs, gamma, S, T, C_new] = ICA_here(clusters, nIC, cent_crit, scal_crit);
    eps_rec_var = abs((eps_rec_new  - eps_rec) / eps_rec_new);
    fprintf('\nReconstruction error variance equal to %d \n', eps_rec_var);
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
    [clusters, nz_idx_clust] = partitionVQ(X, idx, k, isremove);
    [W, nICs, gamma, S, T] = ICA_here(clusters, nIC, cent_crit, scal_crit);
end
% Message about convergence
if (convergence == 0)
    fprintf('\nConvergence not reached in %d iterations \n', iter);
end
fprintf('\n\n');

end

function [W, nICs, gamma, S, T, centroids] = ICA_here(clusters, nIC, cent_crit, scal_crit, is_white, flag, type, TOL, MAX_ITERS)
%% Parameters
% Centering and scaling criterion
if ~exist('cent_crit', 'var') || isempty(cent_crit)
    cent_crit = 1;
end
if ~exist('scal_crit', 'var') || isempty(scal_crit)
    scal_crit = 0;
end
if ~exist('is_white','var') || isempty(is_white)
    % Default is_white
    is_white = true;
end
if ~exist('flag','var') || isempty(flag)
    % Default display flag
    flag = false;
end
if ~exist('type','var') || isempty(type)
    % Default type
    type = 'negentropy';
end
if ~exist('TOL','var') || isempty(TOL)
    % Default TOL
    TOL = 1e-7; % Convergence criteria
end
if ~exist('MAX_ITERS','var') || isempty(MAX_ITERS)
    % Default MAX_ITERS
    MAX_ITERS = 600; % Max iterations
end
%% Main
% Number of clusters
k = length(clusters);
% Initialization of cell arrays
ydim = size(clusters{1}, 2);
centroids = zeros(k, ydim);
W = cell(k, 1);
S = cell(k, 1);
T = cell(k, 1);
nICs = cell(k, 1);
gamma = ones(k, ydim);
% Apply PCA in each cluster
parfor j = 1 : k
    % Center and scale, then do ICA
    [X, centroids(j,:)] = center(clusters{j}, cent_crit);
    [X, gamma(j,:)] = scale(X, clusters{j}, scal_crit);
    [S{j}, W{j}, T{j}] = fastICA(X, nIC, is_white, type, flag, TOL, MAX_ITERS);
end
end

