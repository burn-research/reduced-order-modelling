function [H, n_eig, gamma, W, eigenvalues, centroids, T, eps_rec] = lnmf(nz_X_k, q, cent_crit, scal_crit, alg, replicates, options, idx, varargin)
%% Description

%% Inputs
% Centering and scaling criterion
if ~exist('cent_crit', 'var') || isempty(cent_crit)
    cent_crit = 3;
end
if ~exist('scal_crit', 'var') || isempty(scal_crit)
    scal_crit = 0;
end
% Algorithm
if ~exist('alg','var') || isempty(alg)
    alg = 'mult';
end
% Replicates
if ~exist('replicates','var') || isempty(replicates)
    replicates = 5; 
end
% Idx
if ~exist('idx', 'var') || isempty(idx)
    idx = [];
end
% Options
if ~exist('options','var') || isempty(options)
    options = statset();
    options.Display = 'off';
    options.MaxIter = 300;
    options.TolFun = 1e-6;
end
%% Main
a_tol = 1e-16;
% Number of clusters
k = length(nz_X_k);
n_vars = size(nz_X_k{1}, 2);
% Initialization of cell arrays
H = cell(k, 1);
W = cell(k, 1);
T = cell(k, 1);
n_eig = zeros(k, 1);
gamma = zeros(k, n_vars);
eigenvalues = cell(k, 1);
centroids = zeros(k, n_vars);
% Apply NNMF in each cluster
sq_rec_err = zeros(length(idx), 1);
eps_rec = [];
for j = 1 : k
    % Center and scale, then do NNMF
    [X, centroids(j,:)] = center(nz_X_k{j}, cent_crit);
    [X, gamma(j,:)] = scale(X, nz_X_k{j}, scal_crit);
    if q <= min(size(X))
        n_eig(j) = q;
    else
        n_eig(j) = min(size(X));
    end
    [scores, modes, T{j}] = nnmf(X, n_eig(j), 'algorithm', alg, 'options', options, 'replicates', replicates);
    W{j} = scores;
    H{j} = modes';
    % Rec err
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
        rec_err_os = (nz_X_k{j} - C_mat) - (nz_X_k{j} - C_mat) * D^-1 * H{j} * H{j}' * D; % Get error
        sq_rec_err(idx == j, :) = sum(rec_err_os.^2, 2);
    end
end
if ~isempty(idx)
    eps_rec = mean(sq_rec_err);
end
end



