function [idx_new, eps_rec_new] = lpca_classification(X, Y, idx, nPC, mu, sigma, cent_crit, scal_crit, varargin)
%% Description
%INPUT
% X = old dataset already clustered
% Y = new dataset for classification
% Additional info: the two datasets are supposed to be centered and scaled
% already and correctly.

%OUTPUT
% idx_new = the new idx column with the cluster distribution
% eps_rec_new = global reconstruction error

%% Checks
if size(X,2) ~= size(Y,2)
    error('Input X and Y must have same number of columns.');
end

%% Inputs
% Centering values
if ~exist('mu', 'var') || isempty(mu)
    mu = zeros(1, size(X,2));
elseif isscalar(mu) % User provided a criterion, not the values
    cent = mu;
    [~, mu] = center(X, cent);
end
% Scaling factors
if ~exist('sigma', 'var') || isempty(sigma)
    sigma = ones(1, size(X,2));
elseif isscalar(sigma) % User provided a criterion, not the factors
    scal = sigma;
    [~, sigma] = scale(X, center(X, [], mu), scal);
end
% Center scale
[X, ~] = center(X, [], mu);
[X, ~] = scale(X, [], [], sigma);
[Y, ~] = center(Y, [], mu);
[Y, ~] = scale(Y, [], [], sigma);
% Local centering criterion
if ~exist('cent_crit', 'var') || isempty(cent_crit)
    cent_crit = 1;
end
% Local scaling criterion
if ~exist('scal_crit', 'var') || isempty(scal_crit)
    scal_crit = 0;
end

%% Main
[rows, cols] = size(Y);
a_tol = 1e-16;
% Get non-centered clusters
clusters = get_clusters(X, idx); 
n_clust = length(clusters);
% Perform LPCA 
[A, ~, gamma, ~, ~, C] = lpca(clusters, nPC, cent_crit, scal_crit);
sq_rec_err = zeros(rows, n_clust);
for j = 1 : n_clust
    % Get scaling and centering matrices
    try
        tmp = size(A{j}', 2);
        D = spdiags(gamma{j} + a_tol, 0, tmp, tmp);
    catch
       try
           D = diag(gamma{j} + a_tol);
       catch
           D = 1;
       end
    end
    C_mat = repmat(C(j, :), rows, 1);     
    % Squared mean reconstruction error
    rec_err_os = (Y - C_mat - (Y - C_mat) * D^-1 * A{j} * A{j}' * D);
    sq_rec_err(:, j) = sum(rec_err_os.^2, 2);
end
% Evalutate the recovered minimum error, so assign an object to the
% cluster whose PCA-manifold works best for it
[eps_rec_new, idx_new] = min(sq_rec_err, [], 2);

end


