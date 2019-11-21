function runPca_hiddenFunction(obj, varargin)

% Centered and scaled data (mean and scaling factors)
if ~obj.is_local 
    obj.runCenterScale();
% Do not run if this is a local cluster, the property centered_scaled_data
% has probably been already set.
end

% Run PCA
if ~obj.pca_manual
    % Use Matlab toolbox
    if obj.is_local && obj.is_mesh_variable
        [obj.pca_modes, obj.pca_scores, obj.pca_eigenvalues, obj.correlation_matrix] = local_fun_PCA(obj.training_data);
    else
        [obj.pca_modes, obj.pca_scores, obj.pca_eigenvalues, obj.correlation_matrix] = local_fun_PCA(obj.centered_scaled_data);
    end
else
    % Use cov() and eig()
    if obj.is_local && obj.is_mesh_variable
        [obj.pca_modes, obj.pca_scores, obj.pca_eigenvalues, obj.correlation_matrix] = manual_fun_PCA(obj.training_data);
    else
        [obj.pca_modes, obj.pca_scores, obj.pca_eigenvalues, obj.correlation_matrix] = manual_fun_PCA(obj.centered_scaled_data);
    end
end

% Set approximation order if not set, change it if greater than the size of
% the PCA modes matrix
if isempty(obj.pca_approximation_order)
    obj.pca_approximation_order = size(obj.pca_modes, 2);
elseif obj.pca_approximation_order > size(obj.pca_modes, 2)
    obj.pca_approximation_order = size(obj.pca_modes, 2);
end

end


function [modes, scores, eigenvalues, C] = local_fun_PCA(data, varargin)
%% Inputs
n_args = length(varargin);

% Dimensions
[n_var, n_samples] = size(data);

% Center-scale?
cen = false; 
if n_args > 0
    cen = varargin{1};
end

% Choose algorithm
if n_var > n_samples
    alg = 'svd';
else
    alg = 'svd';
end
if n_args > 1
    alg = varargin{2};
end

%% Correlation matrix
try
    if n_var > n_samples
        C = cov(data);
    else
        C = cov(data');
    end
catch ME
    C = 0;
end

%% Apply PCA (Using Matlab's default toolbox)
[modes, scores, eigenvalues] = pca(data', 'Centered',cen, 'Algorithm',alg);
scores = scores';

% If n_var > n_samples, we can find at most (n_samples - 1) PCA directions,
% but the pca() toolbox does not seem to 'know' that. Doing a test,
% reconstruction errors were higher when using n_samples PCA modes. This
% means that, when using the last PC, error was introduced: this can be
% true only if the last PC should be zero (and thus brining no further
% improvement on the reconstruction error) but is not (it is some very
% small number).
if n_var > n_samples
    eigenvalues(end) = 0;
    modes(:,end) = zeros(n_var,1);
    scores(end,:) = zeros(1, n_samples);
end

end


function [modes, scores, eigenvalues, C] = manual_fun_PCA(data, varargin)
% Dimensions
[n_var, n_samples] = size(data);

% Inputs
n_args = length(varargin);

% Classic PCA or Sirovich
a_tol = 1e-16;
if n_var > n_samples 
    % Sirovich
    C = cov(data); % Correlation among observations
    [modes, eigenvalues] = eig(C);
    [modes, eigenvalues] = sortem(modes, eigenvalues);
    eigenvalues = diag(eigenvalues);
    modes = data * modes / n_samples;
    scaling_factors = sqrt( diag(modes' * modes) );
    for i = 1 : size(modes, 2)
        modes(:,i) = modes(:,i) / (scaling_factors(i) + a_tol);
    end
else
    % Classic
    C = cov(data'); % Correlation among variables
    [modes, eigenvalues] = eig(C);
    [modes, eigenvalues] = sortem(modes, eigenvalues);
    eigenvalues = diag(eigenvalues);
end
% C = cov(data');
% [modes, eigenvalues] = pcacov(C);
scores = modes' * data;

end




