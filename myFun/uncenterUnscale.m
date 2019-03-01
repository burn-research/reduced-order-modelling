function X = uncenterUnscale(X, mean_column, scaling_factors, varargin)
% Deescription:
% X (n_var x n_obs)

% Relavant dimension
N = size(X, 1); 
% Diagonal matrix of scaling factors
D = spdiags(scaling_factors, 0, N, N);
% Unscale
X = D * X;
% Uncenter
M = repmat(mean_column, 1, size(X, 2));
% Decoded data
X = M + X;
end

