function covariance_matrix = covariance_squared_exponential(X, h, lambda)
% This function computes a covariance matrix from a Squared Exponential covariance kernel.
%
% Input:
% ------------
% - X
%           data matrix.
%
% - h
%           hyperparameter: scaling factor.
%
% - lambda
%           hyperparameter: kernel width.
%
% Output:
% ------------
% - covariance_matrix
%           covariance matrix computed from a Squared Exponential covariance kernel.

%% covariance_squared_exponential()
[n_obs, n_vars] = size(X);

covariance_matrix = zeros(n_vars, n_vars);

for i = 1:1:n_obs
    for j = 1:1:n_obs
        covariance_matrix(i,j) = kernel_squared_exponential(X(i,:), X(j,:), h, lambda);
    end
end

end

function covariance_kernel = kernel_squared_exponential(x1, x2, h, lambda)
% This function computes a Squared Exponential covariance kernel.

covariance_function = h^2 * exp(-(x1-x2)^2/(lambda^2));

end
