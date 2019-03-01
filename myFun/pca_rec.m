function [X] = pca_rec(A, Z, c, g, n_pc, varargin)
%% Description
% Input
%       A: matrix of PCs (vars x n_pc)
%       Z: matrix of scores (obs x vars)
%       c, g: centers and scaling factors (vars x 1)
%       n_pc: number of PCs to keep
% Output
%       X: recovered data matrix (obs x vars)

%% Input
x_dim = size(A, 1); % Number of variables
% Data centering factors
if nargin < 3 || isempty(c)
    c = zeros(x_dim, 1);
end
% Data scaling factors
if nargin < 4 || isempty(g)
    g = ones(x_dim, 1);
end
% Number of PCs to keep for the reconstruction of the data
if nargin < 5 || isempty(n_pc)
    n_pc = size(A, 2);
end

%% Main
X = Z(:,1:n_pc) * A(:,1:n_pc)';
X = unscale(X, g);
X = uncenter(X, c);

end



