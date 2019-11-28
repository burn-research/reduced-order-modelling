function [modes, eigenvalues] = kernelPCA(X, type, para)
%% Description
%   X: data matrix, each row is one observation, each column is one feature
%   type: type of kernel, can be 'simple', 'poly', or 'gaussian'
%   para: parameter for computing the 'poly' and 'gaussian' kernel, 
%       for 'simple' it will be ignored
%   modes: eigen-vector

%% Check input
if ~( strcmp(type,'simple') || strcmp(type,'poly') || ...
        strcmp(type,'gaussian') )
    Z = [];
    modes = [];
    fprintf(['\nError: Kernel type ' type ' is not supported. \n']);
    return;
end

%% Main
% Size
N = size(X,1);
% Kernel PCA
K0 = kernel(X, type, para); 
oneN = ones(N, N) / N;
K = K0 - oneN*K0 - K0*oneN + oneN*K0*oneN;
% Eigenvalue analysis
[modes, eigenvalues] = eig(K/N);
eigenvalues = diag(eigenvalues);
[~, IX] = sort(eigenvalues, 'descend');
modes = modes(:,IX);
eigenvalues = eigenvalues(IX);
% Normailization
norm_modes = sqrt(sum(modes.^2));
modes = modes ./ repmat(norm_modes, size(modes,1), 1);

end



