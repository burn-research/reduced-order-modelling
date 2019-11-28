function [idx, C] = init_idx(n, k, method, scal_X, varargin)
%% Description
% Initializes vector of cluster assignments.
% INPUT
% n: number of objects;
% k: number of clusters;
% [optional] method: 0 = uniform; 1 = kmeans; 2 = random;
%

%% Input
% Method
if ~exist('method', 'var') || isempty(method)
    method = 0;
end
% Data
if ~exist('scal_X', 'var') || isempty(scal_X)
    scal_X = [];
end
% Dummy check
if method == 1 && isempty(scal_X)
    error('Cannot use kmeans if no data are provided.');
end

%% Main
C = [];
if method == 0
    % Uniform
    I = round(linspace(1, n, k+1));
    idx = zeros(n, 1);
    for ii = 1 : length(I) - 1
        idx(I(ii):I(ii+1)) = ii;
    end
elseif method == 1
    % kmeans
    [idx, C] = kmeans(scal_X, k, 'MaxIter', 1e8);
else
    % Random
    idx = randi([1 k], 1, n);
end
% Make column vector 
idx = idx(:);

end


