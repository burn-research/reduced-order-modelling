function [Y_ica] = recoverLica(idx, S, W, T, centroids, local_gamma, ave, gamma)
%% Description

%% Input
y_dim = size(W{1}, 1); % Get size from number of weights
if ~exist('T','var') || isempty(T) 
    T = [];
end
if ~exist('centroids','var') || isempty(centroids) 
    centroids = zeros(length(S), y_dim);
end
if ~exist('local_gamma','var') || isempty(local_gamma) 
    local_gamma = ones(length(S), y_dim);
end
if ~exist('ave','var') || isempty(ave) 
    ave = [];
end
if ~exist('gamma','var') || isempty(gamma) 
    gamma = [];
end
% Checks
if length(S) ~= length(W)
    error('Cell-arrays S and W must have same length.');
end

%% Main
nc = length(S);
m = length(idx);
Y_ica = zeros(m, y_dim);
for ii = 1 : nc
    % Recover cluster
    C_rec = (S{ii} * W{ii}'); % [case n x r]
    if ~isempty(T)
        C_rec = C_rec / T{ii};
    end
    C_rec = unscale(C_rec, local_gamma(ii,:));
    C_rec = uncenter(C_rec, centroids(ii,:));
    Y_ica(idx == ii, :) = C_rec; % Assign 
end
% Recover data from ICA
if ~isempty(gamma) && ~isempty(ave)
    Y_ica = unscale(Y_ica, gamma);
    Y_ica = uncenter(Y_ica, ave);
end
end
% Whitening
% Zw = T * Z;
% Independent components
% S = W * Z;
% Z_ica = T \ (W' * S); % [case r x n]
% Z_ica = (S * W') / T; % [case n x r]

