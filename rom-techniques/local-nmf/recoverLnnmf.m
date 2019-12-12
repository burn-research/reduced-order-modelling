function [Y_nnmf] = recoverLnnmf(idx, W, H, centroids, local_gamma, ave, gamma)
%% Description

%% Input
y_dim = size(H{1}, 1); % Get size from number of weights
if ~exist('centroids','var') || isempty(centroids) 
    centroids = zeros(length(W), y_dim);
end
if ~exist('local_gamma','var') || isempty(local_gamma) 
    local_gamma = ones(length(W), y_dim);
end
if ~exist('ave','var') || isempty(ave) 
    ave = [];
end
if ~exist('gamma','var') || isempty(gamma) 
    gamma = [];
end
% Checks
if length(W) ~= length(H)
    error('Cell-arrays S and W must have same length.');
end

%% Main
nc = length(W);
m = length(idx);
Y_nnmf = zeros(m, y_dim);
for ii = 1 : nc
    % Recover cluster
    C_rec = (W{ii} * H{ii}');
    C_rec = unscale(C_rec, local_gamma(ii,:));
    C_rec = uncenter(C_rec, centroids(ii,:));
    try
        Y_nnmf(idx == ii, :) = C_rec; % Assign 
    catch ME
        fprintf('For cluster %i: \n', ii);
        fprintf('Size(C_rec) = %i, %i \n', size(C_rec));
        fprintf('Size(Y_nnmf(idx == ii, :)) = %i, %i', size(Y_nnmf(idx == ii, :)));
        error(ME);
    end
end
% Recover data from NNMF
if ~isempty(gamma) && ~isempty(ave)
    Y_nnmf = unscale(Y_nnmf, gamma);
    Y_nnmf = uncenter(Y_nnmf, ave);
end
end
% Whitening
% Zw = T * Z;
% Independent components
% S = W * Z;
% Z_ica = T \ (W' * S); % [case r x n]
% Z_ica = (S * W') / T; % [case n x r]

