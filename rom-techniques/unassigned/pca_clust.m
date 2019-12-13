function [idx, clusters, clusters_idx, T] = pca_clust(Y, n_pc, varargin)
%% Description
% This function clusters the variables of the input data-matrix around
% their Principal Components (PCs).
%
% INPUT
%       Y: data-matrix (obs x var)
%       n_pc: number of latent variables to use
% OUTPUT
%       idx: integer vector of cluster assignments
%       clusters: cell-array of clusters
%       clusters_idx: cell-arrays of row-indeces of the clusters
%       T: Principal Components
%
%% Input
if ~exist('Y', 'var') || isempty(Y)
    error('Please provide a data-matrix.');
end
if ~exist('n_pc', 'var') || isempty(n_pc)
    n_pc = 2; % Number of PCs used for the splitting
end
%% Main
% Get data dimension and do PCA
y_dim = size(Y,2);
[T, scores, eigvals] = pca(Y, 'Centered', false, 'Algorithm', 'svd');
% Check on n_pc
if n_pc > size(scores, 2)
    n_pc = size(scores, 2);
end
% Evaluate the affinities of each variable
T = abs(T');
T = T(1:n_pc, :);
% Assing each variable to a cluster: assignment depends on the affinities
% to the PCs
[~, idx] = max(T, [], 1);
idx = removeHoles(idx, true);
n_clust = max(idx); % Actual number of clusters found
clusters = cell(n_clust,1);
clusters_idx = cell(n_clust,1);
for ii = 1 : n_clust
    clusters{ii} = Y(:, idx == ii);
    clusters_idx{ii} = idx(idx == ii);
end
end
