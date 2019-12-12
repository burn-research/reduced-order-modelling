function [idx, clusters, clusters_idx, T] = pca_split(Y, n_pc, varargin)
if nargin < 2
    n_pc = 2; % Number of PC used for the splitting
end
y_dim = size(Y,2);
[~, scores, ~] = pca(Y, 'Centered', false, 'Algorithm', 'svd');
n_scores = size(scores,2);
if n_pc > n_scores
%     fprintf('\n n_pc = %i; n_scores = %i; \n', n_pc, n_scores);
    n_pc = n_scores;
end
M = size(scores, 1); % Number of observations
% Evaluate the affinities of each variable to the first N PCs
T = zeros(n_pc, y_dim);
for ii = 1 : n_pc
    temp = (1/M) * (Y' * scores(:,ii)).^2;
    temp = temp / (norm(scores(:,ii))^2 + eps); % Normalize
    T(ii,:) = temp(:)';
end
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


