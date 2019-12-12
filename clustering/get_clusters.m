function [clusters] = get_clusters(X, idx)
% This functions returns a cell array of data set divided into clusters.
%
% Input:
% ------------
% - X
%     the raw data set.
%
% - idx
%     a vector specifying division to clusters.
%
% Output:
% ------------
% - clusters
%     a cell array containing the original data set X divided into clusters.
%     Each cell is linked to a particular cluster and contains all original data set observations in that cluster.

%% get_clusters()
% Checks:
[n_obs, n_vars] = size(X);
if n_obs ~= length(idx)
  error('The number of observations in the dataset X must match the number of elements in idx vector.')
end

% Find the number of clusters:
k = numel(unique(idx));

% Initialize the clusters cell array:
clusters = cell(k, 1);

% Divide the data set into clusters:
for ii = 1 : k
    clusters{ii} = X(idx == ii, :);
end

end
