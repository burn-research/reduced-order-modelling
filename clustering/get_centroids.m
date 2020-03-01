function [centroids] = get_centroids(X, idx)
% This functions returns a matrix of centroids of every cluster.
% Centroids are computed as the mean of all the observations of a specific
% variable in a particular cluster.
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
% - centroids
%     a matrix of centroids of every variable in every cluster.
%     The rows represent a cluster and the columns represent the variable.
%
%            var_1     ...    var_n_vars
%
%       k_1 [                           ]
%       k_2 [                           ]
%       .   [         centroids         ]
%       .   [                           ]
%       .   [                           ]
%       k_k [                           ]
%
%     Each entry in the centroids matrix represents the mean of that
%     variable in that particular cluster.

%% get_centroids()
% Checks:
[n_obs, n_vars] = size(X);

if n_obs ~= length(idx)
  error('The number of observations in the dataset X must match the number of elements in idx vector.')
end

% Find the number of clusters:
k = numel(unique(idx));

% Initialize the centroids matrix:
centroids = zeros(k, n_vars);

% Find the centroids:
for jj = 1:1:k
    centroids(jj, :) = mean(X(idx == jj,:), 1);
end

end
