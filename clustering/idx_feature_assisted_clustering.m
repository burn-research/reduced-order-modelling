function [idx] = idx_feature_assisted_clustering(X, artificial_modes, cent_crit, scal_crit)
% This function partitions the data into `k` clusters according to
% the Feature Assisted Clustering technique. The observation is assigned to
% a cluster `k_i` if the absolute value of the maximum score for that observation
% is for the `mode_i`.
%
% Input:
% ------------
% - X
%     the raw data set.
%
% - artificial_modes
%     a matrix of artificial modes. Each column of this matrix is a mode.
%
% - cent_crit
%     centering criterion (as per `center()` function).
%
% - scal_crit
%     scaling criterion (as per `scale()` function).
%
% Output:
% ------------
% - idx
%       a vector specifying division to clusters.

%% feature_assisted_clustering()
[n_obs, n_vars] = size(X);
[n_weights, n_modes] = size(artificial_modes);

% Checks:
if n_vars ~= n_weights
    error('The number of variables in a data set is not equal to the number of weights in each mode.')
end

% Center and scale data:
[X_cent_scal, ~] = center(X, cent_crit);
[X_cent_scal, ~] = scale(X_cent_scal, X, scal_crit);

% Initialize the output vector:
idx = zeros(n_obs,1);

% Project the data set on the artificial modes:
scores = X_cent_scal*artificial_modes;

% Perform cluster assignments:
for i = 1:1:n_obs
    [~, idx(i)] = max(abs(scores(i,:)));
end

end
