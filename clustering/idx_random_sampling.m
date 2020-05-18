function [idx] = idx_random_sampling(X, k)
% This function partitions the data into `k` clusters randomly. This
% function is mostly useful for benchmarking with other techniques.
%
% Input:
% ------------
% - X
%       the data set to partition.
%
% - k
%       number of bins (clusters) to partition the data set.
%
% Output:
% ------------
% - idx
%       a vector specifying division to clusters.

[n_obs, ~] = size(X);

k_values = 1:1:k;
idx = zeros(n_obs,1);

for j = 1:1:n_obs

    idx(j) = randsample(k_values, 1);
    
end

end
