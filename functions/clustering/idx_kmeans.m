function [idx] = idx_kmeans(X, k)
% This function partitions the data into `k` clusters according to the K-Means algorithm.
%
% Input:
% ------------
% - X
%       the raw data set.
%
% - k
%       number of bins (clusters) to partition the data set.
%
% Output:
% ------------
% - idx
%       a vector specifying division to clusters.

idx = kmeans(X, k)

end
