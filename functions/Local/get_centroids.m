function [centroids] = get_centroids(Y, idx)
%% Description
% Y: (obs x var)
%

%% Main
[~, n] = size(Y);
n_clust = numel(unique(idx));
centroids = zeros(n_clust, n);        
for jj = 1 : n_clust
    mask = (idx == jj);
    centroids(jj, :) = mean(Y(mask,:), 1);
end
end