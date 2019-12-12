function [populations] = get_cluster_populations(idx)
% This functions returns a vector whose each entry is a number of observations in each cluster.
%
% Input:
% ------------
% - idx
%     a vector specifying division to clusters.
%
% Output:
% ------------
% - populations
%     a vector specifying the number of observations present in each of the k clusters.

%% get_cluster_populations()
% Find the number of clusters:
k = numel(unique(idx));

% Initialize the populations vector:
populations = zeros(k,1);

% Find the populations in each clusters:
for ii = 1:1:k
    populations(ii) = sum((idx==ii));
end

end
