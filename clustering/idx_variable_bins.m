function [idx] = idx_variable_bins(variable, k)
% This function partitions data into `k` clusters by partitioning
% a `variable` vector into bins of equal sizes.
%
%   min_variable                                           max_variable
%        |----------|----------|----------|----------|----------|
%            bin 1      bin 2      bin 3      bin 4      bin 5
%
% Input:
% ------------
% - variable
%       variable vector.
%
% - k
%       number of bins (clusters) to partition the mixture fraction space.
%
% Output:
% ------------
% - idx
%       a vector specifying division to clusters.

%% idx_variable_bins()
var_min = min(variable);
var_max = max(variable);
bin_length = (var_max - var_min)/k;
bins_borders = [];

% Create bins borders:
for cluster = 0:1:k
    if cluster == k
        bins_borders = [bins_borders, var_max];
    else
        bins_borders = [bins_borders, (cluster*bin_length + var_min)];
    end
end

% Bin data matrices initialization:
idx_clust = cell(k, 1);
idx = zeros(size(variable, 1), 1);

% Create the cluster division vector:
for bin = 1:1:k

    if bin == k
        idx_clust{bin} = find((variable >= bins_borders(bin)) & (variable <= bins_borders(bin+1)));
    else
        idx_clust{bin} = find((variable >= bins_borders(bin)) & (variable < bins_borders(bin+1)));
    end
    idx(idx_clust{bin}) = bin;

end

end