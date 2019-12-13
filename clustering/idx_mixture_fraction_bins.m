function [idx] = idx_mixture_fraction_bins(Z, k, Z_stoich)
% This function partitions the data into `k` clusters (bins) according to the mixture fraction space `Z`.
% Bins intervals are found starting from the first one which divides the `Z` space
% into two parts at the stoichiometric mixture fraction `Z_stoich`.
%
%      min_Z                      Z_stoich                          max_Z
%        |-------------|-------------|----------|----------|----------|
%             bin 1          bin 2       bin 3       bin 4     bin 5
%
% Input:
% ------------
% - Z
%       mixture fraction vector.
%
% - k
%       number of bins (clusters) to partition the mixture fraction space.
%
% - Z_stoich
%       stoichiometric mixture fraction. If Z_stoich is not specified, equal length intervals are created.
%
%       Selected stoichiometric mixture fractions:
%
%       Z_stoich = 0.055      CH4/air
%       Z_stoich = 0.4375     DNS CO/H2 Flame
%       Z_stoich = 0.351      Sandia Flame F
%       Z_stoich = 0.0579     JHC Flame
%
% Output:
% ------------
% - idx
%       a vector specifying division to clusters.

%% idx_mixture_fraction_bins()
% Checks:
if ~isvector(Z)
    error('The conditioning variable must be a vector.');
end
if ~exist('Z_stoich', 'var') || isempty(Z_stoich)
    Z_stoich = [];
end

% Number of interval borders:
n_bins_borders = k + 1;

% Minimum and maximum mixture fraction:
min_Z = min(Z);
max_Z = max(Z);

% Partition the Z space:
if k == 1
    ints = linspace(min_Z, max_Z, n_bins_borders);
elseif isempty(Z_stoich)

    % If Z_stoich is not specified, equal length intervals are created:
    ints = linspace(min_Z, max_Z, n_bins_borders);
else

    % Z space lower than stoichiometric mixture fraction:
    ints_1 = linspace(min_Z, Z_stoich, ceil(n_bins_borders/2));

    % Z space higher than stoichiometric mixture fraction:
    ints_2 = linspace(Z_stoich, max_Z, ceil((n_bins_borders+1)/2));

    % Combine the two partitions:
    ints = [ints_1(1:ceil(n_bins_borders/2-1)) ints_2];
end

% Bin data matrices initialization:
idx_clust = cell(k, 1);
idx = zeros(size(Z, 1), 1);

% Create the cluster division vector:
for bin = 1:1:k

    idx_clust{bin} = find((Z >= ints(bin)) & (Z <= ints(bin+1)));
    idx(idx_clust{bin}) = bin;

end
end
