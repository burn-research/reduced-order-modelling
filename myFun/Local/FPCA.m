function [bin_data, idx_clust, idx] = FPCA(data, z, n_bins, min_z, max_z, z_stoich)
% Check the inputs
if ~isvector(z)
    error('The conditioning variable MUST be a vector');
end
if ~exist('data', 'var') || isempty(data)
    data = z;
end
if ~exist('min_z', 'var') || isempty(min_z)
    min_z = min(z);
end
if ~exist('max_z', 'var') || isempty(max_z)
    max_z = max(z);
end
if ~exist('z_stoich', 'var') || isempty(z_stoich)
    z_stoich = [];
end
% Number of intervals
n = n_bins + 1;
% Partition of the space
if n_bins == 1
    ints = linspace(min_z, max_z, n);
elseif isempty(z_stoich)
    ints = linspace(min_z, max_z, n);
else 
    ints_1 = linspace(min_z, z_stoich, ceil(n/2));
    ints_2 = linspace(z_stoich, max_z, ceil((n+1)/2));
    ints = [ints_1(1:ceil(n/2-1)) ints_2];
end
% Bin data matrices initialization
bin_data = cell(n_bins, 1);
idx_clust = cell(n_bins, 1);
idx = zeros(size(data,1), 1);
%  Partition
for bin = 1 : n_bins  
    idx_clust{bin} = find((z >= ints(bin)) & (z <= ints(bin+1)));
    bin_data{bin} = data(idx_clust{bin}, :);
    idx(idx_clust{bin}) = bin;
end
end

% F_stoich = 0.4375;% DNS CO/H2 Flame
% F_stoich = 0.351; % Flame F
% F_stoich = 0.0579; %JHC Flame



