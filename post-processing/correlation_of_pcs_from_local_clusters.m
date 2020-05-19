function [correlation] = correlation_of_pcs_from_local_clusters(X, variable, which_pc, scal_crit, idx, weighted)
% This function finds the average correlation between two variables within local clusters.
% By default it computes the average weighted by the number of observations in each cluster.
%
% Input:
% ------------
% - X
%           raw data set.
%
% - variable
%           a vector of variable to find the correlation of PCs with.
%
% - which_pc
%           a scalar number specifying which PC to look at in each cluster.
%
% - scal_crit
%           scaling criteria to scale the raw data set.
%
% - idx
%           a vector specifying division to clusters.
%
% - weighted
%           boolean specifying if the correlation should be weighted by the
%           number of observations in each cluster.
%           It is set to true by default.
%           If set to false, the correlation will be weighted by the number
%           of clusters.
%
% Output:
% ------------
% - correlation
%           a scalar representing the cluster-averaged correlation between
%           the Principal Component score requested and an input variable.

%% average_correlation()
% Checks:
if ~exist('weighted', 'var') || isempty(weighted)
    weighted = true;
end

k = max(idx);
correlation = 0;
n_obs = 0;

% Center and scale data:
[X_cs, ~] = center(X, 1);
[X_cs, ~] = scale(X_cs, X, scal_crit);

% Perform local PCA:
[clusters] = get_clusters(X_cs, idx);
[~, scores, ~, ~, ~] = lpca(clusters);

% Partion the variable to clusters according to idx:
[variable_clusters] = get_clusters(variable, idx);

for i = 1:1:k

    n_obs = n_obs + size(scores{i}, 1);
    cluster_size = size(scores{i}, 1);
    
    if weighted == true
        correlation = correlation + abs(corr(scores{i}(:,which_pc), variable_clusters{i}))*cluster_size;
    else
        correlation = correlation + abs(corr(scores{i}(:,which_pc), variable_clusters{i}));
    end
end

if weighted == true
    correlation = correlation/n_obs;
else
    correlation = correlation/k;
end

end