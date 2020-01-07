function [correlation] = average_correlation(var_1, var_2, weighted)
% This function finds the average correlation between two variables within local clusters.
% By default it computes the average weighted by the number of observations in each cluster.
%
% Input:
% ------------
% - var_1
%           a cell array of the first variable within local clusters.
%
% - var_2
%           a cell array of the second variable within local clusters.
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
%           the average correlation between two input variables.

%% average_correlation()
% Checks:
if size(var_1, 1) ~= size(var_2, 1)
    size(var_1, 1)
    size(var_2, 1)
    error('The number of clusters is different in both variables.');
end

if ~exist('weighted', 'var') || isempty(weighted)
    weighted = true;
end

k = size(var_1, 1);
correlation = 0;
n_obs = 0;

for i = 1:1:k
    
    % Additional check:
    if size(var_1{i}, 1) ~= size(var_2{i}, 1)
        error(['The number of observations in cluster ', num2str(i), ' is not the same in both variables.']);
    end
    
    n_obs = n_obs + size(var_1{i}, 1);
    cluster_size = size(var_1{i}, 1);
    
    if weighted == true
        correlation = correlation + abs(corr(var_1{i}, var_2{i}))*cluster_size;
    else
        correlation = correlation + abs(corr(var_1{i}, var_2{i}));
    end
end

if weighted == true
    correlation = correlation/n_obs;
else
    correlation = correlation/k;
end

end