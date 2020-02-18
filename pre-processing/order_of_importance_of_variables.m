function [order] = order_of_importance_of_variables(X, metric)
% This function returns an order at which variables in a data set
% are important according to a metric specified.
%
% Example:
% ------------
% Suppose your data set is:
%
%      1   2    3  4   5
% X = [10, 100, 2, 56, 20]
% 
% Assigning importance based on maximum values, the returned order vector will be:
%
% order = [2, 4, 5, 1, 3]
% 
% So the most important variable is at index 2, the least important
% variable is at index 3.
%
% Input:
% ------------
% - X
%         raw data set.
%
% - metric
%         option for metric to judge the imporance.
%         By default the importance will be established based on the
%         maximum achieved observation.
%         Available metrics:
%
%         1     Max
%         2     Mean
%
% Output:
% ------------
% - order
%         vector specifying an idx of variables in the order of decreasing importance.

%% importance_order_of_variables()
[n_obs, n_vars] = size(X);

if ~exist('metric', 'var') || isempty(metric)
    metric = 1;
end

order = [];

if metric == 1

    maximum_values = max(X, [], 1);
    [~, order] = sort(maximum_values, 'descend');
    
elseif metric == 2
    
    mean_values = mean(X);
    [~, order] = sort(mean_values, 'descend');
    
end

end
