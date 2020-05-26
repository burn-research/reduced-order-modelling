function [correlation] = global_correlation(var_1, var_2)
% This function finds the global correlation between two vector variables.
%
% Input:
% ------------
% - var_1
%           a vector of the first variable.
%
% - var_2
%           a vector of the second variable.
%
% Output:
% ------------
% - correlation
%           the average correlation between two input variables.

%% global_correlation()
% Checks:
if length(var_1) ~= length(var_2)
    error('The number of elements is different in both variables.');
end

correlation = abs(corr(var_1, var_2));

end