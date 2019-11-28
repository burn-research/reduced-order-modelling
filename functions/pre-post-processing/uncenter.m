function X = uncenter(X, c)
%% Description:
% X: (n_obs x n_vars)
% mean_vals: array of mean values
%

%% Input
if isempty(c)
    return
end

%% Main
[~, cols] = size(X);
if length(c) ~= cols
    error('Array C must be as long as columns of X.');
end
for ii = 1 : cols
    X(:,ii) = X(:,ii) + c(ii);
end

end


