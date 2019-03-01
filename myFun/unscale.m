function X = unscale(X, gamma)
%% Description:
% X: (n_obs x n_vars)
% gamma: array of scaling factors
%

%% Input
if isempty(gamma)
    return
end

%% Main
[~, cols] = size(X);
if length(gamma) ~= cols
    error('Array GAMMA (size: %i) must be as long as columns of X (size: %i).', length(gamma), size(X,2));
end
for ii = 1 : cols
    X(:,ii) = X(:,ii) * gamma(ii);
end

end


