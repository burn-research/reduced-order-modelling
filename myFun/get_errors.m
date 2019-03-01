function y = get_errors(x_original, x_predicted, varargin)

% Sizes
[m, n] = size(x_original);

% Check dimensions are equal
if m ~= size(x_predicted,1) || n ~= size(x_predicted,2)
    error('You are trying to compare matrices with different sizes.');
end

% Normalization factor
norm_factor = zeros(1, n);
parfor i = 1 : n
    norm_factor(:,i) = norm(x_original(:,i)) + eps; % row vector
end

% Get relative errors: err_rel(i) = ||x1 - x2|| / ||x1||
err_rel = zeros(n, 1);
parfor i = 1 : n
    err_rel(i,1) = norm( x_predicted(:,i) - x_original(:,i) ) / norm_factor(i);
end
y = err_rel;

end