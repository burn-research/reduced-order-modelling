function y_err = get_errors(obj, x_original, x_predicted, varargin)
% NOTES:
% The Mean Absolute Percentage Error (MAPE) is not a perfect error
% estimator when the predicted vector contains different variables. This
% function is used both for:
% - evaluate a mean error (among the variables) for each observation
% - evaluate a mean error (among observations) for each variable
% In the first case, averaging MAPEs among different variables may lead to
% the variable with the highest absolute values to dominate.
%

% Input
n_args = length(varargin);

% You should use COND to understand if it's observation errors for all the
% variables, and in that case do an average among the variables per
% observation 
if n_args > 0
    is_scale = varargin{1};
else
    is_scale = false;
end

% Sizes
[m, n] = size(x_original);

% Usually m is the number of variables and n is the number of observations,
% but this function works anyway and will return the error on the second
% dimension, thus per 'observations'

% Check dimensions are equal
if m ~= size(x_predicted,1) || n ~= size(x_predicted,2)
    error('You are trying to compare matrices with different sizes.');
end

% Evaluate errors
a_tol = 1e-16;
if is_scale && false  
    % We should try NRMSE:
    % error = sqrt( 1/s * sum_i^s[ (rom - full)^2 ] ) / (max_full - min_full)
    % We can use the variance though as denominator.
    % This formula usually works for the variable_errors, not the
    % observation errors.
    y_err = NRMSE(n, x_original, x_predicted);
else
    y_err = local_fun(n, x_original, x_predicted);
end

end

function y_err = local_fun(n, x_original, x_predicted, varargin)

n_args = length(varargin);
if n_args > 0
    is_scale = varargin{1};
else
    is_scale = true;
end

a_tol = 1e-16;
% Get errors
y_err = zeros(n, 1);
if is_scale 
    parfor i = 1 : n
        err = abs(x_predicted(:,i) - x_original(:,i));
        den = (norm(x_original(:,i)) + a_tol);
        y_err(i) = norm( err ) / den;
    end
else
    err = abs(x_predicted - x_original);
    y_err = sum(err, 2) / n;
end

y_err(y_err < a_tol) = 0;

end

function y_err = NRMSE(n, A, B)
% error = sqrt( 1/s * sum_i^s[ (rom - full)^2 ] ) / (max_full - min_full)

a_tol = 1e-16;
% Get errors
y_err = zeros(n, 1);
[~, mu, sigma] = zscore(A, 0, 2);
% A = center_scale(A, mu*0, sigma);
% B = center_scale(B, mu*0, sigma);
for i = 1 : n
    delta = (A(:,i) - B(:,i)) ./ (sigma + a_tol);
    y_err(i) = sqrt( (1/n) * sum(delta.^2) );
end

y_err(y_err < a_tol) = 0;
end

function Y = get_var(X, i, N)

if i == 1
    Y = X(i,:);
else
    Y = X(1+N*(i-1):i*N,:);
end

end








