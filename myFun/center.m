% This function provides a tool to center the data
% 
% INPUTS:
% 
% X   = The matrix of uncentered variables
% cent_crit         = centering criterion
%
% OUTPUT
%
% X0     = Matrix of centered data
% X_ave             = Matrix of means to be subtracted
function [X0, ave, X_ave] = center(X, cent_crit, val)
% Description:
% To get the principal components the data shas to be processed in order
% to center the column of X (variables with zero mean). Then 
% we will evaluate the eigenvalues of the matrix 1/(n-1) * X'X, which is
% the sample covariance matrix of X. Beside classical centering other 2
% possibilities are available:
% 1) The column of X are left uncentered
% 2) X is double centered, i.e. both rows and columns have zero mean
% 3) Centering by minimum value
% 4) Centering by maximum value

% Definition of parameters
[rows, columns] = size(X);
mean_var = mean(X, 1);
min_val = min(X, [], 1);
max_val = max(X, [], 1);
mean_obs = mean(X, 2);
grand_mean = mean(mean_var);
% Get user-supplied mean
if exist('val', 'var') && ~isempty(val)
    mean_var = val;
    cent_crit = 1;
end
% Relative tolerance
a_tol = 1e-08;
switch cent_crit
    case 0
        ave = zeros(1, columns);
        X_ave = repmat(ave, rows, 1);
    case 1
        ave = mean_var;
        X_ave = repmat(ave, rows, 1);
    case 2
        ave = mean_var;
        X_ave = repmat(ave, rows, 1) + repmat(mean_obs, 1, columns) - repmat(grand_mean, rows, columns);
    case 3 % Centering by minimum values
        ave = min_val;
        X_ave = repmat(ave, rows, 1);
    case 4 % Centering by maximum values
        ave = max_val;
        X_ave = repmat(ave, rows, 1);
    otherwise
        error('Unknown centering criterion');
end
% Subtract the mean
X0 = X - X_ave;
% We check now that each variable has zero mean
switch cent_crit
    case 1
        mean_vars_data = mean(X0, 1);
        for j = 1 : columns
            if (abs(mean_vars_data(j)) > a_tol)
               % disp(mean_vars_data);
               % error('The mean of each variable does not equal zero');
            end
        end
    case 2
        mean_vars_data = mean(X0, 1);
        mean_obs_data = mean(X0, 2);
        for j = 1 : columns
            if (abs(mean_vars_data(j)) > a_tol)
                error('The mean of each variable does not equal zero');
            end
        end
        for i = 1 : rows
            if (abs(mean_obs_data(i)) > a_tol)
                error('The mean of each row does not equal zero');
            end
        end
end
end



