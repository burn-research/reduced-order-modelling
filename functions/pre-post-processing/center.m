function [centered_data, centerings, centering_name_str] = center(uncentered_data, cent_crit, user_supplied_centering)
% This function centers the data set.
%
% Input:
% ------------
% - uncentered_data
%         an uncentered and unscaled data set.
%
% - cent_crit
%         selected centering criteria. It has to be an integer between 0 and 3.
%
%         Available scalings:
%
%         0   NONE, no centering is applied
%         1   MEAN, centering by mean value
%         2   GRAND-MEAN, both rows and columns have zero mean
%         3   MIN, centering by minimum value
%         4   MAX, centering by maximum value
%
% - user_supplied_centering
%         a user supplied vector of scalings that will be applied to the unscaled data set.
%
% Output:
% ------------
% - centered_data
%         a centered data set.
%
% - centerings
%         a vector of centerings that were applied to the uncentered data set.
%
% - centering_name_str
%         a string specifying the centering option.

% Definition of parameters
[n_obs, n_vars] = size(uncentered_data);
mean_var = mean(uncentered_data, 1);
min_val = min(uncentered_data, [], 1);
max_val = max(uncentered_data, [], 1);
mean_obs = mean(uncentered_data, 2);
grand_mean = mean(mean_var);
% Get user-supplied mean
if exist('user_supplied_centering', 'var') && ~isempty(user_supplied_centering)
    mean_var = user_supplied_centering;
    cent_crit = 1;
end
% Relative tolerance
a_tol = 1e-08;
switch cent_crit
    case 0
        centerings = zeros(1, n_vars);
        centerings_matrix = repmat(centerings, n_obs, 1);
    case 1
        centerings = mean_var;
        centerings_matrix = repmat(centerings, n_obs, 1);
    case 2
        centerings = mean_var;
        centerings_matrix = repmat(centerings, n_obs, 1) + repmat(mean_obs, 1, n_vars) - repmat(grand_mean, n_obs, n_vars);
    case 3 % Centering by minimum values
        centerings = min_val;
        centerings_matrix = repmat(centerings, n_obs, 1);
    case 4 % Centering by maximum values
        centerings = max_val;
        centerings_matrix = repmat(centerings, n_obs, 1);
    otherwise
        error('Unknown centering criterion');
end
% Subtract the mean
centered_data = uncentered_data - centerings_matrix;
% We check now that each variable has zero mean
switch cent_crit
    case 1
        mean_vars_data = mean(centered_data, 1);
        for j = 1 : n_vars
            if (abs(mean_vars_data(j)) > a_tol)
               % disp(mean_vars_data);
               % error('The mean of each variable does not equal zero');
            end
        end
    case 2
        mean_vars_data = mean(centered_data, 1);
        mean_obs_data = mean(centered_data, 2);
        for j = 1 : n_vars
            if (abs(mean_vars_data(j)) > a_tol)
                error('The mean of each variable does not equal zero');
            end
        end
        for i = 1 : n_obs
            if (abs(mean_obs_data(i)) > a_tol)
                error('The mean of each row does not equal zero');
            end
        end
end
end
