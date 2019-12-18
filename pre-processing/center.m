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

%% center()
% Get dimensions:
[n_obs, n_vars] = size(uncentered_data);

a_tol = 1e-08;

% Checks:
if exist('user_supplied_centering', 'var') && ~isempty(user_supplied_centering)
    mean_var = user_supplied_centering;
    cent_crit = 1;
end

% Find the vector of scalings:
if ~exist('user_supplied_centering', 'var') || isempty(user_supplied_centering)
  switch cent_crit
    case 0
        centering_name_str = 'none';
        centerings = zeros(1, n_vars);
        centerings_matrix = repmat(centerings, n_obs, 1);
    case 1
        centering_name_str = 'mean';
        centerings = mean(uncentered_data, 1);
        centerings_matrix = repmat(centerings, n_obs, 1);
    case 2
        centering_name_str = 'grand_mean';
        centerings = mean(uncentered_data, 1);
        centerings_matrix = repmat(centerings, n_obs, 1) + repmat(mean(uncentered_data, 1), 1, n_vars) - repmat(mean(mean_var), n_obs, n_vars);
    case 3
        centering_name_str = 'min';
        centerings = min(uncentered_data, [], 1);
        centerings_matrix = repmat(centerings, n_obs, 1);
    case 4
        centering_name_str = 'max';
        centerings = max(uncentered_data, [], 1);
        centerings_matrix = repmat(centerings, n_obs, 1);
    otherwise
        error('Unknown centering criterion');
else
  % Use user supplied centering:
  centerings = user_supplied_centerings;
  centerings_matrix = repmat(centerings, n_obs, 1);
end

% Subtract the centerings:
centered_data = uncentered_data - centerings_matrix;

end
