function [scaled_data, scalings, scaling_name_str] = scale(unscaled_data, uncentered_data, scal_crit, user_supplied_scaling)
% This function scales the data set.
%
% Input:
% ------------
% - unscaled_data
%         an unscaled data set.
%
% - uncentered_data
%         an uncentered and unscaled data set.
%
% - scal_crit
%         selected scaling criteria. It has to be an integer between 0 and 9.
%
%         Available scalings:
%
%         0   NONE, no scaling is applied
%         1   AUTO (STD), each variable is normalized by its standard
%             deviation
%         2   RANGE, each variable is normalized by its range
%         3   PARETO, each variable is normalized by the square root of its standard
%             deviation
%         4   VAST, each variable is normalized by the standard deviation and
%             coefficient of variation
%         5   LEVEL, each variable is normalized by the mean of the data
%         6   MAX, each variable is normalized by its maximum value
%         7   KURTOSIS-1
%         8   KURTOSIS-2
%         9   KURTOSIS-3
%         10  NORM, each variable is normalized by its norm
%
% - user_supplied_scaling
%         a user supplied vector of scalings that will be applied to the unscaled data set.
%
% Output:
% ------------
% - scaled_data
%         a scaled data set.
%
% - scalings
%         a vector of scalings that were applied to the unscaled data set.
%
% - scaling_name_str
%         a string specifying the scaling option.

%% scale()
% Get dimensions:
[n_obs, n_vars] = size(unscaled_data);

a_tol = 1e-16;

% Checks:
if n_obs ~= size(uncentered_data, 1) || n_vars ~= size(unscaled_data, 2)
  error('Dimensions of the uncentered data must be the same as the dimensions of the unscaled data.')
end

if exist('user_supplied_scaling', 'var')
  if length(user_supplied_scaling) ~= n_vars
    error('The size of the user supplied scalings vector must match the number of variables in the data set.')
  end
end

% Find the vector of scalings:
if ~exist('user_supplied_scaling', 'var') || isempty(user_supplied_scaling)
    switch scal_crit
        case 0
            % NONE
            scalings = ones(1, n_vars);
            scaled_data = unscaled_data;
            return;
        case 1
            % AUTO
            scalings = std(uncentered_data, 1, 1);
            scaling_name_str = 'auto';
        case 2
            % RANGE
            scalings = max(abs(uncentered_data), [], 1) - min(uncentered_data, [], 1);
            scaling_name_str = 'range';
        case 3
            % PARETO
            scalings = std(uncentered_data, 1, 1).^0.5;
            scaling_name_str = 'pareto';
        case 4
            % VAST
            scalings = std(uncentered_data, 1, 1).^2 ./ (mean(uncentered_data, 1) + a_tol);
            scaling_name_str = 'vast';
        case 5
            % LEVEL
            scalings = mean(uncentered_data, 1);
            scaling_name_str = 'level';
        case 6
            % MAX
            scalings = max(abs(uncentered_data), [], 1);
            scaling_name_str = 'max';
        case 7
            % KURTOSIS-1
            kurt_value = kurtosis(uncentered_data);
            std_value = std(uncentered_data, 1, 1);
            mean_value = mean(uncentered_data, 1);
            vast_2 = std_value.^2 .* kurt_value.^2 ./ (mean_value + a_tol);
            scalings = vast_2;
            scaling_name_str = 'std_kurt_1';
        case 8
            % KURTOSIS-2
            kurt_value = kurtosis(uncentered_data);
            vast_3 = std_value.^2 .* kurt_value.^2 ./ max_value;
            scalings = vast_3;
            scaling_name_str = 'std_kurt_2';
        case 9
            % KURTOSIS-3
            kurt_value = kurtosis(uncentered_data);
            vast_4 = std_value.^2 .* kurt_value.^2 ./ range_value;
            scalings = vast_4;
            scaling_name_str = 'std_kurt_3';
        case 10
            % NORM
            scalings = vecnorm(unscaled_data);
            scaling_name_str = 'norm';
        otherwise
            error('Unknown scaling criterion');
    end
else
    % Use user supplied scaling:
    scalings = user_supplied_scaling;
    scaling_name_str = 'custom';
end

% Scale the data set:
scaled_data = zeros(n_obs, n_vars);

for j = 1:1:n_vars
    scaled_data(:, j) = unscaled_data(:, j) / (scalings(j) + a_tol);
end

end
