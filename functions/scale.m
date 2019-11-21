% function [scaled_data, gamma] = scale(unscaled_data, uncentered_data, scal_crit, my_gamma)
% 
% This function provides a tool to scale the data according to several
% criteria
%
% INPUTS:
% 
% unscaled_data     = Data matrix to be scaled. Each row is an observation,
%                     each column a variable
% scal_crit         = Scaling criterion
%
% OUTPUT
%
% scaled_data       = Matrix of scaled data
% gamma             = Scaling parameter

function [scaled_data, gamma] = scale(unscaled_data, uncentered_data, scal_crit, my_gamma) 
% 0) No scaling 
% 1) Auto-scaling (STD), each variable is normalized by its standard
% deviation 
% 2) RANGE each variable is normalized by its range 
% 3) PARETO, each variable is scaled by the square root of its standard
% deviation  
% 4) VAST, each variable is scaled by the standard deviation and
% coefficient of variation 
% 5) LEVEL, each variable is normalized by the mean of the data
% 6) MAX, each variable is scaled by its maximum value

% Definition of parameters
[rows, columns] = size(unscaled_data);
a_tol = 1e-16;
% Main values to estimate
kurt_value = kurtosis(uncentered_data);
mean_value = mean(uncentered_data, 1);
min_value = min(uncentered_data, [], 1);
max_value = max(abs(uncentered_data), [], 1);
range_value = max_value - min_value;
std_value = std(uncentered_data, 1, 1);
pareto_value = std_value.^0.5;
vast_value = std_value.^2 ./ (mean_value + a_tol);
level_value = mean_value;
vast_2 = std_value.^2 .* kurt_value.^2 ./ (mean_value + a_tol);
vast_3 = std_value.^2 .* kurt_value.^2 ./ max_value;
vast_4 = std_value.^2 .* kurt_value.^2 ./ range_value;
% Choose scaling factors
if ~exist('scal_crit', 'var') || isempty(scal_crit)
    scal_crit = 1;
end
if ~exist('my_gamma', 'var') || isempty(my_gamma)
    switch scal_crit
        case 0
            gamma = ones(1, columns);
            scaled_data = unscaled_data;
            return;
        case 1
            gamma = std_value;
        case 2
            gamma = range_value;
        case 3
            gamma = pareto_value;
        case 4
            gamma = vast_value;
        case 5
            gamma = level_value;
        case 6
            gamma = max_value;
        case 7
            gamma = vast_2;
        case 8
            gamma = vast_3;
        case 9
            gamma = vast_4;
        otherwise
            error('Unknown scaling criterion');
    end
else
    % Provided externally
    gamma = my_gamma;
end
% Initialization and scaling
scaled_data = zeros(rows, columns);
for j = 1 : columns
    scaled_data(:, j) = unscaled_data(:, j) / (gamma(j) + a_tol);
end
end
% save gamma.out gamma -ASCII -DOUBLE
%save scaled_data scaled_data