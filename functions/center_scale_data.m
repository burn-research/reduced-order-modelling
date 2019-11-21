function [varargout] = center_scale_data(data, scaling_criterion, varargin)
%% Description:
% Input: 
% - data, obs x var

%% Inputs
n_args = length(varargin);
if n_args > 0
    my_mean = varargin{1};
end
if n_args > 1
    my_gamma = varargin{1};
end

%% Main
if scaling_criterion == 0
    % No centering, no scaling
    centered_scaled_data = data;
    mean_column = 0; 
    scaling_factors = ones(size(data,1), 1);
elseif scaling_criterion == -1
    % Scale according to the scaling_criterion provided but use an
    % externally provided mean vector
    [centered_scaled_data, mean_column, scaling_factors] = ...
        scale_here(data', scaling_criterion, my_mean);
else
    % Center scale according to the scaling_criterion provided
    [centered_scaled_data, mean_column, scaling_factors] = ...
        scale_here(data', scaling_criterion);
end

%% Output
if nargout > 0
    varargout{1} = centered_scaled_data;
end
if nargout > 1
    varargout{2} = mean_column;
end
if nargout > 2
    varargout{3} = scaling_factors;
end

end

function [X, mean_value, gamma] = scale_here(X, scal_crit, my_mean, my_gamma)
% Different choices are available:  
% -1) Mean column provided externally
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
% 
% The input matrix is directly changed in order to avoid creating copies of
% it. You might run out of memory when working with huge matrices
%

% Definition of parameters
[rows, columns] = size(X);
a_tol = 1e-16;

% Get the mean and center the data
if nargin < 3 || isempty(my_mean)
    mean_value = mean(X, 1);
elseif nargin >= 3 && ~isempty(my_mean) && scal_crit == -1
    mean_value = my_mean;
end
X = X - repmat(mean_value, rows, 1);

kurt_value = kurtosis(X);
min_value = min(X, [], 1);
max_value = max(abs(X), [], 1);
range_value = max_value - min_value;
std_value = std(X, 1, 1);
pareto_value = std_value.^0.5;
vast_value = std_value.^2 ./ (mean_value + a_tol);
level_value = mean_value;
vast_2 = std_value.^2 .* kurt_value.^2 ./ (mean_value + a_tol);
vast_3 = std_value.^2 .* kurt_value.^2 ./ max_value;
vast_4 = std_value.^2 .* kurt_value.^2 ./ (range_value + a_tol);


if nargin < 4
    switch scal_crit
        case 0
            gamma = ones(1, columns);
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
    gamma = my_gamma;
end

% Scale
try
    X = X ./ (repmat(gamma(:), 1, columns) + a_tol);
catch
    for j = 1 : columns
        X(:, j) = X(:, j) / (gamma(j) + a_tol);
    end
end

% Must be columns
mean_value = mean_value(:);
gamma = gamma(:);
X = X';

end

